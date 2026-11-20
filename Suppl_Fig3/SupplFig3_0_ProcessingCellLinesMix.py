#!/usr/bin/env python3
"""Process the TrES cell-line mix and write QC, cell lists and merged profiles."""

import csv
import gzip
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib_venn import venn3
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
from upsetplot import from_indicators, plot as upset_plot

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get("DATA_ROOT", SCRIPT_DIR / "inputs")).expanduser().resolve()
WORK_DIR = Path(os.environ.get("WORK_DIR", SCRIPT_DIR / "results")).expanduser().resolve()
BAM_DIR = Path(os.environ.get("BAM_DIR", DATA_ROOT)).expanduser().resolve()
RNA_DIR = Path(os.environ.get("RNA_DIR", DATA_ROOT)).expanduser().resolve()
CELL_CYCLE_GENES = Path(os.environ.get("CELL_CYCLE_GENES", DATA_ROOT / "regev_lab_cell_cycle_genes.txt")).resolve()
DNA_BLACKLIST = Path(os.environ.get("DNA_BLACKLIST", DATA_ROOT / "hg38-blacklist.bed.gz")).resolve()
HUMAN_ANNOTATION = DATA_ROOT / "gencode_v41_GRCh38.gff3.gz"
MOUSE_ANNOTATION = DATA_ROOT / "gencode_vM30_GRCm39.gff3.gz"
MEX_FILES = ("matrix.mtx", "features.tsv", "barcodes.tsv")
MERGED_RNA = "AllCellLines_rna_merged_processed_Filtered.h5ad"
FRAGMENT_THREADS = 64

HUMAN_CELLS = ("GM12878", "JJN3", "Karpas422", "WSU")
MARKS = {"ac": "H3K27ac", "me3": "H3K27me3"}
RNA_CLUSTER_LABELS = {"0": "WSU", "1": "JJN3", "2": "Karpas422", "3": "GM12878"}
DNA_CLUSTER_LABELS = {"0": "WSU", "1": "Karpas422", "2": "JJN3", "3": "GM12878"}
RNA_COUNT_FIELDS = ("log1p_n_genes_by_counts", "log1p_total_counts")
RNA_QC_FIELDS = (*RNA_COUNT_FIELDS, "pct_counts_mt", "pct_counts_in_top_50_genes")
UMI_MIN, UMI_MAX, MAX_MT = 300, 6000, 10
# Minimum and maximum fragments per mark, then maximum percentage in top 50 RNA genes.
CELL_QC = {
    "GM12878": ((300, 6000), (100, 6000), 40),
    "JJN3": ((300, 6000), (300, 6000), 30),
    "Karpas422": ((300, 6000), (1000, 6000), 30),
    "WSU": ((100, 6000), (300, 6000), 30),
    "A20": ((1000, 12000), (400, 10000), 40),
}
DNA_FEATURES = {"ac": (5000, 40000), "me3": (40000, 30000)}


def relabel_ticks(transform, title=None):
    for axis in plt.gcf().axes:
        axis.set_yticklabels([f"{transform(tick):.0f}" for tick in axis.get_yticks()])
        if title is not None:
            axis.set_title(title)


def fragment_violin(adata, title, filename):
    sc.pl.violin(adata, "log10_n_fragment", rotation=0.0000001, show=False)
    relabel_ticks(lambda tick: 10**tick - 1, title)
    plt.savefig(f"./figures/{filename}")
    plt.close()


def import_dna(sample, genome, annotation):
    assays = [f"{sample}_{mark}" for mark in MARKS.values()]
    for assay in assays:
        stats = snap.pp.make_fragment_file(
            str(BAM_DIR / f"Sc_CMD_{assay}_NoDup.bam"),
            str(WORK_DIR / f"{assay}.tsv.gz"),
            is_paired=True, barcode_tag="CB", shift_left=4, shift_right=-5,
            min_mapq=20, chunk_size=500000000, chrM=["chrM", "M"], tempdir=str(WORK_DIR),
        )
        with open(f"{assay}.csv", "w", newline="") as handle:
            csv.writer(handle).writerows(stats.items())

    datasets = {
        assay: snap.pp.import_fragments(
            str(WORK_DIR / f"{assay}.tsv.gz"), genome, min_num_fragments=0,
            sorted_by_barcode=True, chrM=["chrM", "M"], shift_left=0, shift_right=0,
            chunk_size=50000, tempdir=str(WORK_DIR), backend="hdf5", n_jobs=FRAGMENT_THREADS,
        )
        for assay in assays
    }
    gene_anno = annotation if annotation.is_file() else genome
    for adata in datasets.values():
        snap.metrics.tsse(adata, gene_anno=gene_anno, inplace=True)
    for adata in datasets.values():
        snap.pp.filter_cells(adata, min_counts=100, max_counts=80000, min_tsse=0, inplace=True)
    for adata in datasets.values():
        snap.metrics.frag_size_distr(adata, add_key="frag_size_distr", inplace=True)

    for assay, adata in datasets.items():
        counts = adata.obs["n_fragment"]
        adata.obs["log10_n_fragment"] = np.log10(counts)
        fragment_violin(
            adata, f"n_cells: {adata.n_obs}\nMedian n_fragment:{np.median(counts)}",
            f"{assay}_FragViolin.pdf",
        )
    for assay, adata in datasets.items():
        sc.pl.scatter(adata, "log10_n_fragment", "tsse", title=assay, save=f"{assay}_tsse.pdf")
    for assay, adata in datasets.items():
        snap.pl.tsse(
            adata, min_fragment=0, width=750, height=600, interactive=False, show=False,
            out_file=f"./figures/{assay}_tsseDensity_PreFilter.pdf",
        )
        snap.pl.frag_size_distr(
            adata, width=750, height=600, interactive=False, show=False,
            out_file=f"./figures/{assay}_frag_size_distr_PreFilter.pdf",
        )
    for assay, adata in datasets.items():
        adata.obs["Exp"] = [assay] * adata.n_obs
        snap.ex.export_fragments(adata, groupby="Exp", prefix="", suffix=".tsv.gz")
        adata.write(filename=f"{assay}.h5ad")


def read_rna_matrix(source):
    if (source / "genes.tsv").is_file() or all((source / f"{name}.gz").is_file() for name in MEX_FILES):
        return sc.read_10x_mtx(source)

    # Give Scanpy the same gzipped layout when STARsolo files were downloaded plain.
    with TemporaryDirectory(prefix="rna_", dir=WORK_DIR) as scratch:
        directory = Path(scratch)
        for name in MEX_FILES:
            compressed = source / f"{name}.gz"
            target = directory / compressed.name
            if compressed.is_file():
                target.symlink_to(compressed)
            else:
                with (source / name).open("rb") as plain, gzip.open(target, "wb", compresslevel=1) as gzipped:
                    shutil.copyfileobj(plain, gzipped)
        return sc.read_10x_mtx(directory)


def import_rna(sample):
    prefix = f"CellLineMix_{sample}_RNA_"
    if all((RNA_DIR / f"{prefix}{name}.gz").is_file() for name in MEX_FILES):
        adata = sc.read_10x_mtx(RNA_DIR, prefix=prefix)
    else:
        adata = read_rna_matrix(RNA_DIR / f"{sample}.Solo.outGeneFull/filtered")
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^Hb[ab]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    sc.pl.violin(
        adata, ["pct_counts_mt", "pct_counts_in_top_50_genes"],
        jitter=0.4, multi_panel=True, rotation=0.0000001, save=f"{sample}_Pct_PreFilter.pdf",
    )
    sc.pl.violin(
        adata, list(RNA_COUNT_FIELDS),
        jitter=0.4, multi_panel=True, rotation=0.0000001, show=False,
    )
    relabel_ticks(np.expm1)
    plt.savefig(f"./figures/{sample}_PreFilter.pdf")
    plt.close()

    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
        title=f"{sample}\nColor: pct_counts_mt", save=f"{sample}_PreFilter.pdf",
    )
    sc.pl.scatter(
        adata, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",
        title=f"{sample}\nColor: pct_counts_mt", save=f"{sample}_PreFilterlog.pdf",
    )
    adata.write(f"{sample}.h5ad")


def filter_modalities(sample):
    ac_limits, me3_limits, max_top50 = CELL_QC[sample]
    rna = sc.read_h5ad(f"{sample}.h5ad")
    ac = sc.read_h5ad(f"{sample}_H3K27ac.h5ad")
    me3 = sc.read_h5ad(f"{sample}_H3K27me3.h5ad")
    obs = rna.obs
    keep = (
        (obs["total_counts"] >= UMI_MIN)
        & (obs["total_counts"] <= UMI_MAX)
        & (obs["pct_counts_mt"] < MAX_MT)
        & (obs["pct_counts_in_top_50_genes"] < max_top50)
    )
    rna[keep].write(f"{sample}_rna_CountFiltered.h5ad")
    for modality, adata, limits in (("ac", ac, ac_limits), ("me3", me3, me3_limits)):
        keep = snap.pp.filter_cells(
            adata, min_counts=limits[0], max_counts=limits[1], min_tsse=0, inplace=False,
        )
        adata[keep].write(f"{sample}_{modality}_CountFiltered.h5ad")


def write_cell_list(filename, adata):
    Path(filename).write_text("\n".join(adata.obs_names))


def write_rna_cluster_lists(adata):
    for sample in adata.obs["sampleMerge"].cat.categories:
        filename = f"{sample}_CountFiltered_cell_names.lst"
        selected = adata[adata.obs["sampleMerge"] == sample]
        write_cell_list(filename, selected)
        print(f"Wrote {selected.n_obs} cells to {filename}")


def embed_rna(adata, filtered=False):
    suffix = "_Filtered" if filtered else ""
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
    sc.pl.highly_variable_genes(adata, show=False, save=f"AllCellLines_HVG{suffix}.pdf")
    sc.tl.pca(adata, n_comps=30)
    sc.pl.pca_variance_ratio(
        adata, n_pcs=30, log=True, show=False,
        save=f"AllCellLines_pca_variance_ratio{suffix}.pdf",
    )
    if not filtered:
        colors = ["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"]
        for dims, label in (((0, 1), "1-2"), ((2, 3), "3-4")):
            sc.pl.pca(
                adata, color=colors, dimensions=[dims] * len(colors), ncols=3,
                show=False, save=f"AllCellLines_PCA{label}_QC.pdf",
            )
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, metric="cosine")
    sc.tl.umap(adata)


def merge_human_rna():
    datasets = []
    for sample in HUMAN_CELLS:
        adata = sc.read_h5ad(f"{sample}_rna_CountFiltered.h5ad")
        adata.obs["sample"] = sample
        datasets.append(adata)
    adata = ad.concat(datasets, join="inner", label="sampleMerge", keys=HUMAN_CELLS, index_unique=None)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    genes = [line.strip() for line in CELL_CYCLE_GENES.read_text().splitlines()]
    sc.tl.score_genes_cell_cycle(adata, s_genes=genes[:43], g2m_genes=genes[43:])
    embed_rna(adata)
    sc.pl.umap(
        adata, color=["sample", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo",
                      "pct_counts_in_top_50_genes", "log1p_total_counts"],
        ncols=3, show=False, save="AllCellLines_UMAP_QC.pdf",
    )
    sc.pl.umap(
        adata, color=["S_score", "G2M_score", "phase", "sample"],
        ncols=2, show=False, save="AllCellLines_UMAP_QC_CellCycle.pdf",
    )

    qc_colors = [*RNA_COUNT_FIELDS, "sample", "sampleMerge", "leiden"]
    sc.tl.leiden(adata, resolution=1)
    sc.pl.umap(
        adata, color=qc_colors, legend_loc="on data", ncols=3,
        show=False, save="AllCellLines_UMAP_RNA_QC.pdf",
    )
    adata = adata[adata.obs["leiden"] != "5"].copy()
    embed_rna(adata, filtered=True)
    sc.tl.leiden(adata, resolution=0.5)
    adata.obs["sampleMerge"] = adata.obs["leiden"].map(RNA_CLUSTER_LABELS).astype("category")
    sc.pl.umap(
        adata, color=qc_colors, legend_loc="on data", ncols=3,
        show=False, save="AllCellLines_UMAP_RNA_QC_Filtered.pdf",
    )
    return adata


def rank_rna_markers(adata):
    adata.obs["CL"] = adata.obs["sampleMerge"]
    sc.tl.rank_genes_groups(adata, groupby="CL", method="wilcoxon", key_added="CL_rgg")
    sc.tl.dendrogram(adata, groupby="CL")
    sc.pl.rank_genes_groups_dotplot(
        adata, groupby="CL", key="CL_rgg", standard_scale="var", n_genes=5,
        show=False, save="Cell_Lines_DGE_CL.pdf",
    )
    adata.layers["log1p_Norm_counts"] = adata.X.copy()


def plot_qc(sample, rna, ac, me3, shared=False):
    suffix = "Shared" if shared else ""
    counts_name = "Number of shared cells" if shared else "Number of cells"
    mt_label = "pct_counts_mt (shared)" if shared else "pct_counts_mt"

    sc.pl.violin(
        rna, list(RNA_QC_FIELDS),
        jitter=0.4, multi_panel=True, rotation=0.0000001, save=f"{sample}_rna_PostFilter{suffix}.png",
    )
    sc.pl.violin(
        rna, list(RNA_COUNT_FIELDS), jitter=0.4, multi_panel=True,
        rotation=0.0000001, show=False,
    )
    relabel_ticks(np.expm1)
    plt.title(f"{counts_name}: {rna.n_obs}\nMedian UMIs per cell: {np.median(rna.obs['total_counts'])}")
    plt.savefig(f"./figures/{sample}_rna_PostFilter{suffix}NoLog.pdf")
    plt.close()

    for prefix, log in (("", ""), ("log1p_", "log")):
        sc.pl.scatter(
            rna, f"{prefix}total_counts", f"{prefix}n_genes_by_counts", color="pct_counts_mt",
            title=f"{sample}\nColor: {mt_label}", save=f"{sample}_rna_PostFilter{suffix}{log}.png",
        )

    cell_label = "n_shared_cells" if shared else "n_cells"
    for modality, adata in (("me3", me3), ("ac", ac)):
        counts = adata.obs["n_fragment"]
        median = np.median(counts) if modality == "me3" else counts.median()
        fragment_violin(
            adata, f"{cell_label}: {adata.n_obs}\nMedian n_fragment: {median}",
            f"{sample}_{modality}_nfragPostFilter{suffix}.png",
        )


def plot_overlaps(sample, barcodes):
    ac_limits, me3_limits, max_top50 = CELL_QC[sample]
    rna, ac, me3 = barcodes
    plt.figure(figsize=(8, 8))
    venn3(
        barcodes,
        set_labels=(
            f"RNA\n{UMI_MIN}<UMIs<{UMI_MAX}\n<{MAX_MT}% MT\n<{max_top50}% 50TopGenes",
            f"H3K27ac\n{ac_limits[0]}<fragments<{ac_limits[1]}",
            f"H3K27me3\n{me3_limits[0]}<fragments<{me3_limits[1]}",
        ),
    )
    plt.suptitle(
        f"Venn Diagram of Passing Cells\nfor RNA, H3K27ac, and H3K27me3\n{sample}",
        fontsize=16, fontweight="bold",
    )
    plt.tight_layout()
    plt.savefig(f"./figures/{sample}_VennPassingCells.pdf")
    plt.close()

    all_cells = rna.union(ac).union(me3)
    frame = pd.DataFrame(
        {mark: [cell in selected for cell in all_cells]
         for mark, selected in zip(("RNA", "H3K27ac", "H3K27me3"), barcodes)},
        index=list(all_cells),
    )
    fig = plt.figure(figsize=(8, 10))
    upset_plot(
        from_indicators(frame.columns, data=frame), sort_by="cardinality", show_percentages=True,
        min_subset_size="0.2%", facecolor="darkblue", fig=fig, element_size=None,
    )
    plt.title(f"Passing Cells\nfor RNA, H3K27ac, and H3K27me3\n{sample}", fontsize=16, fontweight="bold")
    plt.savefig(f"./figures/{sample}_UpsetPassingCells.pdf")
    plt.close()


def export_shared_cells(sample):
    rna, ac, me3 = (
        sc.read_h5ad(f"{sample}_{modality}_CountFiltered.h5ad")
        for modality in ("rna", "ac", "me3")
    )
    plot_qc(sample, rna, ac, me3)
    write_cell_list(f"{sample}_CountFiltered_cell_names.lst", rna)
    barcodes = [set(adata.obs_names) for adata in (rna, ac, me3)]
    shared = set.intersection(*barcodes)
    plot_overlaps(sample, barcodes)
    if not shared:
        print(f"No shared barcodes for {sample}; skipping shared-cell plots and files.")
        return

    rna, ac, me3 = (adata[adata.obs_names.isin(shared)] for adata in (rna, ac, me3))
    write_cell_list(f"{sample}_CountFilteredShared_cell_names.lst", rna)
    plot_qc(sample, rna, ac, me3, shared=True)
    for modality, adata in (("rna", rna), ("ac", ac), ("me3", me3)):
        adata.write(f"{sample}_{modality}_CountFilteredShared.h5ad")
    for adata in (ac, me3):
        snap.ex.export_fragments(
            adata, groupby="Exp", prefix="", suffix="_FiltShared.bed.gz", compression_level=1,
        )

def read_shared_dna(sample, modality):
    adata = snap.read(f"{sample}_{modality}_CountFilteredShared.h5ad", backed=None)
    adata.obs.index = adata.obs.index.astype(str) + "_" + sample
    return adata


def merge_human_dna():
    for modality, mark in MARKS.items():
        experiment = f"Cell_Lines_{modality}"
        bin_size, n_features = DNA_FEATURES[modality]
        sources = []
        for sample in HUMAN_CELLS:
            adata = read_shared_dna(sample, modality)
            adata.obs["cell_line"] = [sample] * adata.n_obs
            filename = f"{sample}_qcTOREMOVE.h5ad"
            adata.write(filename)
            sources.append((sample, filename))

        dataset = snap.AnnDataSet(sources, filename=f"{experiment}.h5ad", add_key="sample")
        dataset.obs["Exp"] = [experiment] * dataset.n_obs
        snap.ex.export_fragments(dataset, groupby="Exp", prefix="", suffix=".tsv.gz", compression_level=1)
        dataset.close()

        samples = []
        for sample in HUMAN_CELLS:
            adata = read_shared_dna(sample, modality)
            del adata.obsm
            samples.append(adata)
        adata = ad.concat(samples, join="inner", label="sample", keys=HUMAN_CELLS, index_unique=None)
        fragments = snap.pp.import_fragments(
            fragment_file=f"{experiment}.tsv.gz", chrom_sizes=snap.genome.GRCh38,
            sorted_by_barcode=True, min_num_fragments=0, tempdir=".",
        )
        common = adata.obs.index.intersection(fragments.obs.index)
        print(f"{mark}: {len(common)}/{adata.n_obs} cells matched to imported fragments")
        adata = adata[common].copy()
        fragments = fragments[common].copy()
        adata.obsm = fragments.obsm.copy()
        adata.uns["reference_sequences"] = fragments.uns["reference_sequences"].copy()
        del fragments

        snap.pp.add_tile_matrix(
            adata, bin_size=bin_size, counting_strategy="paired-insertion",
            chunk_size=500000, inplace=True,
        )
        snap.pp.select_features(adata, n_features=n_features, inplace=True, blacklist=str(DNA_BLACKLIST))
        snap.tl.spectral(
            adata, n_comps=8, weighted_by_sd=False, chunk_size=80000,
            features="selected", distance_metric="cosine", inplace=True,
        )
        snap.tl.umap(adata, use_rep="X_spectral", key_added="umap", random_state=None)
        snap.pp.knn(adata, n_neighbors=20, use_rep="X_spectral", method="kdtree")
        snap.pl.spectral_eigenvalues(
            adata, width=600, height=400, show=False, interactive=False,
            out_file=f"./figures/AllCellLines_{mark}_Eigenvalues.png",
        )
        sc.pl.umap(
            adata, color=["n_fragment", "log10_n_fragment", "tsse", "sample"], ncols=2,
            save=f"AllCellLines_{mark}_UMAP_QC.pdf", show=False,
        )
        snap.tl.leiden(adata, resolution=0.05)
        sc.pl.umap(
            adata, color=["leiden", "sample"], legend_loc="on data", ncols=2,
            save=f"AllCellLines_{mark}_UMAP_Leiden.pdf", show=False,
        )
        adata.write(f"AllCellLines_{mark}_merged_processed.h5ad")

    for mark in MARKS.values():
        filename = f"AllCellLines_{mark}_merged_processed.h5ad"
        adata = sc.read_h5ad(filename)
        adata.obs["sampleMerge"] = adata.obs["leiden"].map(DNA_CLUSTER_LABELS).astype("category")
        adata.write(filename)


def plot_linked_umaps(
    rna_file, ac_file, me3_file, output_prefix="./figures/LinkedUMAPs_RNA_Ac_Me3_lifted",
):
    rna, ac, me3 = (sc.read_h5ad(path) for path in (rna_file, ac_file, me3_file))
    for adata in (ac, me3):
        adata.obs_names = adata.obs_names.str.rsplit("_", n=1).str[0]

    common = rna.obs_names.intersection(ac.obs_names).intersection(me3.obs_names)
    if not len(common):
        raise ValueError("No shared cell barcodes across RNA, H3K27ac and H3K27me3")
    rna, ac, me3 = (adata[common] for adata in (rna, ac, me3))
    umaps = []
    for i, adata in enumerate((rna, ac, me3), start=1):
        umaps.append(pd.DataFrame(adata.obsm["X_umap"], index=common, columns=[f"x{i}", f"y{i}"]))
    rna_umap, ac_umap, me3_umap = umaps

    gap = 0.20 * max(float(coords.iloc[:, 0].max() - coords.iloc[:, 0].min()) for coords in umaps)
    ac_umap["x2"] += (rna_umap["x1"].max() - ac_umap["x2"].min()) + gap
    me3_umap["x3"] += (ac_umap["x2"].max() - me3_umap["x3"].min()) + gap
    lift = 0.20 * max(float(coords.iloc[:, 1].max() - coords.iloc[:, 1].min()) for coords in umaps)
    rna_umap["y1"] += lift
    me3_umap["y3"] += lift

    coords = pd.concat(umaps, axis=1)
    coords["cell_type"] = rna.obs["sampleMerge"].astype("category")
    colors = {"GM12878": "#EF9649", "JJN3": "#67AB58", "Karpas422": "#C9544F", "WSU": "#9B7EBE"}
    link_colors = {cell: mcolors.to_rgba(color, alpha=0.22) for cell, color in colors.items()}
    fallback = mcolors.to_rgba("#BDBDBD", alpha=0.22)
    line_colors = np.vstack([link_colors.get(cell, fallback) for cell in coords["cell_type"].astype(object)])

    fig, axis = plt.subplots(figsize=(18, 6), dpi=300)
    axis.set_title("RNA → H3K27ac → H3K27me3", fontsize=20, weight="bold")
    for left, right in ((1, 2), (2, 3)):
        segments = np.stack([coords[[f"x{i}", f"y{i}"]].values for i in (left, right)], axis=1)
        axis.add_collection(LineCollection(segments, colors=line_colors, linewidths=0.03, zorder=0))

    for i in (1, 2, 3):
        for cell, color in colors.items():
            rows = coords[coords["cell_type"] == cell]
            if rows.empty:
                continue
            axis.scatter(
                rows[f"x{i}"], rows[f"y{i}"], s=16, c=color, edgecolors="none",
                alpha=0.75, label=cell if i == 1 else None, zorder=1, rasterized=False,
            )
    axis.set_aspect("equal", adjustable="datalim")
    axis.margins(0.03)
    axis.axis("off")
    handles, labels = axis.get_legend_handles_labels()
    entries = dict(zip(labels, handles))
    legend = axis.legend(
        entries.values(), entries.keys(), title="Cell Type", bbox_to_anchor=(1.01, 1),
        loc="upper left", markerscale=1.6, fontsize=9,
    )
    legend.get_title().set_fontsize(10)
    plt.tight_layout()
    for extension in ("pdf", "png"):
        path = f"{output_prefix}.{extension}"
        plt.savefig(path, bbox_inches="tight")
        print(path)
    plt.close(fig)


def process_samples(samples, genome, annotation):
    for sample in samples:
        import_dna(sample, genome, annotation)
    for sample in samples:
        import_rna(sample)
    for sample in samples:
        filter_modalities(sample)


def plot_rna_singlets(adata):
    scores = ["scDblFinder_score", "scDblFinder_class"]
    if not all(column in adata.obs for column in scores):
        print("Skipping AllCellLines_UMAP_RNA_singlet_filtered.pdf: no scDblFinder columns.")
        return
    singlets = adata[adata.obs["leiden"] != "6"].copy()
    singlets.obs["CL"] = singlets.obs["sample"].str.split("_").str[-1]
    sc.pl.umap(
        singlets, color=[*RNA_COUNT_FIELDS, "sample", "sampleMerge", "leiden", *scores],
        legend_loc="on data", ncols=3, show=False, save="AllCellLines_UMAP_RNA_singlet_filtered.pdf",
    )


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    os.chdir(WORK_DIR)
    Path("figures").mkdir(exist_ok=True)
    process_samples(HUMAN_CELLS, snap.genome.GRCh38, HUMAN_ANNOTATION)

    rna = merge_human_rna()
    write_rna_cluster_lists(rna)
    rank_rna_markers(rna)
    rna.write(MERGED_RNA)
    for sample in HUMAN_CELLS:
        export_shared_cells(sample)

    # Keep the original second RNA QC pass and final count lists.
    del rna
    rna = merge_human_rna()
    write_rna_cluster_lists(rna)
    plot_rna_singlets(rna)
    del rna

    merge_human_dna()
    plot_linked_umaps(
        MERGED_RNA, "AllCellLines_H3K27ac_merged_processed.h5ad",
        "AllCellLines_H3K27me3_merged_processed.h5ad",
    )

    process_samples(("A20",), snap.genome.GRCm39, MOUSE_ANNOTATION)
    export_shared_cells("A20")


if __name__ == "__main__":
    main()

