#################################################
#################################################
#########################################
#################################################
#################################################
#########################################
#################################################
#################################################
#################### CL HUMAN #####################
#################################################
#################################################
#########################################
#################################################
#################################################
#################################################
#########################################
#################################################

import matplotlib.pyplot as plt
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
from upsetplot import UpSet, plot, from_indicators
import os
import csv

if not os.path.exists('figures'):
    os.makedirs('figures')

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/CL/"
srcdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"
Experiment='H2R_CellLines'

#################################################
#################### DNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS DNA')
maxFrags = 80000

cell_lines = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

for sample in cell_lines:
    # H3K27ac
    bam_ac = f"Sc_CMD_{sample}_H3K27ac_NoDup.bam"
    exp_ac = sample + '_H3K27ac'
    bam_ac_fgt = exp_ac + '.tsv.gz'
    bam_ac_ad = exp_ac + '.h5ad'

    # H3K27me3
    bam_me3 = f"Sc_CMD_{sample}_H3K27me3_NoDup.bam"
    exp_me3 = sample + '_H3K27me3'
    bam_me3_fgt = exp_me3 + '.tsv.gz'
    bam_me3_ad = exp_me3 + '.h5ad'


    # Make fragment files
    ac_info = snap.pp.make_fragment_file(
        srcdir + bam_ac, workdir + bam_ac_fgt, is_paired=True, barcode_tag='CB',
        shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
        chrM=['chrM', 'M'], tempdir=workdir
    )
    with open(f'{exp_ac}.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        for key, value in ac_info.items():
            writer.writerow([key, value])

    me3_info = snap.pp.make_fragment_file(
        srcdir + bam_me3, workdir + bam_me3_fgt, is_paired=True, barcode_tag='CB',
        shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
        chrM=['chrM', 'M'], tempdir=workdir
    )
    with open(f'{exp_me3}.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        for key, value in me3_info.items():
            writer.writerow([key, value])

    # Import fragments
    ac_ad = snap.pp.import_fragments(
        workdir + bam_ac_fgt, snap.genome.GRCh38, min_num_fragments=0, sorted_by_barcode=True,
        chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64
    )
    me3_ad = snap.pp.import_fragments(
        workdir + bam_me3_fgt, snap.genome.GRCh38, min_num_fragments=0, sorted_by_barcode=True,
        chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64
    )

    # Compute metrics
    snap.metrics.tsse(ac_ad, gene_anno=snap.genome.GRCh38, inplace=True)
    snap.metrics.tsse(me3_ad, gene_anno=snap.genome.GRCh38, inplace=True)

    snap.pp.filter_cells(ac_ad, min_counts=100, max_counts=maxFrags, min_tsse=0, inplace=True)
    snap.pp.filter_cells(me3_ad, min_counts=100, max_counts=maxFrags, min_tsse=0, inplace=True)

    snap.metrics.frag_size_distr(ac_ad, add_key='frag_size_distr', inplace=True)
    snap.metrics.frag_size_distr(me3_ad, add_key='frag_size_distr', inplace=True)

    ac_ad.obs['log10_n_fragment'] = np.log10(ac_ad.obs['n_fragment'])
    me3_ad.obs['log10_n_fragment'] = np.log10(me3_ad.obs['n_fragment'])

    # Violin plot for ac
    sc.pl.violin(ac_ad, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {ac_ad.n_obs}\nMedian n_fragment:{np.median(ac_ad.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp_ac}_FragViolin.pdf')
    plt.close()

    # Violin plot for me3
    sc.pl.violin(me3_ad, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {me3_ad.n_obs}\nMedian n_fragment:{np.median(me3_ad.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp_me3}_FragViolin.pdf')
    plt.close()

    # Scatter plots
    sc.pl.scatter(ac_ad, "log10_n_fragment", "tsse", title=exp_ac, save=exp_ac + '_tsse.pdf')
    sc.pl.scatter(me3_ad, "log10_n_fragment", "tsse", title=exp_me3, save=exp_me3 + '_tsse.pdf')

    # SnapATAC2 QC plots
    snap.pl.tsse(ac_ad, min_fragment=0, width=750, height=600, interactive=False, show=False,
                 out_file=f"./figures/{exp_ac}_tsseDensity_PreFilter.pdf")
    snap.pl.frag_size_distr(ac_ad, width=750, height=600, interactive=False, show=False,
                            out_file=f"./figures/{exp_ac}_frag_size_distr_PreFilter.pdf")

    snap.pl.tsse(me3_ad, min_fragment=0, width=750, height=600, interactive=False, show=False,
                 out_file=f"./figures/{exp_me3}_tsseDensity_PreFilter.pdf")
    snap.pl.frag_size_distr(me3_ad, width=750, height=600, interactive=False, show=False,
                            out_file=f"./figures/{exp_me3}_frag_size_distr_PreFilter.pdf")

    # Add experiment name to obs
    ac_ad.obs['Exp'] = [exp_ac] * ac_ad.n_obs
    me3_ad.obs['Exp'] = [exp_me3] * me3_ad.n_obs

    # Export fragments and AnnData
    snap.ex.export_fragments(ac_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
    ac_ad.write(filename=f'{exp_ac}.h5ad')
    snap.ex.export_fragments(me3_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
    me3_ad.write(filename=f'{exp_me3}.h5ad')


#################################################
#################### RNA QC #######
#################################################

print(f'IMPORT AND QC PLOTS RNA FOR 8 CELL LINES')


cell_line_samples = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

for sample in cell_line_samples:
    print(f"Processing {sample}")

    #adata = sc.read_csv(workdir + cnt_tsv, delimiter='\t', first_column_names=True)
    #adata = adata.T
    mtxPath = f"/mnt/dataFast/ahrmad/triseq_2025052/RNA_RE/{sample}.Solo.outGeneFull/filtered/"
    adata = sc.read_10x_mtx(mtxPath)

    # mitochondrial genes, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    # ribosomal genes (mouse: "Rps", "Rpl")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes (mouse: "Hba", "Hbb")
    adata.var["hb"] = adata.var_names.str.contains("^Hb[ab]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    # Violin plot for percent mt and top 50 genes
    sc.pl.violin(
        adata,
        ["pct_counts_mt", 'pct_counts_in_top_50_genes'],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        save=f'{sample}_Pct_PreFilter.pdf'
    )

    # Violin plot for log1p n_genes and total_counts
    sc.pl.violin(
        adata,
        ["log1p_n_genes_by_counts", "log1p_total_counts"],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        show=False  # Prevent auto-display so we can edit it
    )

    # Get the figure and axes
    fig = plt.gcf()
    axes = fig.axes

    # Loop through the axes to update y-axis tick labels
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]  # reverse log1p
        ax.set_yticklabels(new_labels)

    # Save or show the updated figure
    plt.savefig(f"./figures/{sample}_PreFilter.pdf")
    plt.close()

    # Scatter plots
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts",
        color="pct_counts_mt",
        title=sample + '\nColor: pct_counts_mt',
        save=f'{sample}_PreFilter.pdf'
    )
    sc.pl.scatter(
        adata, "log1p_total_counts", "log1p_n_genes_by_counts",
        color="pct_counts_mt",
        title=sample + '\nColor: pct_counts_mt',
        save=f'{sample}_PreFilterlog.pdf'
    )

    adata.write(f'{sample}.h5ad')


#################################################
############# RNA + DNA FILTERING ###############
#################################################

cell_lines = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

for sample in cell_lines:
    print(f"\n=== Filtering and QC for {sample} ===")
    exp = sample

    # Filtering parameters
    maxFrags_ac = 6000
    minFrags_ac = 300

    maxFrags_me3 = 6000
    minFrags_me3 = 300

    maxUMI = 6000
    minUMI = 300
    max_pct_50g = 30
    max_pct_mt = 10

    # Set minFrags and minUMI per sample
    if "Karpas422" in sample:
        minFrags_me3 = 1000
    if "GM12878" in sample:
        minFrags_me3 = 100
        max_pct_50g = 40
    if "WSU" in sample:
        minFrags_ac = 100

    # File names
    exp_rna = f"{sample}"
    exp_ac = f"{sample}_H3K27ac"
    exp_me3 = f"{sample}_H3K27me3"

    # Load data
    try:
        rna = sc.read_h5ad(f"{exp_rna}.h5ad")
    except Exception as e:
        print(f"Could not load {exp_rna}.h5ad: {e}")
        continue
    try:
        ac = sc.read_h5ad(f"{exp_ac}.h5ad")
    except Exception as e:
        print(f"Could not load {exp_ac}.h5ad: {e}")
        continue
    try:
        me3 = sc.read_h5ad(f"{exp_me3}.h5ad")
    except Exception as e:
        print(f"Could not load {exp_me3}.h5ad: {e}")
        continue

    # --- RNA filtering criteria ---
    filtered_rna = rna[rna.obs['total_counts'] >= minUMI]
    filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
    filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
    filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]

    # --- DNA filtering criteria ---
    filtered_ac = ac[snap.pp.filter_cells(ac, min_counts=minFrags_ac, max_counts=maxFrags_ac, min_tsse=0, inplace=False)]
    filtered_me3 = me3[snap.pp.filter_cells(me3, min_counts=minFrags_me3, max_counts=maxFrags_me3, min_tsse=0, inplace=False)]

    # Save filtered AnnData objects for each modality
    filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')
    filtered_ac.write(f'{exp}_ac_CountFiltered.h5ad')
    filtered_me3.write(f'{exp}_me3_CountFiltered.h5ad')



#################################################
############# RNA ANALYSIS ###############
#################################################


cell_lines = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

# List to hold AnnData objects for each sample
adatas = []
sample_names = []

for sample in cell_lines:
    exp = sample
    try:
        adata_sample = sc.read_h5ad(f"{exp}_rna_CountFiltered.h5ad")
        adata_sample.obs['sample'] = exp  # Add sample name as a column for later coloring
        adatas.append(adata_sample)
        sample_names.append(exp)
    except Exception as e:
        print(f"Could not load {exp}_rna_CountFiltered.h5ad: {e}")

# Concatenate all AnnData objects
if len(adatas) == 0:
    raise ValueError("No RNA AnnData objects loaded for cell lines.")

adata = ad.concat(
    adatas,
    join="inner",
    label="sampleMerge",
    keys=sample_names,
    index_unique=None  # Keep original cell barcodes
)

# Store raw counts in a layer
adata.layers["counts"] = adata.X.copy()


# Normalize and log1p
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Cell cycle scoring
cell_cycle_genes = [x.strip() for x in open('/mnt/dataFast/ahrmad/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Highly variable genes (across all cells)
#sc.pp.highly_variable_genes(adata, batch_key="sample")
sc.pp.highly_variable_genes(
    adata,
    flavor='seurat_v3',
    n_top_genes=2000,
)
sc.pl.highly_variable_genes(adata, show=False, save=f'AllCellLines_HVG.pdf')


# PCA
sc.tl.pca(adata,n_comps=30)
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True, show=False, save=f'AllCellLines_pca_variance_ratio.pdf')

# PCA plots
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"],
    dimensions=[(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)],
    ncols=3,
    show=False, save=f'AllCellLines_PCA1-2_QC.pdf'
)
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"],
    dimensions=[(2, 3), (2, 3), (2, 3), (2, 3), (2, 3)],
    ncols=3,
    show=False, save=f'AllCellLines_PCA3-4_QC.pdf'
)

# Neighbors/UMAP
#sc.pp.neighbors(adata, n_neighbors=30)
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, metric='cosine')  # or adjust as needed
sc.tl.umap(adata)

# UMAPs colored by QC and sample
sc.pl.umap(
    adata,
    color=["sample", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_50_genes", "log1p_total_counts"],
    ncols=3,
    show=False, save=f'AllCellLines_UMAP_QC.pdf'
)

sc.pl.umap(
    adata,
    color=['S_score', 'G2M_score', 'phase', 'sample'],
    ncols=2,
    show=False, save=f'AllCellLines_UMAP_QC_CellCycle.pdf'
)

sc.tl.leiden(adata, resolution=1)

sc.pl.umap(
    adata,
    color=["log1p_n_genes_by_counts", "log1p_total_counts","sample","sampleMerge","leiden"],
    legend_loc="on data",
    ncols=3,
    save=f'AllCellLines_UMAP_RNA_QC.pdf',
    show=False
)

adata = adata[adata.obs["leiden"] != "5"].copy()

sc.pp.highly_variable_genes(
    adata,
    flavor='seurat_v3',
    n_top_genes=2000,
)
sc.pl.highly_variable_genes(adata, show=False, save=f'AllCellLines_HVG_Filtered.pdf')

sc.tl.pca(adata,n_comps=30)
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True, show=False, save=f'AllCellLines_pca_variance_ratio_Filtered.pdf')

sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, metric='cosine')  # or adjust as needed
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5)

cluster2sample = {
    "0": "WSU",
    "1": "JJN2",
    "2": "Karpas422",
    "3": "GM12878",
}

# overwrite sampleMerge from leiden
adata.obs["sampleMerge"] = adata.obs["leiden"].map(cluster2sample).astype("category")


sc.pl.umap(
    adata,
    color=["log1p_n_genes_by_counts", "log1p_total_counts","sample","sampleMerge","leiden"],
    legend_loc="on data",
    ncols=3,
    save=f'AllCellLines_UMAP_RNA_QC_Filtered.pdf',
    show=False
)

# Loop over each cell line in sampleMerge
for exp in adata.obs["sampleMerge"].cat.categories:  # or: for exp in adata.obs["sampleMerge"].unique():
    # subset to that sample
    filtered_rna = adata[adata.obs["sampleMerge"] == exp]

    # get cell names (barcodes)
    cell_names = filtered_rna.obs_names

    # write to file: GM12878_CountFiltered_cell_names.lst, etc.
    out_fn = f"{exp}_CountFiltered_cell_names.lst"
    with open(out_fn, "w") as f:
        f.write("\n".join(cell_names))

    print(f"Wrote {len(cell_names)} cells to {out_fn}")

adata.obs["CL"] = adata.obs["sampleMerge"]

sc.tl.rank_genes_groups(adata, groupby="CL", method="wilcoxon",key_added="CL_rgg")
sc.tl.dendrogram(adata, groupby="CL")
sc.pl.rank_genes_groups_dotplot(adata, groupby="CL",key="CL_rgg", standard_scale="var", n_genes=5,show=False,save=f'Cell_Lines_DGE_CL.pdf')
adata.layers["log1p_Norm_counts"] = adata.X.copy()

adata.write("AllCellLines_rna_merged_processed_Filtered.h5ad")

############################MERGE
############################MERGE
############################MERGE
############################MERGE
############################MERGE
############################MERGE
############################MERGE

cell_lines = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

# List to hold AnnData objects for each sample
adatas = []
sample_names = []

for sample in cell_lines:
    exp = sample

    # Filtering parameters
    maxFrags_ac = 6000
    minFrags_ac = 300

    maxFrags_me3 = 6000
    minFrags_me3 = 300

    maxUMI = 6000
    minUMI = 300
    max_pct_50g = 30
    max_pct_mt = 10

    # Set minFrags and minUMI per sample
    if "Karpas422" in sample:
        minFrags_me3 = 1000
    if "GM12878" in sample:
        minFrags_me3 = 100
        max_pct_50g = 40
    if "WSU" in sample:
        minFrags_ac = 100

    # File names
    exp_rna = f"{sample}"
    exp_ac = f"{sample}_H3K27ac"
    exp_me3 = f"{sample}_H3K27me3"

    try:
        filtered_rna = sc.read_h5ad(f"{exp}_rna_CountFiltered.h5ad")
        lst_path = f"{exp}_CountFiltered_cell_names.lst"
        with open(lst_path) as f:
            keep_cells = [line.strip() for line in f if line.strip()]
        keep_cells = [c for c in keep_cells if c in filtered_rna.obs_names]
        print(f"{len(keep_cells)/len(filtered_rna.obs_names)} cells found in both list and h5ad.")
        adata_filt = filtered_rna[filtered_rna.obs_names.isin(keep_cells)].copy()

        filtered_ac = sc.read_h5ad(f'{exp}_ac_CountFiltered.h5ad')
        filtered_me3 = sc.read_h5ad(f'{exp}_me3_CountFiltered.h5ad')
    except Exception as e:
        print(f"Could not load {exp}_??_CountFiltered.h5ad: {e}")

    # Output plots for each modality separately

    # RNA plots
    sc.pl.violin(
        filtered_rna,
        ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", 'pct_counts_in_top_50_genes'],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        save=f'{exp}_rna_PostFilter.png'
    )
    sc.pl.violin(
        filtered_rna,
        ["log1p_n_genes_by_counts", "log1p_total_counts"],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        show=False
    )
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
    plt.title(f"Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.obs['total_counts'])}")
    plt.savefig(f"./figures/{exp}_rna_PostFilterNoLog.pdf")
    plt.close()

    sc.pl.scatter(filtered_rna, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt', save=f'{exp}_rna_PostFilter.png')
    sc.pl.scatter(filtered_rna, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt', save=f'{exp}_rna_PostFilterlog.png')

    # DNA plots
    sc.pl.violin(filtered_me3, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {filtered_me3.n_obs}\nMedian n_fragment: {np.median(filtered_me3.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_me3_nfragPostFilter.png')
    plt.close()

    sc.pl.violin(filtered_ac, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {filtered_ac.n_obs}\nMedian n_fragment: {filtered_ac.obs['n_fragment'].median()}")
    plt.savefig(f'./figures/{exp}_ac_nfragPostFilter.png')
    plt.close()

    cell_names = filtered_rna.obs_names
    with open(f"{exp}_CountFiltered_cell_names.lst", "w") as f:
        f.write("\n".join(cell_names))

    # Now try to merge the modalities based on shared barcodes
    rna_barcodes = set(filtered_rna.obs.index)
    ac_barcodes = set(filtered_ac.obs.index)
    me3_barcodes = set(filtered_me3.obs.index)
    shared_barcodes = rna_barcodes & ac_barcodes & me3_barcodes

    print(f'{sample}:\nAc:{len(ac_barcodes)}\nMe3:{len(me3_barcodes)}\nRNA:{len(rna_barcodes)}\nShared:{len(shared_barcodes)}')

    # Create a Venn diagram with the sets
    plt.figure(figsize=(8, 8))
    venn = venn3([rna_barcodes, ac_barcodes, me3_barcodes],
                 set_labels=(f'RNA\n{minUMI}<UMIs<{maxUMI}\n<{max_pct_mt}% MT\n<{max_pct_50g}% 50TopGenes',
                             f'H3K27ac\n{minFrags_ac}<fragments<{maxFrags_ac}',
                             f'H3K27me3\n{minFrags_me3}<fragments<{maxFrags_me3}'))
    plt.suptitle(f'Venn Diagram of Passing Cells\nfor RNA, H3K27ac, and H3K27me3\n{sample}',
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'./figures/{exp}_VennPassingCells.pdf')
    plt.close()

    # Create a DataFrame to represent the sets
    all_items = rna_barcodes.union(ac_barcodes).union(me3_barcodes)
    data = {
        'RNA': [item in rna_barcodes for item in all_items],
        'H3K27ac': [item in ac_barcodes for item in all_items],
        'H3K27me3': [item in me3_barcodes for item in all_items],
    }
    df = pd.DataFrame(data, index=list(all_items))

    fig = plt.figure(figsize=(8, 10))
    plot(
        from_indicators(df.columns, data=df),
        sort_by='cardinality',
        show_percentages=True,
        min_subset_size="0.2%",
        facecolor="darkblue",
        fig=fig, element_size=None
    )
    plt.title(f"Passing Cells\nfor RNA, H3K27ac, and H3K27me3\n{sample}", fontsize=16, fontweight='bold')
    plt.savefig(f'./figures/{exp}_UpsetPassingCells.pdf')
    plt.close()

    if len(shared_barcodes) == 0:
        print(f"No shared barcodes for {sample}, skipping shared-barcode analysis.")
        continue

    # Subset the filtered datasets based on shared barcodes
    filtered_ac_shared = filtered_ac[filtered_ac.obs.index.isin(shared_barcodes)]
    filtered_me3_shared = filtered_me3[filtered_me3.obs.index.isin(shared_barcodes)]
    filtered_rna_shared = filtered_rna[filtered_rna.obs.index.isin(shared_barcodes)]

    print(f"Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.obs['total_counts'])}")

    cell_names = filtered_rna_shared.obs_names
    with open(f"{exp}_CountFilteredShared_cell_names.lst", "w") as f:
        f.write("\n".join(cell_names))

    # Output plots for each modality, for the shared barcodes

    # RNA plots (shared)
    sc.pl.violin(
        filtered_rna_shared,
        ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", 'pct_counts_in_top_50_genes'],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        save=f'{exp}_rna_PostFilterShared.png'
    )
    sc.pl.violin(
        filtered_rna_shared,
        ["log1p_n_genes_by_counts", "log1p_total_counts"],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        show=False
    )
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
    plt.title(f"Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.obs['total_counts'])}")
    plt.savefig(f"./figures/{exp}_rna_PostFilterSharedNoLog.pdf")
    plt.close()

    sc.pl.scatter(filtered_rna_shared, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt (shared)', save=f'{exp}_rna_PostFilterShared.png')
    sc.pl.scatter(filtered_rna_shared, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt (shared)', save=f'{exp}_rna_PostFilterSharedlog.png')

    # DNA plots (shared)
    sc.pl.violin(filtered_me3_shared, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_shared_cells: {filtered_me3_shared.n_obs}\nMedian n_fragment: {np.median(filtered_me3_shared.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_me3_nfragPostFilterShared.png')
    plt.close()

    sc.pl.violin(filtered_ac_shared, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_shared_cells: {filtered_ac_shared.n_obs}\nMedian n_fragment: {np.median(filtered_ac_shared.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_ac_nfragPostFilterShared.png')
    plt.close()

    # Save filtered AnnData objects for shared barcodes for each modality
    filtered_rna_shared.write(f'{exp}_rna_CountFilteredShared.h5ad')
    filtered_ac_shared.write(f'{exp}_ac_CountFilteredShared.h5ad')
    filtered_me3_shared.write(f'{exp}_me3_CountFilteredShared.h5ad')

    snap.ex.export_fragments(filtered_ac_shared, groupby='Exp', prefix='', suffix='_FiltShared.bed.gz',compression_level=1)
    snap.ex.export_fragments(filtered_me3_shared, groupby='Exp', prefix='', suffix='_FiltShared.bed.gz',compression_level=1)



import matplotlib.pyplot as plt
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
from upsetplot import UpSet, plot, from_indicators
import os
import csv

if not os.path.exists('figures'):
    os.makedirs('figures')

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/CL/"
srcdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"
Experiment='H2R_CellLines'

#################################################
############# RNA ANALYSIS CELL LINES ###############
#################################################


cell_lines = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

# List to hold AnnData objects for each sample
adatas = []
sample_names = []

for sample in cell_lines:
    exp = sample
    try:
        adata_sample = sc.read_h5ad(f"{exp}_rna_CountFiltered.h5ad")
        adata_sample.obs['sample'] = exp  # Add sample name as a column for later coloring
        adatas.append(adata_sample)
        sample_names.append(exp)
    except Exception as e:
        print(f"Could not load {exp}_rna_CountFiltered.h5ad: {e}")

# Concatenate all AnnData objects
if len(adatas) == 0:
    raise ValueError("No RNA AnnData objects loaded for cell lines.")

adata = ad.concat(
    adatas,
    join="inner",
    label="sampleMerge",
    keys=sample_names,
    index_unique=None  # Keep original cell barcodes
)

# Store raw counts in a layer
adata.layers["counts"] = adata.X.copy()


# Normalize and log1p
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Cell cycle scoring
cell_cycle_genes = [x.strip() for x in open('/mnt/dataFast/ahrmad/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Highly variable genes (across all cells)
#sc.pp.highly_variable_genes(adata, batch_key="sample")
sc.pp.highly_variable_genes(
    adata,
    flavor='seurat_v3',
    n_top_genes=2000,
)
sc.pl.highly_variable_genes(adata, show=False, save=f'AllCellLines_HVG.pdf')

#sc.pp.scale(adata,zero_center=False, max_value=20) # optional clipping

# PCA
sc.tl.pca(adata,n_comps=30)
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True, show=False, save=f'AllCellLines_pca_variance_ratio.pdf')

# PCA plots
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"],
    dimensions=[(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)],
    ncols=3,
    show=False, save=f'AllCellLines_PCA1-2_QC.pdf'
)
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"],
    dimensions=[(2, 3), (2, 3), (2, 3), (2, 3), (2, 3)],
    ncols=3,
    show=False, save=f'AllCellLines_PCA3-4_QC.pdf'
)

# Neighbors/UMAP
#sc.pp.neighbors(adata, n_neighbors=30)
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, metric='cosine')  # or adjust as needed
sc.tl.umap(adata)

# UMAPs colored by QC and sample
sc.pl.umap(
    adata,
    color=["sample", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_50_genes", "log1p_total_counts"],
    ncols=3,
    show=False, save=f'AllCellLines_UMAP_QC.pdf'
)

sc.pl.umap(
    adata,
    color=['S_score', 'G2M_score', 'phase', 'sample'],
    ncols=2,
    show=False, save=f'AllCellLines_UMAP_QC_CellCycle.pdf'
)

sc.tl.leiden(adata, resolution=1)

sc.pl.umap(
    adata,
    color=["log1p_n_genes_by_counts", "log1p_total_counts","sample","sampleMerge","leiden"],
    legend_loc="on data",
    ncols=3,
    save=f'AllCellLines_UMAP_RNA_QC.pdf',
    show=False
)

adata = adata[adata.obs["leiden"] != "5"].copy()




sc.pp.highly_variable_genes(
    adata,
    flavor='seurat_v3',
    n_top_genes=2000,
)
sc.pl.highly_variable_genes(adata, show=False, save=f'AllCellLines_HVG_Filtered.pdf')

sc.tl.pca(adata,n_comps=30)
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True, show=False, save=f'AllCellLines_pca_variance_ratio_Filtered.pdf')

sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, metric='cosine')  # or adjust as needed
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5)

cluster2sample = {
    "0": "WSU",
    "1": "JJN2",
    "2": "Karpas422",
    "3": "GM12878",
}

# overwrite sampleMerge from leiden
adata.obs["sampleMerge"] = adata.obs["leiden"].map(cluster2sample).astype("category")


sc.pl.umap(
    adata,
    color=["log1p_n_genes_by_counts", "log1p_total_counts","sample","sampleMerge","leiden"],
    legend_loc="on data",
    ncols=3,
    save=f'AllCellLines_UMAP_RNA_QC_Filtered.pdf',
    show=False
)

# Loop over each cell line in sampleMerge
for exp in adata.obs["sampleMerge"].cat.categories:  # or: for exp in adata.obs["sampleMerge"].unique():
    # subset to that sample
    filtered_rna = adata[adata.obs["sampleMerge"] == exp]

    # get cell names (barcodes)
    cell_names = filtered_rna.obs_names

    # write to file: GM12878_CountFiltered_cell_names.lst, etc.
    out_fn = f"{exp}_CountFiltered_cell_names.lst"
    with open(out_fn, "w") as f:
        f.write("\n".join(cell_names))

    print(f"Wrote {len(cell_names)} cells to {out_fn}")


adata.write("AllCellLines_rna_merged_processed_Filtered.h5ad")





adata = adata[adata.obs["leiden"] != "6"].copy()

adata.obs["CL"] = adata.obs["sample"].str.split("_").str[-1]

sc.pl.umap(
    adata,
    color=["log1p_n_genes_by_counts", "log1p_total_counts","sample","sampleMerge","leiden","scDblFinder_score","scDblFinder_class"],
    legend_loc="on data",
    ncols=3,
    save=f'AllCellLines_UMAP_RNA_singlet_filtered.pdf',
    show=False
)

adata.write("AllCellLines_rna_merged_processed_Filtered.h5ad")








#################################################
############# DNA ANALYSIS HUMAN CELL LINES ###############
#################################################

#### Histone Mark Analysis 

import matplotlib.pyplot as plt
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
from upsetplot import UpSet, plot, from_indicators
import os
import csv

if not os.path.exists('figures'):
    os.makedirs('figures')

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/CL"

cell_lines = [
    "GM12878",
    "JJN2",
    "Karpas422",
    "WSU"
]

BL = '/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz'
GG = snap.genome.GRCh38

for moda in ['ac','me3']:
    
    qc_ext=f'_{moda}_CountFilteredShared.h5ad'
    Experiment=f"Cell_Lines_{moda}"

    if moda == 'ac':
        n=40000
        bs=5000
    elif moda == 'me3':
        n=30000
        bs=40000
    

    # Concatenation
    #Create concatenation list
    adata_list = []
    print(f'Loading AnnData...')
    for sample in cell_lines:
        print(f'Reading {sample}')
        a = snap.read(os.path.join(workdir,sample+qc_ext), backed=None)

        # Add metadata
        a.obs['cell_line'] = [sample] * a.n_obs

        # Change index
        a.obs.index = a.obs.index.astype(str) + '_' + sample

        # Export
        a_fname = f'{sample}_qcTOREMOVE.h5ad'
        a.write(a_fname)
        # Append to list for concatenation
        adata_list.append((sample,a_fname))
        del a

    # Exporting fragments
    print(f'Exporting fragments from all samples in {Experiment}...')
    adata = snap.AnnDataSet(adata_list,filename=f'{Experiment}.h5ad', add_key='sample')

    adata.obs['Exp'] = [Experiment]*adata.n_obs
    snap.ex.export_fragments(adata, groupby='Exp', prefix='', suffix=f'.tsv.gz',compression_level=1)
    adata.close()


    print(f'Concatenating')
    print(f'Samples: {cell_lines}')


    #Concatenation
    # building concat list
    adata_list = []
    for sample in cell_lines:
        print(f'Reading {sample}')
        a = snap.read(os.path.join(workdir,sample+qc_ext), backed=None)
        a.obs.index = a.obs.index.astype(str) + '_' + sample

        del a.obsm
        adata_list.append(a)
        del a

    print(f'Concatenating...')
    adata = ad.concat(adata_list, join='inner',label='sample',keys=cell_lines,index_unique=None)

    print(f'Analysis...')
    b = snap.pp.import_fragments(fragment_file=f'{Experiment}.tsv.gz',
        chrom_sizes=GG,
        sorted_by_barcode=True,min_num_fragments=0,
        tempdir='.')

    # Match cells
    common = adata.obs.index.intersection(b.obs.index)

    # check fraction of cells that survive the intersection
    frac_rna = len(common) / adata.shape[0]
    frac_atac = len(common) / b.shape[0]
    print("Fraction kept from RNA:", frac_rna)
    print("Fraction kept from ATAC:", frac_atac)

    # sanity check: does the sample column agree with the index suffix?
    suffix = pd.Index(adata.obs_names).str.rsplit('_', n=1).str[-1]
    # ideally -> all True

    adata = adata[common].copy()
    b = b[common].copy()

    # Get fragments,ref from b
    adata.obsm = b.obsm.copy()
    adata.uns['reference_sequences'] = b.uns['reference_sequences'].copy()

    # Build cell-by-peak matrix and store in pm
    #pm = snap.pp.make_peak_matrix(adata, peak_file=f'peaks/combined_H3K27ac_ITMPeaks.bed',counting_strategy='paired-insertion')
    snap.pp.add_tile_matrix(adata, bin_size=bs ,counting_strategy='paired-insertion', chunk_size=500000, inplace=True)
    # Copy fragments,ref from data
    del b


    # select features with mark-specific parameters
    snap.pp.select_features(adata, n_features=n, inplace=True, blacklist=BL)
    # Spectral embedding and UMAP
    snap.tl.spectral(adata, n_comps=8, weighted_by_sd=False, chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)
    snap.tl.umap(adata, use_rep='X_spectral', key_added='umap', random_state=None)
    snap.pp.knn(adata, n_neighbors=20, use_rep='X_spectral', method='kdtree')
    # Plots
    mark_name = f"H3K27{moda}"
    snap.pl.spectral_eigenvalues(adata, width=600, height=400, show=False, interactive=False, out_file=f'./figures/AllCellLines_{mark_name}_Eigenvalues.png')

    sc.pl.umap(
        adata,
        color=["n_fragment", "log10_n_fragment", "tsse", "sample"],
        ncols=2,
        save=f'AllCellLines_{mark_name}_UMAP_QC.pdf',
        show=False
    )
    snap.tl.leiden(adata, resolution=0.05)
    sc.pl.umap(
        adata,
        color=["leiden", "sample"],
        legend_loc="on data",
        ncols=2,
        save=f'AllCellLines_{mark_name}_UMAP_Leiden.pdf',
        show=False
    )
    # Save the merged AnnData object for this mark
    adata.write(f"AllCellLines_{mark_name}_merged_processed.h5ad")


adata = sc.read_h5ad(f"AllCellLines_H3K27ac_merged_processed.h5ad")

cluster2sample = {
    "0": "WSU",
    "1": "Karpas422",
    "2": "JJN2",
    "3": "GM12878",
}

# overwrite sampleMerge from leiden
adata.obs["sampleMerge"] = adata.obs["leiden"].map(cluster2sample).astype("category")
adata.write(f"AllCellLines_H3K27ac_merged_processed.h5ad")



adata = sc.read_h5ad(f"AllCellLines_H3K27me3_merged_processed.h5ad")

cluster2sample = {
    "0": "WSU",
    "1": "Karpas422",
    "2": "JJN2",
    "3": "GM12878",
}

# overwrite sampleMerge from leiden
adata.obs["sampleMerge"] = adata.obs["leiden"].map(cluster2sample).astype("category")
adata.write(f"AllCellLines_H3K27me3_merged_processed.h5ad")





#################################################
############# Linked UMAPs (scaled correctly) ###
#################################################

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from anndata import AnnData
from typing import Union

def plot_linked_umaps_procrustes(
    adata_rna: Union[str, AnnData],
    adata_ac:  Union[str, AnnData],
    adata_me3: Union[str, AnnData],
    ct_type: str,
    output_prefix: str = "./figures/LinkedUMAPs_RNA_Ac_Me3_lifted",
    fig_size=(18, 6),
    fig_dpi=300,
    lift_rna_frac: float = 0.20,
    lift_me3_frac: float = 0.20,
    gap_frac: float = 0.20,
):

    # ---------- helpers ----------
    def _ensure_adata(x):
        if isinstance(x, AnnData):
            return x
        elif isinstance(x, str):
            return sc.read_h5ad(x)
        else:
            raise TypeError("adata inputs must be AnnData or file paths to .h5ad")

    def _strip_suffix_after_last_underscore(adata, col_name="barcode_orig"):
        """
        Remove everything after the last '_' in obs_names.
        E.g. 'AAAC-1_GM12878_ac' -> 'AAAC-1'
        """
        adata = adata.copy()
        idx = pd.Index(adata.obs_names)
        adata.obs[col_name] = idx  # keep original if you ever need it
        new_idx = idx.str.rsplit("_", n=1).str[0]
        adata.obs_names = new_idx
        return adata

    def procrustes_align_to_base(target, base):
        X = np.asarray(base, dtype=float)
        Y = np.asarray(target, dtype=float)

        Xmean = X.mean(0, keepdims=True)
        Ymean = Y.mean(0, keepdims=True)
        Xc = X - Xmean
        Yc = Y - Ymean

        M = Yc.T @ Xc
        U, S, Vt = np.linalg.svd(M, full_matrices=False)
        R = U @ Vt
        s = S.sum() / (Yc**2).sum()

        Y_aligned = s * (Yc @ R) + Xmean
        return Y_aligned

    def x_extent(arr):
        return float(arr[:, 0].max() - arr[:, 0].min())

    def y_extent(arr):
        return float(arr[:, 1].max() - arr[:, 1].min())

    # ---------- I/O ----------
    a1 = _ensure_adata(adata_rna)  # RNA
    a2 = _ensure_adata(adata_ac)   # H3K27ac
    a3 = _ensure_adata(adata_me3)  # H3K27me3

    # REMOVE the "_something" suffix from ac & me3 obs_names
    a2 = _strip_suffix_after_last_underscore(a2)
    a3 = _strip_suffix_after_last_underscore(a3)

    # ---------- intersect cells ----------
    common = a1.obs_names.intersection(a2.obs_names).intersection(a3.obs_names)
    if len(common) == 0:
        raise ValueError(
            "No shared cell barcodes across RNA, H3K27ac and H3K27me3 "
            "even after stripping suffixes. Check that obs_names really match."
        )

    a1 = a1[common]
    a2 = a2[common]
    a3 = a3[common]

    umap1 = pd.DataFrame(a1.obsm['X_umap'], index=common, columns=['x1','y1'])
    umap2 = pd.DataFrame(a2.obsm['X_umap'], index=common, columns=['x2','y2'])
    umap3 = pd.DataFrame(a3.obsm['X_umap'], index=common, columns=['x3','y3'])

    # ---------- align ac & me3 to RNA ----------
    #umap2[['x2','y2']] = procrustes_align_to_base(
    #    umap2[['x2','y2']].values, umap1[['x1','y1']].values
    #)
    #umap3[['x3','y3']] = procrustes_align_to_base(
    #    umap3[['x3','y3']].values, umap1[['x1','y1']].values
    #)

    # ---------- spread panels horizontally ----------
    w = max(
        x_extent(umap1[['x1','y1']].values),
        x_extent(umap2[['x2','y2']].values),
        x_extent(umap3[['x3','y3']].values),
    )
    gap = gap_frac * w

    x1min = umap1['x1'].min()
    x2min = umap2['x2'].min()
    x3min = umap3['x3'].min()

    umap2['x2'] += (umap1['x1'].max() - x2min) + gap
    umap3['x3'] += (umap2['x2'].max() - x3min) + gap

    # ---------- vertical lift ----------
    H = max(
        y_extent(umap1[['x1','y1']].values),
        y_extent(umap2[['x2','y2']].values),
        y_extent(umap3[['x3','y3']].values),
    )
    umap1['y1'] += lift_rna_frac * H
    umap3['y3'] += lift_me3_frac * H

    # ---------- combine + colors ----------
    df = pd.concat([umap1, umap2, umap3], axis=1)
    df['cell_type'] = a1.obs[ct_type].astype('category')

    # colors matching your violin plot
    cell_type_colors = {
        "GM12878":  "#EF9649",
        "JJN2":     "#67AB58",
        "Karpas422":"#C9544F",
        "WSU":      "#9B7EBE",
    }

    to_rgba = {k: mcolors.to_rgba(v, alpha=0.22) for k, v in cell_type_colors.items()}
    fallback = mcolors.to_rgba("#BDBDBD", alpha=0.22)
    colors_arr = np.vstack([to_rgba.get(ct, fallback) for ct in df['cell_type'].astype(object)])

    # ---------- plot ----------
    fig, ax = plt.subplots(figsize=fig_size, dpi=fig_dpi)
    ax.set_title("RNA → H3K27ac → H3K27me3", fontsize=20, weight='bold')

    seg12 = np.stack([df[['x1','y1']].values, df[['x2','y2']].values], axis=1)
    seg23 = np.stack([df[['x2','y2']].values, df[['x3','y3']].values], axis=1)
    ax.add_collection(LineCollection(seg12, colors=colors_arr, linewidths=0.03, zorder=0))
    ax.add_collection(LineCollection(seg23, colors=colors_arr, linewidths=0.03, zorder=0))

    for (x, y, label_panel) in [('x1','y1','RNA'), ('x2','y2','H3K27ac'), ('x3','y3','H3K27me3')]:
        for ct, col in cell_type_colors.items():
            sub = df[df['cell_type'] == ct]
            if not len(sub):
                continue
            ax.scatter(
                sub[x], sub[y],
                s=16, c=col, edgecolors='none',
                alpha=0.75,
                label=ct if label_panel == 'RNA' else None,
                zorder=1, rasterized=False
            )

    ax.set_aspect('equal', adjustable='datalim')
    ax.margins(0.03)
    ax.axis('off')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    leg = ax.legend(
        by_label.values(), by_label.keys(), title='Cell Type',
        bbox_to_anchor=(1.01, 1), loc='upper left', markerscale=1.6, fontsize=9
    )
    leg.get_title().set_fontsize(10)

    plt.tight_layout()

    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pdf_path = f"{output_prefix}.pdf"
    png_path = f"{output_prefix}.png"
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, bbox_inches='tight')
    plt.close()

    return {"pdf": pdf_path, "png": png_path, "coords": df}


# ---- Example call (same files as your script) ----
results = plot_linked_umaps_procrustes(
     "AllCellLines_rna_merged_processed_Filtered.h5ad",
     "AllCellLines_H3K27ac_merged_processed.h5ad",
     "AllCellLines_H3K27me3_merged_processed.h5ad",
     ct_type="sampleMerge"  # or your column name
)

print(results["pdf"], results["png"])



#################################################
#################################################
#########################################
#################################################
#################################################
#########################################
#################################################
#################################################
#################### A20 MOUSE #####################
#################################################
#################################################
#########################################
#################################################
#################################################
#################################################
#########################################
#################################################

import matplotlib.pyplot as plt
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
from upsetplot import UpSet, plot, from_indicators
import os
import csv

if not os.path.exists('figures'):
    os.makedirs('figures')

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/CL/"
srcdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"

Experiment='H2R_CellLines'

#################################################
#################### DNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS DNA')
maxFrags = 80000

cell_lines = [
    "A20",
]

for sample in cell_lines:
    # H3K27ac
    bam_ac = f"Sc_CMD_{sample}_H3K27ac_NoDup.bam"
    exp_ac = sample + '_H3K27ac'
    bam_ac_fgt = exp_ac + '.tsv.gz'
    bam_ac_ad = exp_ac + '.h5ad'

    # H3K27me3
    bam_me3 = f"Sc_CMD_{sample}_H3K27me3_NoDup.bam"
    exp_me3 = sample + '_H3K27me3'
    bam_me3_fgt = exp_me3 + '.tsv.gz'
    bam_me3_ad = exp_me3 + '.h5ad'


    # Make fragment files
    ac_info = snap.pp.make_fragment_file(
        srcdir + bam_ac, workdir + bam_ac_fgt, is_paired=True, barcode_tag='CB',
        shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
        chrM=['chrM', 'M'], tempdir=workdir
    )
    with open(f'{exp_ac}.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        for key, value in ac_info.items():
            writer.writerow([key, value])

    me3_info = snap.pp.make_fragment_file(
        srcdir + bam_me3, workdir + bam_me3_fgt, is_paired=True, barcode_tag='CB',
        shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
        chrM=['chrM', 'M'], tempdir=workdir
    )
    with open(f'{exp_me3}.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        for key, value in me3_info.items():
            writer.writerow([key, value])

    # Import fragments
    ac_ad = snap.pp.import_fragments(
        workdir + bam_ac_fgt, snap.genome.GRCm39, min_num_fragments=0, sorted_by_barcode=True,
        chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64
    )
    me3_ad = snap.pp.import_fragments(
        workdir + bam_me3_fgt, snap.genome.GRCm39, min_num_fragments=0, sorted_by_barcode=True,
        chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64
    )

    # Compute metrics
    snap.metrics.tsse(ac_ad, gene_anno=snap.genome.GRCm39, inplace=True)
    snap.metrics.tsse(me3_ad, gene_anno=snap.genome.GRCm39, inplace=True)

    snap.pp.filter_cells(ac_ad, min_counts=100, max_counts=maxFrags, min_tsse=0, inplace=True)
    snap.pp.filter_cells(me3_ad, min_counts=100, max_counts=maxFrags, min_tsse=0, inplace=True)

    snap.metrics.frag_size_distr(ac_ad, add_key='frag_size_distr', inplace=True)
    snap.metrics.frag_size_distr(me3_ad, add_key='frag_size_distr', inplace=True)

    ac_ad.obs['log10_n_fragment'] = np.log10(ac_ad.obs['n_fragment'])
    me3_ad.obs['log10_n_fragment'] = np.log10(me3_ad.obs['n_fragment'])

    # Violin plot for ac
    sc.pl.violin(ac_ad, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {ac_ad.n_obs}\nMedian n_fragment:{np.median(ac_ad.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp_ac}_FragViolin.pdf')
    plt.close()

    # Violin plot for me3
    sc.pl.violin(me3_ad, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {me3_ad.n_obs}\nMedian n_fragment:{np.median(me3_ad.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp_me3}_FragViolin.pdf')
    plt.close()

    # Scatter plots
    sc.pl.scatter(ac_ad, "log10_n_fragment", "tsse", title=exp_ac, save=exp_ac + '_tsse.pdf')
    sc.pl.scatter(me3_ad, "log10_n_fragment", "tsse", title=exp_me3, save=exp_me3 + '_tsse.pdf')

    # SnapATAC2 QC plots
    snap.pl.tsse(ac_ad, min_fragment=0, width=750, height=600, interactive=False, show=False,
                 out_file=f"./figures/{exp_ac}_tsseDensity_PreFilter.pdf")
    snap.pl.frag_size_distr(ac_ad, width=750, height=600, interactive=False, show=False,
                            out_file=f"./figures/{exp_ac}_frag_size_distr_PreFilter.pdf")

    snap.pl.tsse(me3_ad, min_fragment=0, width=750, height=600, interactive=False, show=False,
                 out_file=f"./figures/{exp_me3}_tsseDensity_PreFilter.pdf")
    snap.pl.frag_size_distr(me3_ad, width=750, height=600, interactive=False, show=False,
                            out_file=f"./figures/{exp_me3}_frag_size_distr_PreFilter.pdf")

    # Add experiment name to obs
    ac_ad.obs['Exp'] = [exp_ac] * ac_ad.n_obs
    me3_ad.obs['Exp'] = [exp_me3] * me3_ad.n_obs

    # Export fragments and AnnData
    snap.ex.export_fragments(ac_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
    ac_ad.write(filename=f'{exp_ac}.h5ad')
    snap.ex.export_fragments(me3_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
    me3_ad.write(filename=f'{exp_me3}.h5ad')


#################################################
#################### RNA QC #######
#################################################

print(f'IMPORT AND QC PLOTS RNA FOR 8 CELL LINES')


cell_line_samples = ["A20"]

for sample in cell_line_samples:
    print(f"Processing {sample}")

    #adata = sc.read_csv(workdir + cnt_tsv, delimiter='\t', first_column_names=True)
    #adata = adata.T
    mtxPath = f"/mnt/dataFast/ahrmad/triseq_2025052/RNA_RE/{sample}.Solo.outGeneFull/filtered/"
    adata = sc.read_10x_mtx(mtxPath)

    # mitochondrial genes, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    # ribosomal genes (mouse: "Rps", "Rpl")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes (mouse: "Hba", "Hbb")
    adata.var["hb"] = adata.var_names.str.contains("^Hb[ab]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    # Violin plot for percent mt and top 50 genes
    sc.pl.violin(
        adata,
        ["pct_counts_mt", 'pct_counts_in_top_50_genes'],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        save=f'{sample}_Pct_PreFilter.pdf'
    )

    # Violin plot for log1p n_genes and total_counts
    sc.pl.violin(
        adata,
        ["log1p_n_genes_by_counts", "log1p_total_counts"],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        show=False  # Prevent auto-display so we can edit it
    )

    # Get the figure and axes
    fig = plt.gcf()
    axes = fig.axes

    # Loop through the axes to update y-axis tick labels
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]  # reverse log1p
        ax.set_yticklabels(new_labels)

    # Save or show the updated figure
    plt.savefig(f"./figures/{sample}_PreFilter.pdf")
    plt.close()

    # Scatter plots
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts",
        color="pct_counts_mt",
        title=sample + '\nColor: pct_counts_mt',
        save=f'{sample}_PreFilter.pdf'
    )
    sc.pl.scatter(
        adata, "log1p_total_counts", "log1p_n_genes_by_counts",
        color="pct_counts_mt",
        title=sample + '\nColor: pct_counts_mt',
        save=f'{sample}_PreFilterlog.pdf'
    )

    adata.write(f'{sample}.h5ad')


#################################################
############# RNA + DNA FILTERING ###############
#################################################

# Filtering parameters (max values unchanged)
maxFrags_me3 = 10000
maxFrags_ac = 12000
maxUMI = 6000
max_pct_50g = 40

cell_lines = [
    "A20",
]

for sample in cell_lines:
    print(f"\n=== Filtering and QC for {sample} ===")
    exp = sample

    # Set minFrags and minUMI per sample
    if "A20" in sample:
        minFrags_ac = 1000
        minFrags_me3 = 400
        max_pct_mt = 10
    else:
        max_pct_mt = 10
        minFrags_ac = 0
        minFrags_me3 = 0

    minUMI = 300  # For all samples

    # File names
    exp_rna = f"{sample}"
    exp_ac = f"{sample}_H3K27ac"
    exp_me3 = f"{sample}_H3K27me3"

    # Load data
    try:
        rna = sc.read_h5ad(f"{exp_rna}.h5ad")
    except Exception as e:
        print(f"Could not load {exp_rna}.h5ad: {e}")
        continue
    try:
        ac = sc.read_h5ad(f"{exp_ac}.h5ad")
    except Exception as e:
        print(f"Could not load {exp_ac}.h5ad: {e}")
        continue
    try:
        me3 = sc.read_h5ad(f"{exp_me3}.h5ad")
    except Exception as e:
        print(f"Could not load {exp_me3}.h5ad: {e}")
        continue

    # --- RNA filtering criteria ---
    filtered_rna = rna[rna.obs['total_counts'] >= minUMI]
    filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
    filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
    filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]

    # --- DNA filtering criteria ---
    filtered_ac = ac[snap.pp.filter_cells(ac, min_counts=minFrags_ac, max_counts=maxFrags_ac, min_tsse=0, inplace=False)]
    filtered_me3 = me3[snap.pp.filter_cells(me3, min_counts=minFrags_me3, max_counts=maxFrags_me3, min_tsse=0, inplace=False)]

    # Save filtered AnnData objects for each modality
    filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')
    filtered_ac.write(f'{exp}_ac_CountFiltered.h5ad')
    filtered_me3.write(f'{exp}_me3_CountFiltered.h5ad')


    # Output plots for each modality separately

    # RNA plots
    sc.pl.violin(
        filtered_rna,
        ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", 'pct_counts_in_top_50_genes'],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        save=f'{exp}_rna_PostFilter.png'
    )
    sc.pl.violin(
        filtered_rna,
        ["log1p_n_genes_by_counts", "log1p_total_counts"],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        show=False
    )
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
    plt.title(f"Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.obs['total_counts'])}")
    plt.savefig(f"./figures/{exp}_rna_PostFilterNoLog.pdf")
    plt.close()

    sc.pl.scatter(filtered_rna, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt', save=f'{exp}_rna_PostFilter.png')
    sc.pl.scatter(filtered_rna, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt', save=f'{exp}_rna_PostFilterlog.png')

    # DNA plots
    sc.pl.violin(filtered_me3, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {filtered_me3.n_obs}\nMedian n_fragment: {np.median(filtered_me3.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_me3_nfragPostFilter.png')
    plt.close()

    sc.pl.violin(filtered_ac, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {filtered_ac.n_obs}\nMedian n_fragment: {filtered_ac.obs['n_fragment'].median()}")
    plt.savefig(f'./figures/{exp}_ac_nfragPostFilter.png')
    plt.close()

    cell_names = filtered_rna.obs_names
    with open(f"{exp}_CountFiltered_cell_names.lst", "w") as f:
        f.write("\n".join(cell_names))

    # Now try to merge the modalities based on shared barcodes
    rna_barcodes = set(filtered_rna.obs.index)
    ac_barcodes = set(filtered_ac.obs.index)
    me3_barcodes = set(filtered_me3.obs.index)
    shared_barcodes = rna_barcodes & ac_barcodes & me3_barcodes

    print(f'Ac:{len(ac_barcodes)}\nMe3:{len(me3_barcodes)}\nRNA:{len(rna_barcodes)}\nShared:{len(shared_barcodes)}')

    # Create a Venn diagram with the sets
    plt.figure(figsize=(8, 8))
    venn = venn3([rna_barcodes, ac_barcodes, me3_barcodes],
                 set_labels=(f'RNA\n{minUMI}<UMIs<{maxUMI}\n<{max_pct_mt}% MT\n<{max_pct_50g}% 50TopGenes',
                             f'H3K27ac\n{minFrags_ac}<fragments<{maxFrags_ac}',
                             f'H3K27me3\n{minFrags_me3}<fragments<{maxFrags_me3}'))
    plt.suptitle(f'Venn Diagram of Passing Cells\nfor RNA, H3K27ac, and H3K27me3\n{sample}',
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'./figures/{exp}_VennPassingCells.pdf')
    plt.close()

    # Create a DataFrame to represent the sets
    all_items = rna_barcodes.union(ac_barcodes).union(me3_barcodes)
    data = {
        'RNA': [item in rna_barcodes for item in all_items],
        'H3K27ac': [item in ac_barcodes for item in all_items],
        'H3K27me3': [item in me3_barcodes for item in all_items],
    }
    df = pd.DataFrame(data, index=list(all_items))

    fig = plt.figure(figsize=(8, 10))
    plot(
        from_indicators(df.columns, data=df),
        sort_by='cardinality',
        show_percentages=True,
        min_subset_size="0.2%",
        facecolor="darkblue",
        fig=fig, element_size=None
    )
    plt.title(f"Passing Cells\nfor RNA, H3K27ac, and H3K27me3\n{sample}", fontsize=16, fontweight='bold')
    plt.savefig(f'./figures/{exp}_UpsetPassingCells.pdf')
    plt.close()

    if len(shared_barcodes) == 0:
        print(f"No shared barcodes for {sample}, skipping shared-barcode analysis.")
        continue

    # Subset the filtered datasets based on shared barcodes
    filtered_ac_shared = filtered_ac[filtered_ac.obs.index.isin(shared_barcodes)]
    filtered_me3_shared = filtered_me3[filtered_me3.obs.index.isin(shared_barcodes)]
    filtered_rna_shared = filtered_rna[filtered_rna.obs.index.isin(shared_barcodes)]

    print(f"Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.obs['total_counts'])}")

    cell_names = filtered_rna_shared.obs_names
    with open(f"{exp}_CountFilteredShared_cell_names.lst", "w") as f:
        f.write("\n".join(cell_names))

    # Output plots for each modality, for the shared barcodes

    # RNA plots (shared)
    sc.pl.violin(
        filtered_rna_shared,
        ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", 'pct_counts_in_top_50_genes'],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        save=f'{exp}_rna_PostFilterShared.png'
    )
    sc.pl.violin(
        filtered_rna_shared,
        ["log1p_n_genes_by_counts", "log1p_total_counts"],
        jitter=0.4,
        multi_panel=True,
        rotation=0.0000001,
        show=False
    )
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
    plt.title(f"Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.obs['total_counts'])}")
    plt.savefig(f"./figures/{exp}_rna_PostFilterSharedNoLog.pdf")
    plt.close()

    sc.pl.scatter(filtered_rna_shared, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt (shared)', save=f'{exp}_rna_PostFilterShared.png')
    sc.pl.scatter(filtered_rna_shared, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",
                  title=exp+'\nColor: pct_counts_mt (shared)', save=f'{exp}_rna_PostFilterSharedlog.png')

    # DNA plots (shared)
    sc.pl.violin(filtered_me3_shared, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_shared_cells: {filtered_me3_shared.n_obs}\nMedian n_fragment: {np.median(filtered_me3_shared.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_me3_nfragPostFilterShared.png')
    plt.close()

    sc.pl.violin(filtered_ac_shared, 'log10_n_fragment', rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_shared_cells: {filtered_ac_shared.n_obs}\nMedian n_fragment: {np.median(filtered_ac_shared.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_ac_nfragPostFilterShared.png')
    plt.close()

    # Save filtered AnnData objects for shared barcodes for each modality
    filtered_rna_shared.write(f'{exp}_rna_CountFilteredShared.h5ad')
    filtered_ac_shared.write(f'{exp}_ac_CountFilteredShared.h5ad')
    filtered_me3_shared.write(f'{exp}_me3_CountFilteredShared.h5ad')

    snap.ex.export_fragments(filtered_ac_shared, groupby='Exp', prefix='', suffix='_FiltShared.bed.gz',compression_level=1)
    snap.ex.export_fragments(filtered_me3_shared, groupby='Exp', prefix='', suffix='_FiltShared.bed.gz',compression_level=1)


