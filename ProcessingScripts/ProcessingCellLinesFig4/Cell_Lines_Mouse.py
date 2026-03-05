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

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"

Experiment='H2R_CellLines'

#################################################
#################### DNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS DNA')
maxFrags = 80000

# List of mouse cell line sample names (prefixes)
cell_lines = [
    "Sc_CMD_A20",
]

for sample in cell_lines:
    # H3K27ac
    bam_ac = f"{sample}_H3K27ac_NoDup.bam"
    exp_ac = bam_ac.split('_NoDup.bam')[0]
    bam_ac_fgt = exp_ac + '.tsv.gz'
    bam_ac_ad = exp_ac + '.h5ad'

    # H3K27me3
    bam_me3 = f"{sample}_H3K27me3_NoDup.bam"
    exp_me3 = bam_me3.split('_NoDup.bam')[0]
    bam_me3_fgt = exp_me3 + '.tsv.gz'
    bam_me3_ad = exp_me3 + '.h5ad'

    # Make fragment files
    ac_info = snap.pp.make_fragment_file(
        workdir + bam_ac, workdir + bam_ac_fgt, is_paired=True, barcode_tag='CB',
        shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
        chrM=['chrM', 'M'], tempdir=workdir
    )
    with open(f'{exp_ac}.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        for key, value in ac_info.items():
            writer.writerow([key, value])

    me3_info = snap.pp.make_fragment_file(
        workdir + bam_me3, workdir + bam_me3_fgt, is_paired=True, barcode_tag='CB',
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
#################### RNA QC: Cell Lines #######
#################################################

print(f'IMPORT AND QC PLOTS RNA FOR CELL LINES')

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"

# List of cell line sample names and their count tsv files

cell_line_samples = [
    ("Sc_CMr_A20", "Sc_CMr_A20_NoDupUMI_count.tsv"),
]

for sample, cnt_tsv in cell_line_samples:
    print(f"Processing {sample}")

    adata = sc.read_csv(workdir + cnt_tsv, delimiter='\t', first_column_names=True)
    adata = adata.T

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
maxFrags_me3 = 12000
maxFrags_ac = 12000
maxUMI = 12000
max_pct_50g = 30

# List of mouse cell line sample names (prefixes)
cell_lines = [
    "Sc_CMD_A20",
]

for sample in cell_lines:
    print(f"\n=== Filtering and QC for {sample} ===")
    exp = sample

    # Set minFrags and minUMI per sample

    max_pct_mt = 10
    minFrags_ac = 1000
    minFrags_me3 = 500

    minUMI = 300  # For all samples

    # File names
    exp_rna = f"{sample.replace('CMD','CMr')}"
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
    plt.title(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.X.sum(axis=1))}')
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
        ax.set_title(f"n_cells: {filtered_ac.n_obs}\nMedian n_fragment: {np.median(filtered_ac.obs['n_fragment'])}")
    plt.savefig(f'./figures/{exp}_ac_nfragPostFilter.png')
    plt.close()

    # Output cell names for each modality
    for i, e in [(filtered_rna, "rna"), (filtered_ac, "ac"), (filtered_me3, "me3")]:
        cell_names = i.obs_names
        with open(f"{exp}_{e}_CountFiltered_cell_names.lst", "w") as f:
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

    print(f'Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.X.sum(axis=1))}')

    # Output cell names for shared set
    for i, e in [(filtered_rna_shared, "rna"), (filtered_ac_shared, "ac"), (filtered_me3_shared, "me3")]:
        cell_names = i.obs_names
        with open(f"{exp}_{e}_CountFilteredShared_cell_names.lst", "w") as f:
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
    plt.title(f'Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.X.sum(axis=1))}')
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


#################################################
############# RNA ANALYSIS MOUSE CELL LINES ###############
#################################################


cell_lines = [
    "Sc_CMD_A20",
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

# Initialize empty lists
s_genes = []
g2m_genes = []

# Read the CSV file
with open('/mnt/dataFast/ahrmad/Mouse_G2M_S_genes_Tirosh2016.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip header
    for row in reader:
        # Handle possible NA values
        if row[0] != "NA":
            s_genes.append(row[0])
        if row[1] != "NA":
            g2m_genes.append(row[1])

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Highly variable genes (across all cells)
sc.pp.highly_variable_genes(adata, batch_key="sample")
sc.pl.highly_variable_genes(adata, show=False, save=f'AllCellLines_Mouse_HVG.pdf')

# PCA
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True, show=False, save=f'AllCellLines_Mouse_pca_variance_ratio.pdf')

# PCA plots
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"],
    dimensions=[(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)],
    ncols=3,
    show=False, save=f'AllCellLines_Mouse_PCA1-2_QC.pdf'
)
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase", "sample"],
    dimensions=[(2, 3), (2, 3), (2, 3), (2, 3), (2, 3)],
    ncols=3,
    show=False, save=f'AllCellLines_Mouse_PCA3-4_QC.pdf'
)

# Neighbors/UMAP
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)

# UMAPs colored by QC and sample
sc.pl.umap(
    adata,
    color=["sample", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_50_genes", "log1p_total_counts"],
    ncols=3,
    show=False, save=f'AllCellLines_Mouse_UMAP_QC.pdf'
)
sc.pl.umap(
    adata,
    color=["log1p_total_counts", "log1p_n_genes_by_counts", "sample"],
    ncols=3,
    show=False, save=f'AllCellLines_Mouse_UMAP_QC_Log.pdf'
)
sc.pl.umap(
    adata,
    color=['S_score', 'G2M_score', 'phase', 'sample'],
    ncols=2,
    show=False, save=f'AllCellLines_Mouse_UMAP_QC_CellCycle.pdf'
)

sc.tl.leiden(adata, resolution=0.6)
sc.pl.umap(
    adata,
    color=['leiden', 'sample'],
    legend_loc=  "on data",
    ncols=2,
    show=False, save=f'AllCellLines_Mouse_UMAP_QC_RNA_Leiden.pdf'
)

# Save the resulting AnnData object
adata.write("AllCellLines_Mouse_rna_merged_processed.h5ad")

#################################################
############# VIOLIN RNA MOUSE CELL LINES ###############
#################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Use log1p counts for plotting
count_key = 'log1p_total_counts'

# Prepare DataFrame
df = adata.obs[['sample', count_key]].copy()
df.rename(columns={count_key: 'log1p_counts'}, inplace=True)

# Sample labels with cell counts
sample_counts = df['sample'].value_counts().sort_index()
df['sample_label'] = df['sample'].map(lambda s: f"{s}\n(n={sample_counts[s]})")
df['sample_label'] = pd.Categorical(df['sample_label'], categories=sorted(df['sample_label'].unique()), ordered=True)

# Plot
plt.figure(figsize=(1.5 * df['sample_label'].nunique(), 6))
sns.violinplot(data=df, x='sample_label', y='log1p_counts', inner=None, color='lightgray')
sns.boxplot(data=df, x='sample_label', y='log1p_counts', whis=1.5, width=0.2, fliersize=0, boxprops=dict(alpha=0.6))

# Set real-number Y-axis labels
yticks = plt.gca().get_yticks()
plt.gca().set_yticklabels([f"{int(np.expm1(y)):,}" for y in yticks])

# Aesthetics
plt.xticks(rotation=45, ha='right')
plt.xlabel("Sample")
plt.ylabel("Total Counts per Cell")
plt.title("Library Size per Cell (RNA)")
plt.tight_layout()
plt.savefig('./figures/ViolinRNA.pdf')




#################################################
############# DNA ANALYSIS MOUSE CELL LINES ###############
#################################################

#### Histone Mark Analysis for All Mouse Cell Lines (Joint UMAP per mark)

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

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"

cell_lines = [
    "Sc_CMD_A20",
]

BL = '/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz'
GG = snap.genome.GRCm39

# Set parameters for each mark as previously used
mark_params = {
    "ac": {"n": 200000, "bs": 1000},
    "me3": {"n": 200000, "bs": 5000}
}

for mark in ['me3','ac']:
    print(mark)
    n = mark_params[mark]["n"]
    bs = mark_params[mark]["bs"]
    adatas = []
    for sample in cell_lines:
        exp = sample
        ad_fn = f"{exp}_ac_CountFiltered.h5ad" if mark == 'ac' else f"{exp}_me3_CountFiltered.h5ad"
        try:
            adata_sample = sc.read_h5ad(ad_fn)
        except Exception as e:
            print(f"Could not load {ad_fn}: {e}")
            continue
        adata_sample.obs['sample'] = sample
        snap.pp.add_tile_matrix(adata_sample, bin_size=bs, counting_strategy='paired-insertion', chunk_size=500000, inplace=True)
        adatas.append(adata_sample)
    
    # Concatenate all samples for this mark
    adata = ad.concat(
        adatas,
        join="outer",
        label="sampleMerge",
        keys=cell_lines,
        index_unique=None
    )

    # select features with mark-specific parameters
    snap.pp.select_features(adata, n_features=n, inplace=True, blacklist=BL)
    # Spectral embedding and UMAP
    snap.tl.spectral(adata, n_comps=10, weighted_by_sd=True, chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)
    snap.tl.umap(adata, use_rep='X_spectral', key_added='umap', random_state=None)
    snap.pp.knn(adata, n_neighbors=10, use_rep='X_spectral', method='kdtree')
    # Plots
    mark_name = "H3K27ac" if mark == "ac" else "H3K27me3"
    snap.pl.spectral_eigenvalues(adata, width=600, height=400, show=False, interactive=False, out_file=f'./figures/AllCellLines_Mouse_{mark_name}_Eigenvalues.png')
    sc.pl.umap(
        adata,
        color=["n_fragment", "log10_n_fragment", "tsse", "sample"],
        ncols=2,
        save=f'AllCellLines_Mouse_{mark_name}_UMAP_QC.pdf',
        show=False
    )
    snap.tl.leiden(adata, resolution=0.5)
    sc.pl.umap(
        adata,
        color=["leiden", "sample"],
        legend_loc="on data",
        ncols=2,
        save=f'AllCellLines_Mouse_{mark_name}_UMAP_Leiden.pdf',
        show=False
    )
    # Save the merged AnnData object for this mark
    adata.write(f"AllCellLines_Mouse_{mark_name}_merged_processed.h5ad")
    jaccard_types(adata,mark_name,Experiment='Cell_Lines_H2R_202505')

#################################################
############# VIOLIN DNA MOUSE CELL LINES ###############
#################################################

for pp in ["H3K27ac","H3K27me3"]:
    adata = sc.read_h5ad(f"AllCellLines_Mouse_{pp}_merged_processed.h5ad")

    # Choose count column
    count_key = 'log10_n_fragment'

    # Prepare DataFrame
    df = adata.obs[['sample', count_key]].copy()
    df['sample'] = df['sample'].str.split("Sc_CMD_").str[1]

    df.rename(columns={count_key: 'log1p_counts'}, inplace=True)

    # Sample labels with cell counts
    sample_counts = df['sample'].value_counts().sort_index()
    df['sample_label'] = df['sample'].map(lambda s: f"{s}\n(n={sample_counts[s]})")
    df['sample_label'] = pd.Categorical(df['sample_label'], categories=sorted(df['sample_label'].unique()), ordered=True)

    # Plot
    plt.figure(figsize=(1.5 * df['sample_label'].nunique(), 6))
    sns.violinplot(data=df, x='sample_label', y='log1p_counts', inner=None, color='lightgray')
    sns.boxplot(data=df, x='sample_label', y='log1p_counts', whis=1.5, width=0.2, fliersize=0, boxprops=dict(alpha=0.6))

    # Set real-number Y-axis labels
    yticks = plt.gca().get_yticks()
    if 'H3' in pp:
        plt.gca().set_yticklabels([f"{int(10**y - 1):,}" for y in yticks])
    else:
        plt.gca().set_yticklabels([f"{int(np.expm1(y)):,}" for y in yticks])

    # Aesthetics
    plt.xticks(rotation=45, ha='right')
    plt.xlabel("Sample")
    plt.ylabel("Total Counts per Cell")
    plt.title(f"Library Size per Sample ({pp})")
    plt.tight_layout()
    plt.savefig(f'./figures/AllCellLines_Mouse_VIOLIN_{pp}.pdf')
    plt.close()

#################################################
############# DNA SHARED ###############
#################################################
adata1 = sc.read_h5ad(f"AllCellLines_Mouse_H3K27ac_merged_processed.h5ad")
adata2 = sc.read_h5ad(f"AllCellLines_Mouse_H3K27me3_merged_processed.h5ad")
adata3 = sc.read_h5ad("AllCellLines_Mouse_rna_merged_processed.h5ad")

exp="Cell_Lines_H2R_202505"
filtered_rna = adata3.copy()
filtered_ac = adata1.copy()
filtered_me3 = adata2.copy()

# Define the set of barcodes for each dataset after filtering
rna_barcodes = set(filtered_rna.obs.index)
ac_barcodes = set(filtered_ac.obs.index)
me3_barcodes = set(filtered_me3.obs.index)

# Get the intersection of barcodes across all datasets
shared_barcodes = rna_barcodes & ac_barcodes & me3_barcodes

print(f'Ac:{len(ac_barcodes)}\nMe3:{len(me3_barcodes)}\nRNA:{len(rna_barcodes)}\nShared:{len(shared_barcodes)}')

# Create a Venn diagram with the sets
plt.figure(figsize=(8, 8))
venn = venn3([rna_barcodes, ac_barcodes, me3_barcodes],
                 set_labels=(f'RNA', 
                         f'H3K27ac', 
                         f'H3K27me3'))
plt.suptitle('Venn Diagram of Passing Cells\nfor RNA, H3K27ac, and H3K27me3', 
          fontsize=16, fontweight='bold')
plt.tight_layout()

# Show the plot
plt.savefig(f'./figures/{exp}_VennPassingCells.pdf')
plt.close()

# Create a DataFrame to represent the sets
all_items = rna_barcodes.union(ac_barcodes).union(me3_barcodes)

data = {
    'RNA': [True if item in rna_barcodes else False for item in all_items],
    'H3K27ac': [True if item in ac_barcodes else False for item in all_items],
    'H3K27me3': [True if item in me3_barcodes else False for item in all_items],
}

# Create a DataFrame
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
plt.title(f"Passing Cells\nfor RNA, H3K27ac, and H3K27me3", fontsize=16, fontweight='bold')
plt.savefig(f'./figures/{exp}_UpsetPassingCells.pdf')
plt.close()

shared_barcodes
for adata,pp in [(adata1,'H3K27ac'),(adata2,'H3K27me3'),(adata3,'RNA')]:
    plwefd = adata[adata.obs.index.isin(shared_barcodes).copy()]
    sc.pl.umap(
        plwefd,
        color=['leiden', 'sample'],
        legend_loc=  "on data",
        ncols=2,
        show=False, save=f'AllCellLines_Mouse_UMAP_QC_{pp}_Leiden_SHARED.pdf'
    )


for adataPP,pp in [(adata1,'H3K27ac'),(adata2,'H3K27me3'),(adata3,'RNA')]:
    adata = adataPP[adataPP.obs.index.isin(shared_barcodes).copy()]

    # Choose count column
    if 'H3' in pp:
        count_key = 'log10_n_fragment'
    else:
        count_key = 'log1p_total_counts'

    # Prepare DataFrame
    df = adata.obs[['sample', count_key]].copy()
    df['sample'] = df['sample'].str.split("Sc_CMD_").str[1]

    df.rename(columns={count_key: 'log1p_counts'}, inplace=True)

    # Sample labels with cell counts
    sample_counts = df['sample'].value_counts().sort_index()
    df['sample_label'] = df['sample'].map(lambda s: f"{s}\n(n={sample_counts[s]})")
    df['sample_label'] = pd.Categorical(df['sample_label'], categories=sorted(df['sample_label'].unique()), ordered=True)

    # Plot
    plt.figure(figsize=(1.5 * df['sample_label'].nunique(), 6))
    sns.violinplot(data=df, x='sample_label', y='log1p_counts', inner=None, color='lightgray')
    sns.boxplot(data=df, x='sample_label', y='log1p_counts', whis=1.5, width=0.2, fliersize=0, boxprops=dict(alpha=0.6))

    # Set real-number Y-axis labels
    yticks = plt.gca().get_yticks()
    if 'H3' in pp:
        plt.gca().set_yticklabels([f"{int(10**y - 1):,}" for y in yticks])
    else:
        plt.gca().set_yticklabels([f"{int(np.expm1(y)):,}" for y in yticks])

    # Aesthetics
    plt.xticks(rotation=45, ha='right')
    plt.xlabel("Sample")
    plt.ylabel("Total Counts per Cell")
    plt.title(f"Library Size per Sample ({pp})")
    plt.tight_layout()
    plt.savefig(f'./figures/AllCellLines_Mouse_VIOLIN_{pp}_SHARED.pdf')
    plt.close()
