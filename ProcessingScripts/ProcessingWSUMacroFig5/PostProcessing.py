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
    
workdir = "/mnt/dataFast/ahrmad/triseq_202506/processed/"
Experiment='H2R_CD47i'

#################################################
#################### DNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS DNA')
maxFrags = 80000

# List of 6 human cell line sample names (prefixes)
cell_lines = [
    "Sc_MWD_WSU_CD47_Inhibitor",
    "Sc_MWD_WSU_Untreated",
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
        workdir + bam_ac_fgt, snap.genome.hg38, min_num_fragments=0, sorted_by_barcode=True,
        chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64
    )
    me3_ad = snap.pp.import_fragments(
        workdir + bam_me3_fgt, snap.genome.hg38, min_num_fragments=0, sorted_by_barcode=True,
        chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64
    )

    # Compute metrics
    snap.metrics.tsse(ac_ad, gene_anno=snap.genome.hg38, inplace=True)
    snap.metrics.tsse(me3_ad, gene_anno=snap.genome.hg38, inplace=True)

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
                 out_file=f"./figures/{exp_ac}_tsseDensity_PreFilter.png")
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
cell_lines = [
    "Sc_MWD_WSU_CD47_Inhibitor",
    "Sc_MWD_WSU_Untreated",
]

# List of 8 cell line sample names and their count tsv files
cell_line_samples = [
    ("Sc_MWD_WSU_CD47_Inhibitor", "Sc_MWr_WSU_CD47_Inhibitor_NoDupUMI_count.tsv"),
    ("Sc_MWD_WSU_Untreated", "Sc_MWr_WSU_Untreated_NoDupUMI_count.tsv"),
]

for sample, cnt_tsv in cell_line_samples:
    print(f"Processing {sample}")

    #adata = sc.read_csv(workdir + cnt_tsv, delimiter='\t', first_column_names=True)
    #adata = adata.T
    mtxPath = f"/mnt/dataFast/ahrmad/triseq_202506/RNA_reproc/{sample}.Solo.outGeneFull/filtered"
    adata = sc.read_10x_mtx(mtxPath)

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

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
maxFrags_me3 = 6000
maxFrags_ac = 6000
maxUMI = 6000
max_pct_50g = 30

# List of 6 human cell line sample names (prefixes)
cell_lines = [
    "Sc_MWD_WSU_CD47_Inhibitor",
    "Sc_MWD_WSU_Untreated",
]

for sample in cell_lines:
    print(f"\n=== Filtering and QC for {sample} ===")
    exp = sample

    # Set minFrags and minUMI per sample
    if "CD47_Inhibitor" in sample:
        minFrags_ac = 100
        minFrags_me3 = 150
        max_pct_mt = 12
        minUMI = 400  

    else:
        minFrags_ac = 250
        minFrags_me3 = 300
        max_pct_mt = 4
        minUMI = 400  

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
    filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_ribo'] < 3]


    # --- DNA filtering criteria ---
    filtered_ac = ac[snap.pp.filter_cells(ac, min_counts=minFrags_ac, max_counts=maxFrags_ac, min_tsse=0, inplace=False)]
    filtered_me3 = me3[snap.pp.filter_cells(me3, min_counts=minFrags_me3, max_counts=maxFrags_me3, min_tsse=0, inplace=False)]

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


    # Save filtered AnnData objects for each modality
    filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')
    filtered_ac.write(f'{exp}_ac_CountFiltered.h5ad')
    filtered_me3.write(f'{exp}_me3_CountFiltered.h5ad')


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

    # Subset the filtered datasets based on shared barcodes
    filtered_ac_shared = filtered_ac[filtered_ac.obs.index.isin(shared_barcodes)]
    filtered_me3_shared = filtered_me3[filtered_me3.obs.index.isin(shared_barcodes)]
    filtered_rna_shared = filtered_rna[filtered_rna.obs.index.isin(shared_barcodes)]

    print(f'Number of shared cells: {filtered_rna_shared.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna_shared.X.sum(axis=1))}')

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



