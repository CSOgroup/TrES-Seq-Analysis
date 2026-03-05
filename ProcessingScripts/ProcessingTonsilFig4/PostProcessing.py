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
    
workdir = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/"

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
#################### RNA QC IMPORT #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

sample='Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil'

#gzip everything
mtxPath = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/RNA_ONLY/STARsoloFILT/"
#mtxPath = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/RNA_ONLY/STARsoloRAW/"
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

sc.pl.violin(
    adata,
    ["pct_counts_mt", 'pct_counts_in_top_50_genes'],
    jitter=0.4,
    multi_panel=True
    ,rotation=0.0000001,save=f'{sample}_Pct_PreFilter.pdf'
)

# Create the violin plot
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

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",title=sample+'\nColor: pct_counts_mt',save=f'{sample}_PreFilter.pdf')
sc.pl.scatter(adata, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",title=sample+'\nColor: pct_counts_mt',save=f'{sample}_PreFilterlog.pdf')

adata.write(f'{sample}.h5ad')

#################################################
############# RNA FILTERING ###############
#################################################

exp ='Sc_VTHumanTonsil'

exp_rna='Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil'

#RNA
rna = sc.read_h5ad(f'{exp_rna}.h5ad')

#RNA
minUMI = 500
maxUMI = 8000
max_pct_mt = 2.5
max_pct_50g = 25


# --- RNA filtering criteria ---
filtered_rna = rna[rna.obs['total_counts'] >= minUMI]
filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]

sc.pl.violin(
    filtered_rna,
    ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", 'pct_counts_in_top_50_genes'],
    jitter=0.4,
    multi_panel=True
    ,rotation=0.0000001,save=f'{exp}_rna_PostFilter.png'
)

# Create the violin plot
sc.pl.violin(
    filtered_rna,
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
plt.suptitle(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.obs["total_counts"])}')
plt.savefig(f"./figures/{exp}_rna_PostFilterNoLog.pdf")
plt.close()


sc.pl.scatter(filtered_rna, "total_counts", "n_genes_by_counts", color="pct_counts_mt",title=exp+'\nColor: pct_counts_mt',save=f'{exp}_rna_PostFilter.png')
sc.pl.scatter(filtered_rna, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",title=exp+'\nColor: pct_counts_mt',save=f'{exp}_rna_PostFilterlog.png')


filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')


