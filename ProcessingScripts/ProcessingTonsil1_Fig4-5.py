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

#################################################
#################### RNA QC IMPORT #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"
cnt_tsv = 'Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil_NoDupUMI_count.tsv'
sample='Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil'

#adata = sc.read_csv(workdir+cnt_tsv,delimiter='\t',first_column_names=True)
#adata = adata.T

adata = sc.read_10x_h5("RNA_CB_filtered.h5")


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


cellbender remove-background \
     --input Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil.h5ad \
     --output RNA_CB.h5 \
     --expected-cells 6000 \
     --total-droplets-included 15000 \
     --cpu-threads 30


#################################################
############# RNA FILTERING ###############
#################################################

exp ='Sc_VTHumanTonsil'

exp_rna='Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil'

workdir = "/mnt/dataFast/ahrmad/triseq_2025052/processed/"

#RNA
rna = sc.read_h5ad(f'{exp_rna}.h5ad')

#RNA
minUMI = 500
maxUMI = 12000
max_pct_mt = 5
max_pct_50g = 30


# --- RNA filtering criteria ---
# Filter RNA dataset based on >100 UMI counts and <20% MT DNA
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
plt.title(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.X.sum(axis=1))}')
plt.savefig(f"./figures/{exp}_rna_PostFilterNoLog.pdf")
plt.close()


sc.pl.scatter(filtered_rna, "total_counts", "n_genes_by_counts", color="pct_counts_mt",title=exp+'\nColor: pct_counts_mt',save=f'{exp}_rna_PostFilter.png')
sc.pl.scatter(filtered_rna, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",title=exp+'\nColor: pct_counts_mt',save=f'{exp}_rna_PostFilterlog.png')



filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')