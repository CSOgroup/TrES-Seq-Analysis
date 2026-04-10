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
    
workdir = "/mnt/dataFast/ahrmad/Tonsil2/"

#################################################
#################### DNA QC #####################
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
    
workdir = "/mnt/dataFast/ahrmad/Tonsil2/"

#################################################
#################### DNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS DNA')
maxFrags = 80000

bam_ac = 'Sc_VTD_HumanTonsil2_H3K27ac_NoDup.bam'
exp_ac = bam_ac.split('_NoDup.bam')[0]
bam_ac_fgt = exp_ac + '.tsv.gz'
bam_ac_ad = exp_ac + '.h5ad'

bam_me3 = 'Sc_VTD_HumanTonsil2_H3K27me3_NoDup.bam'
exp_me3 = bam_me3.split('_NoDup.bam')[0]
bam_me3_fgt = exp_me3 + '.tsv.gz'
bam_me3_ad = exp_me3 + '.h5ad'

ac_info = snap.pp.make_fragment_file(workdir+bam_ac,workdir+bam_ac_fgt, is_paired=True, barcode_tag='CB',
                            shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
                            chrM=['chrM', 'M'], tempdir=workdir)

with open(f'{exp_ac}.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write data
    for key, value in ac_info.items():
        writer.writerow([key, value])

me3_info = snap.pp.make_fragment_file(workdir+bam_me3, workdir+bam_me3_fgt, is_paired=True, barcode_tag='CB',
                            shift_left=4, shift_right=-5, min_mapq=20, chunk_size=500000000,
                           chrM=['chrM', 'M'], tempdir=workdir)

with open(f'{exp_me3}.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write data
    for key, value in me3_info.items():
        writer.writerow([key, value])


ac_ad = snap.pp.import_fragments(workdir+bam_ac_fgt, snap.genome.hg38, min_num_fragments=0, sorted_by_barcode=True,
 chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64)

me3_ad = snap.pp.import_fragments(workdir+bam_me3_fgt, snap.genome.hg38, min_num_fragments=0, sorted_by_barcode=True,
 chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir=workdir, backend='hdf5', n_jobs=64)

snap.metrics.tsse(ac_ad, gene_anno=snap.genome.hg38, inplace=True)
snap.metrics.tsse(me3_ad, gene_anno=snap.genome.hg38, inplace=True)

snap.pp.filter_cells(ac_ad, min_counts=100,max_counts=maxFrags, min_tsse=0, inplace=True)
snap.pp.filter_cells(me3_ad, min_counts=100,max_counts=maxFrags, min_tsse=0, inplace=True)

snap.metrics.frag_size_distr(ac_ad,add_key='frag_size_distr', inplace=True)
snap.metrics.frag_size_distr(me3_ad,add_key='frag_size_distr', inplace=True)

ac_ad.obs['log10_n_fragment'] = np.log10(ac_ad.obs['n_fragment'])
me3_ad.obs['log10_n_fragment'] = np.log10(me3_ad.obs['n_fragment'])

# Create the violin plot
sc.pl.violin(ac_ad, 'log10_n_fragment',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {ac_ad.n_obs}\nMedian n_fragment:{np.median(ac_ad.obs['n_fragment'])}")


# Save or show the updated figure
plt.savefig(f'./figures/{exp_ac}_FragViolin.pdf')
plt.close()



# Create the violin plot
sc.pl.violin(me3_ad, 'log10_n_fragment',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {me3_ad.n_obs}\nMedian n_fragment:{np.median(me3_ad.obs['n_fragment'])}")


# Save or show the updated figure
plt.savefig(f'./figures/{exp_me3}_FragViolin.pdf')
plt.close()



sc.pl.scatter(ac_ad, "log10_n_fragment", "tsse",title=exp_ac,save=exp_ac+'_tsse.pdf')
sc.pl.scatter(me3_ad, "log10_n_fragment", "tsse",title=exp_me3,save=exp_me3+'_tsse.pdf')

snap.pl.tsse(ac_ad, min_fragment=0, width=750, height=600, interactive=False, show=False, out_file=f"./figures/{exp_ac}_tsseDensity_PreFilter.pdf")
snap.pl.frag_size_distr(ac_ad, width=750, height=600, interactive=False, show=False, out_file=f"./figures/{exp_ac}_frag_size_distr_PreFilter.pdf")

snap.pl.tsse(me3_ad, min_fragment=0, width=750, height=600, interactive=False, show=False, out_file=f"./figures/{exp_me3}_tsseDensity_PreFilter.pdf")
snap.pl.frag_size_distr(me3_ad, width=750, height=600, interactive=False, show=False, out_file=f"./figures/{exp_me3}_frag_size_distr_PreFilter.pdf")

ac_ad.obs['Exp'] = [exp_ac]*ac_ad.n_obs
me3_ad.obs['Exp'] = [exp_me3]*me3_ad.n_obs


snap.ex.export_fragments(ac_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
ac_ad.write(filename=f'{exp_ac}.h5ad')
snap.ex.export_fragments(me3_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
me3_ad.write(filename=f'{exp_me3}.h5ad')


#################################################
#################### RNA QC (OLD) #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

workdir = "/mnt/dataFast/ahrmad/Tonsil2/"
cnt_tsv = 'Sc_VTR_HumanTonsil2_NoDupUMI_count.tsv'
sample='Sc_VTR_HumanTonsil2'

adata = sc.read_csv(workdir+cnt_tsv,delimiter='\t',first_column_names=True)
adata = adata.T

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
#################### RNA QC 2 #####################
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
    
workdir = "/mnt/dataFast/ahrmad/Tonsil2/"

#################################################
#################### RNA QC IMPORT #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

sample='Sc_VTR_HumanTonsil2'

#workdir = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/"
#cnt_tsv = 'Sc_VTr10_S2_Sc_VTr12_S3_VTHumanTonsil_NoDupUMI_count.tsv'
#adata = sc.read_csv(workdir+cnt_tsv,delimiter='\t',first_column_names=True)
#adata = adata.T

#if cellbender
#adata = sc.read_10x_h5("RNA_CB_filtered.h5")

#gzip everything
mtxPath = "/mnt/dataFast/ahrmad/Tonsil2/Sc_VTR_HumanTonsil2.Solo.outGeneFull/filtered/"
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

exp ='Sc_HumanTonsil2'

exp_rna='Sc_VTR_HumanTonsil2'

#RNA
rna = sc.read_h5ad(f'{exp_rna}.h5ad')

#DNA
minFrags_me3 = 400
maxFrags_me3 = 10000
minFrags_ac = 400
maxFrags_ac = 10000
#RNA
minUMI = 500
maxUMI = 8000
max_pct_mt = 5
max_pct_hb = 5
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


#filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')



#################################################
############# RNA + DNA FILTERING ###############
#################################################

exp ='Sc_HumanTonsil2'

exp_rna='Sc_VTR_HumanTonsil2'
exp_ac = 'Sc_VTD_HumanTonsil2_H3K27ac'
exp_me3 = 'Sc_VTD_HumanTonsil2_H3K27me3'

workdir = "/mnt/dataFast/ahrmad/Tonsil2/"

#RNA
rna = sc.read_h5ad(f'{exp_rna}.h5ad')

#DNA
ac = sc.read_h5ad(f'{exp_ac}.h5ad') 
me3 = sc.read_h5ad(f'{exp_me3}.h5ad')




# --- RNA filtering criteria ---
# Filter RNA dataset based on >100 UMI counts and <20% MT DNA
filtered_rna = rna[rna.obs['total_counts'] >= minUMI]
filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_hb'] < max_pct_hb]

# --- DNA filtering criteria ---
filtered_ac = ac[snap.pp.filter_cells(ac, min_counts=minFrags_ac,max_counts=maxFrags_ac, min_tsse=0, inplace=False)]
filtered_me3 = me3[snap.pp.filter_cells(me3, min_counts=minFrags_me3,max_counts=maxFrags_me3, min_tsse=0, inplace=False)]

filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')
filtered_ac.write(f'{exp}_ac_CountFiltered.h5ad')
filtered_me3.write(f'{exp}_me3_CountFiltered.h5ad')


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
             set_labels=(f'RNA\n{minUMI}<UMIs<{maxUMI}\n<{max_pct_mt}% MT\n<{max_pct_50g}% 50TopGenes', 
                         f'H3K27ac\n{minFrags_ac}<fragments<{maxFrags_ac}', 
                         f'H3K27me3\n{minFrags_me3}<fragments<{maxFrags_me3}'))
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


# Subset the filtered datasets based on shared barcodes
filtered_ac = filtered_ac[filtered_ac.obs.index.isin(shared_barcodes)]
filtered_me3 = filtered_me3[filtered_me3.obs.index.isin(shared_barcodes)]
filtered_rna = filtered_rna[filtered_rna.obs.index.isin(shared_barcodes)]

print(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.X.sum(axis=1))}')

for i,e in  [(filtered_rna,"rna"),(filtered_ac,"ac"),(filtered_me3,"me3")]:
    cell_names = i.obs_names
    # Write directly to a text file, one name per line
    with open(f"{exp}_{e}_CountFiltered_cell_names.lst", "w") as f:
        f.write("\n".join(cell_names))

merged_ad = filtered_rna.copy()
merged_ad.obs = merged_ad.obs.reindex(index=sorted(merged_ad.obs.index))
ac_c = filtered_ac.copy()
ac_c.obs = ac_c.obs.reindex(index=sorted(ac_c.obs.index))
me3_c = filtered_me3.copy()
me3_c.obs = me3_c.obs.reindex(index=sorted(me3_c.obs.index))
assert me3_c.obs_names.equals(ac_c.obs_names) and me3_c.obs_names.equals(merged_ad.obs_names), "The obs_names do not match across me3_c, ac_c, and merged_ad"
ac_obs = ac_c.obs.add_suffix('_ac')
me3_obs = me3_c.obs.add_suffix('_me3')
merged_ad.obs = merged_ad.obs.join([ac_obs, me3_obs])

#sc.pl.scatter(merged_ad, "log10_n_fragment_me3", "log10_n_fragment_ac",color='log1p_total_counts',title="Log10 Unique Fragments in Histone Marks\nColored by log1p RNA UMIs",save=f'{exp}_Counts3.png')

sc.pl.scatter(merged_ad, "log10_n_fragment_me3", "log10_n_fragment_ac",color='log1p_total_counts',title="Log10 Unique Fragments in Histone Marks\nColored by log1p RNA UMIs",show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)

# Loop through the axes to update y-axis tick labels
for ax in axes:
    xticks = ax.get_xticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in xticks]
    ax.set_xticklabels(new_labels)


# Save or show the updated figure
plt.savefig(f"./figures/{exp}_Counts3.png")
plt.close()

#sc.pl.scatter(merged_ad, "tsse_ac", "tsse_me3",color='log1p_total_counts',title="TSSe in Histone Marks\nColored by log1p RNA Counts",save=f'{exp}_TSSe_RNACounts.png')

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

sc.pl.violin(filtered_me3, 'log10_n_fragment',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {filtered_me3.n_obs}\nMedian n_fragment: {np.median(filtered_me3.obs['n_fragment'])}")


# Save or show the updated figure
plt.savefig(f'./figures/{exp}_me3_nfragPostFilter.png')
plt.close()


sc.pl.violin(filtered_ac, 'log10_n_fragment',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {filtered_ac.n_obs}\nMedian n_fragment: {np.median(filtered_ac.obs['n_fragment'])}")


# Save or show the updated figure
plt.savefig(f'./figures/{exp}_ac_nfragPostFilter.png')
plt.close()


filtered_rna.write(f'{exp}_rna_CountFilteredShared.h5ad')
filtered_ac.write(f'{exp}_ac_CountFilteredShared.h5ad')
filtered_me3.write(f'{exp}_me3_CountFilteredShared.h5ad')


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def _backtransform(vals, modality):
    """Convert plotted values back to raw counts."""
    if 'H3' in modality:
        # assuming log10(n + 1)
        return (10 ** vals) - 1
    else:
        # assuming natural log1p
        return np.expm1(vals)

def vio(adata, modality):
    """
    Violin plot per sample for a given AnnData and modality.
    Title is just the modality. Shows per-sample median (raw counts) as a dot + label.
    """
    # choose count key
    if 'H3' in modality:
        count_key = 'log10_n_fragment'
        y_label = "Fragments per cell"
    else:
        count_key = 'log1p_total_counts'
        y_label = "UMI counts per cell"

    if count_key not in adata.obs.columns:
        raise KeyError(f"Expected '{count_key}' in adata.obs for modality '{modality}'")
    adata.obs['sample']="Tonsil2"
    # prepare dataframe
    df = adata.obs[['sample', count_key]].copy()
    # make sample names compact if prefixed like "Sc_MWD_WSU_XXX"
    df['sample'] = df['sample'].astype(str)
    if df['sample'].str.contains("Sc_MWD_WSU_").any():
        df['sample'] = df['sample'].str.split("Sc_MWD_WSU_").str[-1]

    df.rename(columns={count_key: 'metric'}, inplace=True)

    # sample labels with n
    sample_counts = df['sample'].value_counts().sort_index()
    df['sample_label'] = df['sample'].map(lambda s: f"{s}\n(n={sample_counts[s]})")
    cats = sorted(df['sample_label'].unique(), key=lambda x: (x.split("\n")[0]))
    df['sample_label'] = pd.Categorical(df['sample_label'], categories=cats, ordered=True)

    # figure
    width = max(4, 1.3 * df['sample_label'].nunique())
    plt.figure(figsize=(width, 4))

    # violin + box
    sns.violinplot(data=df, x='sample_label', y='metric', inner=None, cut=0, linewidth=0.5, color='lightgray')
    sns.boxplot(data=df, x='sample_label', y='metric', whis=1.5, width=0.18, fliersize=0, boxprops=dict(alpha=0.6))

    ax = plt.gca()

    # annotate medians (dot at median on plotted scale, text with raw median)
    med_plot = df.groupby('sample_label')['metric'].median()
    med_raw = _backtransform(med_plot.values, modality)
    for i, (cat, y_med_plot, y_med_raw) in enumerate(zip(med_plot.index, med_plot.values, med_raw)):
        # dot at median
        ax.plot(i, y_med_plot, marker='o', markersize=5, mec='black', mfc='white', zorder=5)
        # label slightly above the dot
        ax.text(i, y_med_plot, f"{int(round(y_med_raw)):,}", ha='center', va='bottom', fontsize=8, rotation=0, zorder=6)

    # y-axis ticks to raw counts
    yticks = ax.get_yticks()
    if 'H3' in modality:
        ax.set_yticklabels([f"{int(max(0, round(10**y - 1))):,}" for y in yticks])
    else:
        ax.set_yticklabels([f"{int(max(0, round(np.expm1(y)))):,}" for y in yticks])

    # aesthetics
    plt.xticks(rotation=45, ha='right')
    plt.xlabel("")
    plt.ylabel(y_label)
    plt.title(modality)
    plt.tight_layout()
    outpath = f'./figures/Type_VIOLIN_{modality}.pdf'
    plt.savefig(outpath)
    plt.close()
    print(f"Saved: {outpath}")

# --- calls for the three subsetted AnnData objects ---
# Expecting you already subset to common barcodes:
# adata1 = adata1[common_barcodes].copy()
# adata2 = adata2[common_barcodes].copy()
# adata3 = adata3[common_barcodes].copy()

vio(filtered_rna, "RNA")
vio(filtered_ac, "H3K27ac")
vio(filtered_me3, "H3K27me3")

