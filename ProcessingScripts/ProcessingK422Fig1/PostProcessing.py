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
    
workdir = "/mnt/dataFast/ahrmad/ScH2R_202409/processed/res/"

#################################################
#################### DNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS DNA')
maxFrags = 80000

bam_ac = 'H2RD_Human_H3K27ac_NoDup.bam'
exp_ac = bam_ac.split('_NoDup.bam')[0]
bam_ac_fgt = exp_ac + '.tsv.gz'
bam_ac_ad = exp_ac + '.h5ad'

bam_me3 = 'H2RD_Human_H3K27me3_NoDup.bam'
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
#################### RNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

workdir = "/mnt/dataFast/ahrmad/ScH2R_202409/processed/res/"
cnt_tsv = 'H2RR_Human_NoDupUMI_count.tsv'
sample='H2RRNA_Human'

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
############# RNA + DNA FILTERING ###############
#################################################

exp ='H2R_K422'

exp_rna='H2RRNA_Human'
exp_ac = 'H2RD_Human_H3K27ac'
exp_me3 = 'H2RD_Human_H3K27me3'

workdir = "/mnt/dataFast/ahrmad/ScH2R_202409/processed/res/"

#RNA
rna = sc.read_h5ad(f'{exp_rna}.h5ad')

#DNA
ac = sc.read_h5ad(f'{exp_ac}.h5ad')
me3 = sc.read_h5ad(f'{exp_me3}.h5ad')

#DNA
minFrags_me3 = 500
maxFrags_me3 = 15000
minFrags_ac = 250
maxFrags_ac = 12000
#RNA
minUMI = 750
maxUMI = 15000
max_pct_mt = 4
max_pct_50g = 30


# --- RNA filtering criteria ---
# Filter RNA dataset based on >100 UMI counts and <20% MT DNA
filtered_rna = rna[rna.obs['total_counts'] >= minUMI]
filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]

# --- DNA filtering criteria ---
filtered_ac = ac[snap.pp.filter_cells(ac, min_counts=minFrags_ac,max_counts=maxFrags_ac, min_tsse=0, inplace=False)]
filtered_me3 = me3[snap.pp.filter_cells(me3, min_counts=minFrags_me3,max_counts=maxFrags_me3, min_tsse=0, inplace=False)]

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


filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')
filtered_ac.write(f'{exp}_ac_CountFiltered.h5ad')
filtered_me3.write(f'{exp}_me3_CountFiltered.h5ad')





# BASH Filtering and Peak Calling
exps=("me3" "ac")
for e in "${exps[@]}"; do
    samtools view --threads 32 --with-header --read-group-file H2R_K422_${e}_CountFiltered_cell_names.lst --output H2R_K422_${e}_CountFiltered.bam H2RD_Human_H3K27${e}_NoDup.bam
    samtools index --threads 32 H2R_K422_${e}_CountFiltered.bam
done

#!/usr/bin/env bash
set -euo pipefail

export THREADS_PER_TASK=32

do_one() {
  e="$1"
  samtools view -@ "$THREADS_PER_TASK" \
    --with-header \
    --read-group-file "H2R_K422_${e}_CountFiltered_cell_names.lst" \
    -o "H2R_K422_${e}_CountFiltered.bam" \
    "H2RD_Human_H3K27${e}_NoDup.bam"

  samtools index -@ "$THREADS_PER_TASK" "H2R_K422_${e}_CountFiltered.bam"
}
export -f do_one

# Adjust -j to how many you want to run concurrently
parallel -j 70 do_one ::: me3 ac


bamCoverage -p 32 -bs 100 --extendReads --centerReads -b /mnt/dataFast/divyanshu/scnanocnt_deepseq/filteredBams/ScKDMA_S1_AGGCTATA_filtered.bam -o 10x_H3K27me3.bw -of bigwig --effectiveGenomeSize 2913022398
bamCoverage -p 32 -bs 100 --extendReads --centerReads -b /mnt/dataFast/divyanshu/scnanocnt_deepseq/filteredBams/ScKDMA_S1_GCCTCTAT_filtered.bam -o 10x_H3K27ac.bw -of bigwig --effectiveGenomeSize 2913022398


bamCoverage -p 32 -bs 100 --extendReads --centerReads -b /mnt/dataFast/ahrmad/ScH2R_202409/processed/res/H2R_K422_me3_CountFiltered.bam -o TrES-Seq_H3K27me3_rep1.bw -of bigwig --effectiveGenomeSize 2913022398
bamCoverage -p 32 -bs 100 --extendReads --centerReads -b /mnt/dataFast/ahrmad/ScH2R_202409/processed/res/H2R_K422_ac_CountFiltered.bam -o TrES-Seq_H3K27ac_rep1.bw -of bigwig --effectiveGenomeSize 2913022398


#################################################
############# COMPARE H3K27 Peaks ###############
#################################################

filtered_ac_acPeaks = sc.read_h5ad(f'{exp_ac}_PostQC.h5ad')
filtered_ac_me3Peaks = sc.read_h5ad(f'{exp_ac}_PostQC.h5ad')
filtered_me3_acPeaks = sc.read_h5ad(f'{exp_me3}_PostQC.h5ad')
filtered_me3_me3Peaks = sc.read_h5ad(f'{exp_me3}_PostQC.h5ad')

me3peak_file = f'./me3_macs3_q005_broad/me3_peaks_noBL.broadPeak'
me3peak_file = f'./me3_macs3_prev/me3_peaks_noBL.narrowPeak'
acpeak_file = f'./ac_macs3/ac_peaks_noBL.narrowPeak'

def remake_peak_matrix_keep_frag_ref(adata, peak_file, chunk_size=500000):
    # Save only the required metadata
    frag = adata.obsm.get("fragment_paired")
    refseq = adata.uns.get("reference_sequences")

    # Recreate peak matrix
    new_adata = snap.pp.make_peak_matrix(
        adata,
        peak_file=peak_file,
        counting_strategy="paired-insertion",
        chunk_size=chunk_size,
        inplace=False
    )

    # Restore metadata
    if frag is not None:
        new_adata.obsm["fragment_paired"] = frag
    if refseq is not None:
        new_adata.uns["reference_sequences"] = refseq

    return new_adata


# --- Apply to all four ---
filtered_me3_acPeaks = remake_peak_matrix_keep_frag_ref(filtered_me3_acPeaks, acpeak_file)
filtered_me3_me3Peaks = remake_peak_matrix_keep_frag_ref(filtered_me3_me3Peaks, me3peak_file)
filtered_ac_acPeaks = remake_peak_matrix_keep_frag_ref(filtered_ac_acPeaks, acpeak_file)
filtered_ac_me3Peaks = remake_peak_matrix_keep_frag_ref(filtered_ac_me3Peaks, me3peak_file)


snap.metrics.frip(filtered_me3_acPeaks, {"FRiP": acpeak_file}, inplace=True)
snap.metrics.frip(filtered_me3_me3Peaks, {"FRiP": me3peak_file}, inplace=True)
me3_me3Peaks = np.array(filtered_me3_me3Peaks.X.sum(axis=0)).flatten()
me3_acPeaks = np.array(filtered_me3_acPeaks.X.sum(axis=0)).flatten()

snap.metrics.frip(filtered_ac_acPeaks, {"FRiP": acpeak_file}, inplace=True)
snap.metrics.frip(filtered_ac_me3Peaks, {"FRiP": me3peak_file}, inplace=True)
ac_me3Peaks = np.array(filtered_ac_me3Peaks.X.sum(axis=0)).flatten()
ac_acPeaks = np.array(filtered_ac_acPeaks.X.sum(axis=0)).flatten()

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Normalize counts to counts per million (CPM)
def normalize_cpm_ac(data):
    return (data / np.sum(sc.read_h5ad(f'{exp_ac}_PostQC.h5ad').obs['n_fragment'].sum())) * 1e6
def normalize_cpm_me3(data):
    return (data / np.sum(sc.read_h5ad(f'{exp_me3}_PostQC.h5ad').obs['n_fragment'].sum())) * 1e6

me3_me3Peaks_cpm = normalize_cpm_me3(me3_me3Peaks)
ac_me3Peaks_cpm = normalize_cpm_ac(ac_me3Peaks)
me3_acPeaks_cpm = normalize_cpm_me3(me3_acPeaks)
ac_acPeaks_cpm = normalize_cpm_ac(ac_acPeaks)

# Combine all points for an overall correlation
x_combined = np.concatenate([me3_me3Peaks_cpm, me3_acPeaks_cpm])
y_combined = np.concatenate([ac_me3Peaks_cpm, ac_acPeaks_cpm])

# Calculate the overall correlation coefficient
overall_corr, _ = pearsonr(x_combined, y_combined)

# Plot setup
plt.figure(figsize=(8, 8))
sns.set(style="white")

# Scatter plot with pastel colors and transparency
plt.scatter(ac_me3Peaks_cpm, me3_me3Peaks_cpm, color='lightblue', edgecolor=None, alpha=0.6, label="H3K27me3 Peak")
plt.scatter(ac_acPeaks_cpm, me3_acPeaks_cpm, color='lightcoral', edgecolor=None, alpha=0.6, label="H3K27ac Peak")
plt.plot([0, 1], [0, 1], transform=plt.gca().transAxes, linestyle="--", color="gray", linewidth=1)

# Labels and titles
plt.ylabel("H3K27me3 (CPM)")
plt.xlabel("H3K27ac (CPM)")
plt.xlim(-5,400)
plt.ylim(-5,400)

plt.title(f"H3K27me3 and H3K27ac Peaks Correlation (r={overall_corr:.2f})")
plt.legend()
plt.tight_layout()

# Save plot
plt.savefig('./figures/P3_noBL.pdf')
plt.close()