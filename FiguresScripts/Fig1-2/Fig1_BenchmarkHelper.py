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
    
workdir = "./"

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


#snap.ex.export_fragments(ac_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
ac_ad.write(filename=f'{exp_ac}.h5ad')
#snap.ex.export_fragments(me3_ad, groupby='Exp', prefix='', suffix='.tsv.gz')
me3_ad.write(filename=f'{exp_me3}.h5ad')

#################################################
#################### RNA QC #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

workdir = "./"
sample='H2RRNA_Human'

#cnt_tsv = 'H2RR_Human_NoDupUMI_count.tsv'
#adata = sc.read_csv(workdir+cnt_tsv,delimiter='\t',first_column_names=True)
#adata = adata.T

mtxPath = "/mnt/dataFast/ahrmad/ScH2R_202409/RNA2/H2RRNA_Human.Solo.outGeneFull/filtered"
adata = sc.read_10x_mtx(mtxPath)
adata.obs_names = adata.obs_names.str[-4:] + adata.obs_names.str[:-4]

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
maxFrags_ac = 15000
#RNA
minUMI = 1000
maxUMI = 25000
max_pct_mt = 10
max_pct_hb = 5
max_pct_50g = 30


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





###################
###################
###################
#MTC HELPER DNA/RNA####
###################
###################
###################

# snap_import.py
import os
from glob import glob
import snapatac2 as snap

# ---- user knobs ----
N_JOBS   = 64                     # CPU cores to use for parallel import
MIN_FRAGS = 100                  # keep all cells for now; adjust later if you want
# ---------------------

def sample_key_from_path(p: str) -> str:
    """Make a readable per-sample key, e.g. GSM8529544_S1"""
    b = os.path.basename(p)
    parts = b.split("_")
    # Expecting: GSM8529544_HIST_S1_H3K27ac.fragments.tsv.gz
    return f"{parts[0]}_{parts[2]}" if len(parts) >= 3 else os.path.splitext(b)[0]

def import_modality(pattern: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    frag_files = sorted(glob(pattern))
    assert frag_files, f"No files matched: {pattern}"

    # one output .h5ad per fragment file
    out_files = [
        os.path.join(out_dir, os.path.basename(f).replace(".fragments.tsv.gz", ".h5ad"))
        for f in frag_files
    ]

    # Import unsorted fragments -> backed h5ad files (one per sample)
    # IMPORTANT when passing a list: protect with if __name__ == "__main__" when running as a script.
    snap.pp.import_fragments(
        fragment_file=frag_files,
        chrom_sizes=snap.genome.hg38,      # hg38 / GRCh38
        file=out_files,                     # write backed AnnData files
        min_num_fragments=MIN_FRAGS,
        sorted_by_barcode=False,            # your files are NOT sorted by barcode
        n_jobs=N_JOBS
    )  # returns backed AnnData objects on disk

    # Build an AnnDataSet (one per modality) that links the sample .h5ad files
    keys = [sample_key_from_path(p) for p in frag_files]
    adatas = [(k, snap.read(f)) for k, f in zip(keys, out_files)]
    ads = snap.AnnDataSet(adatas=adatas, filename=f"{dataset_name}.h5ads")
    # Optional: tag the modality
    
    return ads

if __name__ == "__main__":
    # H3K27ac
    ac_ads = import_modality(
        pattern="dna/*_H3K27ac.fragments.tsv.gz",
        out_dir="h5ad/H3K27ac"    )

    # H3K27me3
    me3_ads = import_modality(
        pattern="dna/*_H3K27me3.fragments.tsv.gz",
        out_dir="h5ad/H3K27me3"
    )

    # You now have:
    #   - per-sample backed .h5ad files in h5ad/H3K27ac and h5ad/H3K27me3
    #   - one AnnDataSet per modality: H3K27ac.h5ads and H3K27me3.h5ads

    # (Optional) If you truly need a single AnnData in memory:
big_ac  = ac_ads.to_adata()   # may use substantial RAM
dataset_name="H3K27ac"
big_ac.obs["modality"] = dataset_name
big_ac.obs["barcode_short"] = big_ac.obs_names.str.replace(r":S\d+$", "", regex=True)

big_me3 = me3_ads.to_adata()
dataset_name="H3K27me3"
big_me3.obs["modality"] = dataset_name
big_me3.obs["barcode_short"] = big_me3.obs_names.str.replace(r":S\d+$", "", regex=True)



barcodes = set(df.iloc[:, 0].dropna().astype(str))

big_ac_sub  = big_ac[big_ac.obs["barcode_short"].isin(barcodes) ].copy()
big_me3_sub = big_me3[big_me3.obs["barcode_short"].isin(barcodes) ].copy()




import pandas as pd
df = pd.read_csv("barcode_short.txt",header=None)
df.iloc[:, 0] = df.iloc[:, 0].str.split(":").str[:-1].str.join(":")

import os
import snapatac2 as snap

def rebuild_ads(path, out_name):
    # path points to the folder containing sample .h5ad files
    ad_files = sorted([os.path.join(path, f) for f in os.listdir(path) if f.endswith(".h5ad")])
    # extract keys like GSM8529544_S1 from filenames
    keys = [os.path.splitext(os.path.basename(f))[0] for f in ad_files]
    # create AnnDataSet object and save to .h5ads
    ads = snap.AnnDataSet(
        adatas=[(k, f) for k, f in zip(keys, ad_files)],
        filename=out_name,
    )
    return ads

# H3K27ac
ac_ads = rebuild_ads("subset_ads/H3K27ac/anndatas", "H3K27ac.h5ads")

# H3K27me3
me3_ads = rebuild_ads("subset_ads/H3K27me3/anndatas", "H3K27me3.h5ads")


import os, re, gzip, numpy as np, pandas as pd, snapatac2 as snap

# 2) barcodes to keep (short form)
keep = set(df.iloc[:, 0].dropna().astype(str))

def subset_by_short(ads, out_path):
    # ensure a fresh/empty destination
    if os.path.exists(out_path):
        shutil.rmtree(out_path)
    os.makedirs(out_path, exist_ok=True)

    names = pd.Index(ads.obs_names)
    short = names.str.replace(r":S\d+$", "", regex=True)
    mask  = np.fromiter((s in keep for s in short), dtype=bool)
    idx   = np.flatnonzero(mask)

    # write subset to NEW location
    sub_ads, _ = ads.subset(obs_indices=idx, out=out_path)
    return sub_ads

# 3) choose NEW output dirs (NOT your source dirs)
ac_ads_sub  = subset_by_short(ac_ads,  "subset_ads_filtered/H3K27ac")
me3_ads_sub = subset_by_short(me3_ads, "subset_ads_filtered/H3K27me3")


# per-cell fragment counts -> CSV (barcode_short + n_fragment)
def write_counts_from_ads(ads, outfile):
    obs = ads.obs
    obs["barcode_short"] = pd.Index(ads.obs_names).str.replace(r":S\d+$", "", regex=True)
    obs[["barcode_short", "n_fragment"]].to_csv(outfile, index=False)

os.makedirs("fragments_out", exist_ok=True)
#write_counts_from_ads(ac_ads_sub,  "fragments_out/H3K27ac_n_fragment.csv")
#write_counts_from_ads(me3_ads_sub, "fragments_out/H3K27me3_n_fragment.csv")

# export fragments per-sample, then concatenate into one gz per modality
def export_concat(ads_sub, prefix, out_gz):
    # per-sample exports
    out_dir = f"fragments_out/{prefix}_tmp"
    os.makedirs(out_dir, exist_ok=True)
    files = snap.ex.export_fragments(
        ads_sub, groupby="sample", out_dir=out_dir, prefix=f"{prefix}.", suffix=".bed.gz"
    )
    # concatenate all .bed.gz into a single gzip (single member)
    with gzip.open(out_gz, "wb") as w:
        for fn in sorted(files.values()):
            with gzip.open(os.path.join(out_dir, os.path.basename(fn)), "rb") as r:
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    w.write(chunk)

export_concat(ac_ads_sub,  "H3K27ac",  "fragments_out/H3K27ac.fragments.bed.gz")
export_concat(me3_ads_sub, "H3K27me3", "fragments_out/H3K27me3.fragments.bed.gz")



import scanpy as sc

a = sc.read_h5ad("scRNA_D0_withCellNames.h5ad")
b = sc.read_h5ad("../TrES/Technical/PostQC/H2RRNA_Human_PostQC.h5ad")




import numpy as np
import scanpy as sc

# --- 0) Make sure genes are in var_names and unique ---
if 'features' in a.var.columns and not a.var_names.equals(a.var['features'].astype(str)):
    a.var_names = a.var['features'].astype(str).values
a.var_names_make_unique()

# If raw counts live in a layer (common after conversions), use that for QC
qc_layer = 'counts' if 'counts' in getattr(a, 'layers', {}) else None

# --- 1) Define gene categories in .var (case-insensitive) ---
names_upper = a.var_names.str.upper()

# Mitochondrial genes (e.g., MT-CO1; works for human/mouse as we upper-case)
a.var['mt'] = names_upper.str.startswith(('MT-', 'MT.'))

# Ribosomal protein genes (cytosolic; RPS*, RPL*)
a.var['ribo'] = names_upper.str.match(r'^(RPS|RPL)\d+')

# Hemoglobin genes (HBA, HBB, HBD, HBE, HBG, HBM, HBZ, etc.)
a.var['hb'] = names_upper.str.match(r'^HB[ABDEGMZ]')

# --- 2) Compute Scanpy-style QC metrics ---
sc.pp.calculate_qc_metrics(
    a,
    qc_vars=['mt', 'ribo', 'hb'],
    percent_top=[50, 100, 200, 500],
    layer=qc_layer,
    inplace=True
)

a.obs.index = a.obs.index.str.split('_').str[1]

a.write("scRNA_D0_withCellNames_likeScanpy.h5ad")












###################
###################
###################
#MTC HELPER RNA#### -> in R !!!
###################
###################
###################




# --- Minimal deps ---
library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(stringr)
library(Seurat)
# For h5ad export
suppressPackageStartupMessages(library(SeuratDisk))
# --- I/O root ---
in.file.dir <- "/mnt/dataFast/ahrmad/ScH2R_202409/FullBench/MTC"

# --- find sample directories under ".../mtx" (same idea, just robust) ---
dir <- list.dirs(file.path(in.file.dir, "rna"), full.names = TRUE, recursive = FALSE)
stopifnot(length(dir) > 0)

# keep your naming convention line (no behavioral change)
names(dir) <- str_split(basename(dir), "_", simplify = TRUE)[, 1]

# --- read each 10x matrix directory and create Seurat objects (same params) ---
scRNAlist <- vector("list", length(dir))
for (i in seq_along(dir)) {
  message("Reading: ", dir[i])
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features = 100)
}

# --- merge all runs (same as your code, just generalized to N samples) ---
scRNA2 <- Reduce(function(x, y) merge(x, y), scRNAlist)

# --- QC subset (same thresholds) ---
#scRNA2 <- subset(scRNA2, subset = nCount_RNA >= 1000 & nFeature_RNA >= 500)

# ---- save merged FIRST (RDS + h5ad), exactly as you requested ----
saveRDS(scRNA2, file = file.path(in.file.dir, "scRNA_merged.seurat.rds"))
h5s_all  <- file.path(in.file.dir, "scRNA_merged.h5Seurat")
h5ad_all <- file.path(in.file.dir, "scRNA_merged.h5ad")
SaveH5Seurat(scRNA2, filename = h5s_all, overwrite = TRUE)
Convert(h5s_all, dest = "h5ad", filename = h5ad_all, overwrite = TRUE)
# --- keep your downstream strategy exactly the same (D0 only) ---
rm(scRNAlist)

scRNA2@meta.data$sample_bio.id <- stringr::str_split(rownames(scRNA2@meta.data), ":", simplify = TRUE)[, 4]
scRNA2@meta.data$batch <- "Batch1"

# keep only D0 (code "01")
scRNA2 <- subset(scRNA2, subset = sample_bio.id == "01")

foo <- scRNA2@meta.data %>% tibble::rownames_to_column("CB") %>% as.data.table()
foo[sample_bio.id %in% c("01"), anno := "Def.endo.D0"]
foo[sample_bio.id %in% c("02"), anno := "Def.endo.D1"]
foo[sample_bio.id %in% c("03"), anno := "Def.endo.D2"]
foo[sample_bio.id %in% c("04"), anno := "Def.endo.D3"]
foo <- foo[, .(CB, anno)] %>% as.data.frame() %>%
  tibble::column_to_rownames("CB") %>% .[rownames(scRNA2@meta.data), ]

scRNA2@meta.data <- cbind(scRNA2@meta.data, foo)

# put the human-readable label into sample_bio.id
scRNA2@meta.data$sample_bio.id <- scRNA2@meta.data$anno
# (optionally drop anno)
# scRNA2@meta.data$anno <- NULL

scRNA2@meta.data$CB <- rownames(scRNA2@meta.data)
scRNA2@meta.data$sub.lib <- scRNA2$orig.ident
scRNA2[["percent.mt"]] <- PercentageFeatureSet(scRNA2, pattern = "^MT-")

# save D0-only object (RDS + h5ad)
saveRDS(scRNA2, file = file.path(in.file.dir, "scRNA_D0.seurat.rds"))
h5s_d0  <- file.path(in.file.dir, "scRNA_D0.h5Seurat")
h5ad_d0 <- file.path(in.file.dir, "scRNA_D0.h5ad")
SeuratDisk::SaveH5Seurat(scRNA2, filename = h5s_d0, overwrite = TRUE)
SeuratDisk::Convert(h5s_d0, dest = "h5ad", filename = h5ad_d0, overwrite = TRUE)


mm install -y r-base r-seurat r-seuratdisk r-stringr r-data.table r-purrr r-ggplot2 r-viridis r-tibble r-reticulate r-hdf5r python=3.11 anndata h5py numpy pandas

#!/usr/bin/env bash
# organize_10x_files.sh
# Run this inside your "rna" directory

for prefix in $(ls *.matrix.mtx.gz | sed 's/\.matrix\.mtx\.gz//' ); do
    echo "Organizing $prefix ..."
    mkdir -p "$prefix"
    mv "${prefix}".matrix.mtx.gz "$prefix/matrix.mtx.gz"
    mv "${prefix}".barcodes.tsv.gz "$prefix/barcodes.tsv.gz"
    mv "${prefix}".features.tsv.gz "$prefix/features.tsv.gz"
done


#################################################
#################################################
#################################################
#In-house nanocnt 
#################################################
#################################################
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

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
if not os.path.exists('figures'):
    os.makedirs('figures')

workdir = "/mnt/dataFast/ahrmad/ScH2R_202409/FullBench/InHouse_NanoCNT/"

#################################################
#################### DNA QC #####################
#################################################

print('IMPORT AND QC PLOTS DNA')
maxFrags = 80000

# --- H3K9me3 ---
bam_k9me3 = 'ScKDMA_H3K27ac.bam'
exp_k9me3 = bam_k9me3.split('.bam')[0]
frags_k9me3_tsv = exp_k9me3 + '.tsv.gz'
ad_k9me3_path = exp_k9me3 + '.h5ad'

# --- H3K27me3 ---
bam_k27me3 = 'ScKDMA_H3K27me3.bam'
exp_k27me3 = bam_k27me3.split('.bam')[0]
frags_k27me3_tsv = exp_k27me3 + '.tsv.gz'
ad_k27me3_path = exp_k27me3 + '.h5ad'

# Create fragment files (paired-end, 20 MAPQ cutoff, with shifts)
k9me3_info = snap.pp.make_fragment_file(
    workdir + bam_k9me3,
    workdir + frags_k9me3_tsv,
    is_paired=True,
    barcode_tag='RG',
    shift_left=4,
    shift_right=-5,
    min_mapq=20,
    chunk_size=500000000,
    chrM=['chrM', 'M'],
    tempdir=workdir,
)

with open(f'{exp_k9me3}.csv', mode='w', newline='') as f:
    writer = csv.writer(f)
    for key, value in k9me3_info.items():
        writer.writerow([key, value])

k27me3_info = snap.pp.make_fragment_file(
    workdir + bam_k27me3,
    workdir + frags_k27me3_tsv,
    is_paired=True,
    barcode_tag='RG',
    shift_left=4,
    shift_right=-5,
    min_mapq=20,
    chunk_size=500000000,
    chrM=['chrM', 'M'],
    tempdir=workdir,
)

with open(f'{exp_k27me3}.csv', mode='w', newline='') as f:
    writer = csv.writer(f)
    for key, value in k27me3_info.items():
        writer.writerow([key, value])

# Import fragments as AnnData
ad_k9me3 = snap.pp.import_fragments(
    workdir + frags_k9me3_tsv,
    snap.genome.hg38,
    min_num_fragments=0,
    sorted_by_barcode=True,
    chrM=['chrM', 'M'],
    shift_left=0,
    shift_right=0,
    chunk_size=50000,
    tempdir=workdir,
    backend='hdf5',
    n_jobs=64,
)

ad_k27me3 = snap.pp.import_fragments(
    workdir + frags_k27me3_tsv,
    snap.genome.hg38,
    min_num_fragments=0,
    sorted_by_barcode=True,
    chrM=['chrM', 'M'],
    shift_left=0,
    shift_right=0,
    chunk_size=50000,
    tempdir=workdir,
    backend='hdf5',
    n_jobs=64,
)



# Tag and export
ad_k9me3.obs['Exp'] = [exp_k9me3] * ad_k9me3.n_obs
ad_k27me3.obs['Exp'] = [exp_k27me3] * ad_k27me3.n_obs

snap.ex.export_fragments(ad_k9me3, groupby='Exp', prefix='', suffix='.tsv.gz')
snap.ex.export_fragments(ad_k27me3, groupby='Exp', prefix='', suffix='.tsv.gz')

