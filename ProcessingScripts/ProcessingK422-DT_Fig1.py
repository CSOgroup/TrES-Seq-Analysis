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

workdir = "/mnt/dataFast/ahrmad/triseq_202510/processed/"

#################################################
#################### DNA QC #####################
#################################################

print('IMPORT AND QC PLOTS DNA')
maxFrags = 80000

# --- H3K9me3 ---
bam_k9me3 = 'Sc_K9D_H3K9me3_NoDup.bam'
exp_k9me3 = bam_k9me3.split('_NoDup.bam')[0]
frags_k9me3_tsv = exp_k9me3 + '.tsv.gz'
ad_k9me3_path = exp_k9me3 + '.h5ad'

# --- H3K27me3 ---
bam_k27me3 = 'Sc_K9D_H3K27me3_NoDup.bam'
exp_k27me3 = bam_k27me3.split('_NoDup.bam')[0]
frags_k27me3_tsv = exp_k27me3 + '.tsv.gz'
ad_k27me3_path = exp_k27me3 + '.h5ad'

# Create fragment files (paired-end, 20 MAPQ cutoff, with shifts)
k9me3_info = snap.pp.make_fragment_file(
    workdir + bam_k9me3,
    workdir + frags_k9me3_tsv,
    is_paired=True,
    barcode_tag='CB',
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
    barcode_tag='CB',
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

# Metrics
snap.metrics.tsse(ad_k9me3, gene_anno=snap.genome.hg38, inplace=True)
snap.metrics.tsse(ad_k27me3, gene_anno=snap.genome.hg38, inplace=True)

snap.pp.filter_cells(ad_k9me3, min_counts=100, max_counts=maxFrags, min_tsse=0, inplace=True)
snap.pp.filter_cells(ad_k27me3, min_counts=100, max_counts=maxFrags, min_tsse=0, inplace=True)

snap.metrics.frag_size_distr(ad_k9me3, add_key='frag_size_distr', inplace=True)
snap.metrics.frag_size_distr(ad_k27me3, add_key='frag_size_distr', inplace=True)

ad_k9me3.obs['log10_n_fragment'] = np.log10(ad_k9me3.obs['n_fragment'])
ad_k27me3.obs['log10_n_fragment'] = np.log10(ad_k27me3.obs['n_fragment'])

# Violin for H3K9me3
sc.pl.violin(ad_k9me3, 'log10_n_fragment', rotation=1e-7, show=False)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {ad_k9me3.n_obs}\nMedian n_fragment:{np.median(ad_k9me3.obs['n_fragment'])}")
plt.savefig(f'./figures/{exp_k9me3}_FragViolin.pdf'); plt.close()

# Violin for H3K27me3
sc.pl.violin(ad_k27me3, 'log10_n_fragment', rotation=1e-7, show=False)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {ad_k27me3.n_obs}\nMedian n_fragment:{np.median(ad_k27me3.obs['n_fragment'])}")
plt.savefig(f'./figures/{exp_k27me3}_FragViolin.pdf'); plt.close()

# Scatter log10_n_fragment vs tsse
sc.pl.scatter(ad_k9me3,  "log10_n_fragment", "tsse", title=exp_k9me3,  save=exp_k9me3 + '_tsse.pdf')
sc.pl.scatter(ad_k27me3, "log10_n_fragment", "tsse", title=exp_k27me3, save=exp_k27me3 + '_tsse.pdf')

# Density/frag-size plots
snap.pl.tsse(ad_k9me3, min_fragment=0, width=750, height=600, interactive=False, show=False,
             out_file=f"./figures/{exp_k9me3}_tsseDensity_PreFilter.pdf")
snap.pl.frag_size_distr(ad_k9me3, width=750, height=600, interactive=False, show=False,
                        out_file=f"./figures/{exp_k9me3}_frag_size_distr_PreFilter.pdf")

snap.pl.tsse(ad_k27me3, min_fragment=0, width=750, height=600, interactive=False, show=False,
             out_file=f"./figures/{exp_k27me3}_tsseDensity_PreFilter.pdf")
snap.pl.frag_size_distr(ad_k27me3, width=750, height=600, interactive=False, show=False,
                        out_file=f"./figures/{exp_k27me3}_frag_size_distr_PreFilter.pdf")

# Tag and export
ad_k9me3.obs['Exp'] = [exp_k9me3] * ad_k9me3.n_obs
ad_k27me3.obs['Exp'] = [exp_k27me3] * ad_k27me3.n_obs

snap.ex.export_fragments(ad_k9me3, groupby='Exp', prefix='', suffix='.tsv.gz')
ad_k9me3.write(filename=ad_k9me3_path)
snap.ex.export_fragments(ad_k27me3, groupby='Exp', prefix='', suffix='.tsv.gz')
ad_k27me3.write(filename=ad_k27me3_path)

#################################################
#################### RNA QC #####################
#################################################

print('IMPORT AND QC PLOTS RNA')

cnt_tsv = 'Sc_K9r_NoDupUMI_count.tsv'
sample_rna = 'Sc_K9r'

adata = sc.read_csv(workdir + cnt_tsv, delimiter='\t', first_column_names=True)
adata = adata.T

# mitochondrial genes (human), ribosomal, and hemoglobin
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.str.contains('^HB[^ (P)]')  # avoid HBP/ HBPx

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True, log1p=True)

sc.pl.violin(
    adata,
    ['pct_counts_mt', 'pct_counts_in_top_50_genes'],
    jitter=0.4,
    multi_panel=True,
    rotation=1e-7,
    save=f'{sample_rna}_Pct_PreFilter.pdf',
)

# Violin (log scale shown, then relabel ticks back to raw)
sc.pl.violin(
    adata,
    ['log1p_n_genes_by_counts', 'log1p_total_counts'],
    jitter=0.4,
    multi_panel=True,
    rotation=1e-7,
    show=False,
)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
plt.savefig(f"./figures/{sample_rna}_PreFilter.pdf"); plt.close()

sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_mt',
              title=sample_rna + '\nColor: pct_counts_mt', save=f'{sample_rna}_PreFilter.pdf')
sc.pl.scatter(adata, 'log1p_total_counts', 'log1p_n_genes_by_counts', color='pct_counts_mt',
              title=sample_rna + '\nColor: pct_counts_mt', save=f'{sample_rna}_PreFilterlog.pdf')

adata.write(f'{sample_rna}.h5ad')


#################################################
############# RNA + DNA FILTERING ###############
#################################################

exp_base = 'Sc_K9'

exp_rna = 'Sc_K9r'
exp_k9me3 = 'Sc_K9D_H3K9me3'
exp_k27me3 = 'Sc_K9D_H3K27me3'

# Load AnnData
rna_ad = sc.read_h5ad(f'{exp_rna}.h5ad')
k9me3_ad = sc.read_h5ad(f'{exp_k9me3}.h5ad')
k27me3_ad = sc.read_h5ad(f'{exp_k27me3}.h5ad')

# Mapping tables
IDX_HIST = {'GAT': '1', 'AGT': '2', 'TCA': '3', 'ACG': '4', 'GTC': '5'}
IDX_RNA  = {'GCAT': '1', 'CGAT': '2', 'GACT': '3', 'AGCT': '4', 'CAGT': '5'}

INLINE_OFFSET = 1      # starts 1 nt from beginning (0-based index 1)
LEN_HIST = 3
LEN_RNA  = 4

def harmonize_obs_names(adata, assay, inline_offset=INLINE_OFFSET):
    """
    Create harmonized cell IDs from inline indices embedded in the barcode.
    New obs_names format: "{group}_{core_barcode_without_inline_tag}"

    Columns added:
      - barcode_raw: original obs_names
      - inline_index_seq: extracted inline tag sequence
      - inline_group: '1'..'5' as string
      - barcode_core: barcode with inline tag removed
      - barcode_harmonized: final "{group}_{core}" (also set as obs_names)
    """
    if assay in ('k9me3', 'k27me3'):
        tag_len = LEN_HIST
        mapping = IDX_HIST
    elif assay == 'rna':
        tag_len = LEN_RNA
        mapping = IDX_RNA
    else:
        raise ValueError("assay must be one of: 'k9me3', 'k27me3', 'rna'")

    # Work as a Series for vectorized slicing/mapping
    s = pd.Series(adata.obs_names, index=adata.obs_names, dtype=object)

    # Extract inline index and map to group
    inline_seq = s.str.slice(inline_offset, inline_offset + tag_len)
    group = inline_seq.map(mapping)

    # Validate all tags are recognized
    if group.isna().any():
        bad = s[group.isna()]
        bad_tag = inline_seq[group.isna()]
        example = f"{bad.iloc[0]} (tag='{bad_tag.iloc[0]}')" if len(bad) else "n/a"
        raise ValueError(
            f"Found {group.isna().sum()} barcodes with unknown inline tag. "
            f"First example: {example}"
        )

    # Remove the inline tag from the barcode to form the core
    core = s.str.slice(0, inline_offset) + s.str.slice(inline_offset + tag_len, None)

    # New harmonized ID (prefix with group to preserve uniqueness across groups)
    new_ids = group + "_" + core

    # Ensure uniqueness; if collisions exist, append numeric suffix
    if new_ids.duplicated().any():
        dup_idx = new_ids.duplicated(keep=False)
        # add _1, _2... to duplicates (stable)
        new_ids.loc[dup_idx] = (
            new_ids[dup_idx] + "_" +
            new_ids[dup_idx].groupby(new_ids[dup_idx]).cumcount().astype(str)
        )

    # Attach columns and set obs_names
    adata.obs['barcode_raw'] = s.values
    adata.obs['inline_index_seq'] = inline_seq.values
    adata.obs['inline_group'] = group.values
    adata.obs['barcode_core'] = core.values
    adata.obs['barcode_harmonized'] = new_ids.values
    adata.obs_names = pd.Index(new_ids.astype(str).to_numpy(dtype=object),name=adata.obs_names.name)

    return adata

# --- Apply to your three AnnData objects (in-place) ---
# k27me3_ad, k9me3_ad, rna_ad already exist as in your session
harmonize_obs_names(k27me3_ad, assay='k27me3')
harmonize_obs_names(k9me3_ad,  assay='k9me3')
harmonize_obs_names(rna_ad,    assay='rna')



# DNA thresholds
minFrags_k27me3 = 400
maxFrags_k27me3 = 6000
minFrags_k9me3 = 400
maxFrags_k9me3 = 6000

# RNA thresholds
minUMI = 500
maxUMI = 8000
max_pct_mt = 3
max_pct_hb = 5
max_pct_50g = 25

# --- RNA filtering criteria ---
filtered_rna = rna_ad[rna_ad.obs['total_counts'] >= minUMI]
filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_hb'] < max_pct_hb]

# --- DNA filtering criteria ---
filtered_k9me3 = k9me3_ad[snap.pp.filter_cells(k9me3_ad, min_counts=minFrags_k9me3,
                                              max_counts=maxFrags_k9me3, min_tsse=0, inplace=False)]
filtered_k27me3 = k27me3_ad[snap.pp.filter_cells(k27me3_ad, min_counts=minFrags_k27me3,
                                                max_counts=maxFrags_k27me3, min_tsse=0, inplace=False)]

# Barcode sets
rna_barcodes = set(filtered_rna.obs.index)
k9me3_barcodes = set(filtered_k9me3.obs.index)
k27me3_barcodes = set(filtered_k27me3.obs.index)

shared_barcodes = rna_barcodes & k9me3_barcodes & k27me3_barcodes

print(f'K9me3:{len(k9me3_barcodes)}\nK27me3:{len(k27me3_barcodes)}\nRNA:{len(rna_barcodes)}\nShared:{len(shared_barcodes)}')

# Venn
plt.figure(figsize=(8, 8))
venn = venn3(
    [rna_barcodes, k9me3_barcodes, k27me3_barcodes],
    set_labels=(
        f'RNA\n{minUMI}<UMIs<{maxUMI}\n<{max_pct_mt}% MT\n<{max_pct_50g}% 50TopGenes',
        f'H3K9me3\n{minFrags_k9me3}<fragments<{maxFrags_k9me3}',
        f'H3K27me3\n{minFrags_k27me3}<fragments<{maxFrags_k27me3}',
    ),
)
plt.suptitle('Venn Diagram of Passing Cells\nfor RNA, H3K9me3, and H3K27me3', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f'./figures/{exp_base}_VennPassingCells.pdf')
plt.close()

# UpSet
all_items = rna_barcodes.union(k9me3_barcodes).union(k27me3_barcodes)

data = {
    'RNA': [item in rna_barcodes for item in all_items],
    'H3K9me3': [item in k9me3_barcodes for item in all_items],
    'H3K27me3': [item in k27me3_barcodes for item in all_items],
}

df = pd.DataFrame(data, index=list(all_items))

fig = plt.figure(figsize=(8, 10))
plot(
    from_indicators(df.columns, data=df),
    sort_by='cardinality',
    show_percentages=True,
    min_subset_size="0.2%",
    facecolor="darkblue",
    fig=fig,
    element_size=None,
)
plt.title('Passing Cells\nfor RNA, H3K9me3, and H3K27me3', fontsize=16, fontweight='bold')
plt.savefig(f'./figures/{exp_base}_UpsetPassingCells.pdf')
plt.close()

# Subset to shared barcodes
filtered_k9me3 = filtered_k9me3[filtered_k9me3.obs.index.isin(shared_barcodes)]
filtered_k27me3 = filtered_k27me3[filtered_k27me3.obs.index.isin(shared_barcodes)]
filtered_rna = filtered_rna[filtered_rna.obs.index.isin(shared_barcodes)]

print(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.X.sum(axis=1))}')

for adata_obj, tag in [(filtered_rna, 'rna'), (filtered_k9me3, 'ac'), (filtered_k27me3, 'me3')]:
    cell_names = adata_obj.obs['barcode_raw']
    with open(f"{exp_base}_{tag}_CountFiltered_cell_names.lst", "w") as f:
        f.write("\n".join(cell_names))

# Merge obs for plotting
merged_ad = filtered_rna.copy()
merged_ad.obs = merged_ad.obs.reindex(index=sorted(merged_ad.obs.index))

k9me3_c = filtered_k9me3.copy(); k9me3_c.obs = k9me3_c.obs.reindex(index=sorted(k9me3_c.obs.index))
k27me3_c = filtered_k27me3.copy(); k27me3_c.obs = k27me3_c.obs.reindex(index=sorted(k27me3_c.obs.index))

assert k27me3_c.obs_names.equals(k9me3_c.obs_names) and k27me3_c.obs_names.equals(merged_ad.obs_names), \
    "The obs_names do not match across k27me3_c, k9me3_c, and merged_ad"

k9me3_obs = k9me3_c.obs.add_suffix('_ac')      # keep original suffix to avoid changing filenames downstream
k27me3_obs = k27me3_c.obs.add_suffix('_me3')
merged_ad.obs = merged_ad.obs.join([k9me3_obs, k27me3_obs])

# Scatter of log10 fragments across marks vs RNA UMIs
sc.pl.scatter(merged_ad, 'log10_n_fragment_me3', 'log10_n_fragment_ac', color='log1p_total_counts',
              title='Log10 Unique Fragments in Histone Marks\nColored by log1p RNA UMIs', show=False)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks(); ax.set_yticklabels([f"{10**tick - 1:.0f}" for tick in yticks])
for ax in axes:
    xticks = ax.get_xticks(); ax.set_xticklabels([f"{10**tick - 1:.0f}" for tick in xticks])
plt.savefig(f"./figures/{exp_base}_Counts3.png"); plt.close()

# RNA post-filter plots
sc.pl.violin(
    filtered_rna,
    ['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt', 'pct_counts_in_top_50_genes'],
    jitter=0.4,
    multi_panel=True,
    rotation=1e-7,
    save=f'{exp_base}_rna_PostFilter.png',
)

sc.pl.violin(
    filtered_rna,
    ['log1p_n_genes_by_counts', 'log1p_total_counts'],
    jitter=0.4,
    multi_panel=True,
    rotation=1e-7,
    show=False,
)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks()
    ax.set_yticklabels([f"{np.expm1(tick):.0f}" for tick in yticks])
plt.title(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.X.sum(axis=1))}')
plt.savefig(f"./figures/{exp_base}_rna_PostFilterNoLog.pdf"); plt.close()

sc.pl.scatter(filtered_rna, 'total_counts', 'n_genes_by_counts', color='pct_counts_mt',
              title=exp_base + '\nColor: pct_counts_mt', save=f'{exp_base}_rna_PostFilter.png')
sc.pl.scatter(filtered_rna, 'log1p_total_counts', 'log1p_n_genes_by_counts', color='pct_counts_mt',
              title=exp_base + '\nColor: pct_counts_mt', save=f'{exp_base}_rna_PostFilterlog.png')

# Post-filter DNA violins (with relabeled y axis)
sc.pl.violin(filtered_k27me3, 'log10_n_fragment', rotation=1e-7, show=False)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks()
    ax.set_yticklabels([f"{10**tick - 1:.0f}" for tick in yticks])
    ax.set_title(f"n_cells: {filtered_k27me3.n_obs}\nMedian n_fragment: {np.median(filtered_k27me3.obs['n_fragment'])}")
plt.savefig(f'./figures/{exp_base}_me3_nfragPostFilter.png'); plt.close()

sc.pl.violin(filtered_k9me3, 'log10_n_fragment', rotation=1e-7, show=False)
fig = plt.gcf(); axes = fig.axes
for ax in axes:
    yticks = ax.get_yticks()
    ax.set_yticklabels([f"{10**tick - 1:.0f}" for tick in yticks])
    ax.set_title(f"n_cells: {filtered_k9me3.n_obs}\nMedian n_fragment: {np.median(filtered_k9me3.obs['n_fragment'])}")
plt.savefig(f'./figures/{exp_base}_ac_nfragPostFilter.png'); plt.close()

# Save filtered AnnData (filenames kept unchanged for downstream compatibility)
filtered_rna.write(f'{exp_base}_rna_CountFiltered.h5ad')
filtered_k9me3.write(f'{exp_base}_k9_CountFiltered.h5ad')
filtered_k27me3.write(f'{exp_base}_me3_CountFiltered.h5ad')

# -----------------------------------------------------------------------------
# Per-sample violin helper (no logic change to calculations)
# -----------------------------------------------------------------------------
import seaborn as sns

def _backtransform(vals, modality):
    """Convert plotted values back to raw counts."""
    if 'H3' in modality:
        # assuming log10(n + 1)
        return (10 ** vals) - 1
    else:
        # assuming natural log1p
        return np.expm1(vals)


def vio(adata_in, modality):
    """
    Violin plot per sample for a given AnnData and modality.
    Title is just the modality. Shows per-sample median (raw counts) as a dot + label.
    """
    # choose count key
    if 'H3' in modality:
        count_key = 'log10_n_fragment'
        y_label = 'Fragments per cell'
    else:
        count_key = 'log1p_total_counts'
        y_label = 'UMI counts per cell'

    if count_key not in adata_in.obs.columns:
        raise KeyError(f"Expected '{count_key}' in adata.obs for modality '{modality}'")

    # (Retain original placeholder sample assignment)
    adata_in.obs['sample'] = 'K422_K9_DT'

    dfv = adata_in.obs[['sample', count_key]].copy()
    dfv['sample'] = dfv['sample'].astype(str)
    dfv.rename(columns={count_key: 'metric'}, inplace=True)

    sample_counts = dfv['sample'].value_counts().sort_index()
    dfv['sample_label'] = dfv['sample'].map(lambda s: f"{s}\n(n={sample_counts[s]})")
    cats = sorted(dfv['sample_label'].unique(), key=lambda x: (x.split('\n')[0]))
    dfv['sample_label'] = pd.Categorical(dfv['sample_label'], categories=cats, ordered=True)

    width = max(4, 1.3 * dfv['sample_label'].nunique())
    plt.figure(figsize=(width, 4))

    sns.violinplot(data=dfv, x='sample_label', y='metric', inner=None, cut=0, linewidth=0.5, color='lightgray')
    sns.boxplot(data=dfv, x='sample_label', y='metric', whis=1.5, width=0.18, fliersize=0, boxprops=dict(alpha=0.6))

    ax = plt.gca()

    med_plot = dfv.groupby('sample_label')['metric'].median()
    med_raw = _backtransform(med_plot.values, modality)
    for i, (cat, y_med_plot, y_med_raw) in enumerate(zip(med_plot.index, med_plot.values, med_raw)):
        ax.plot(i, y_med_plot, marker='o', markersize=5, mec='black', mfc='white', zorder=5)
        ax.text(i, y_med_plot, f"{int(round(y_med_raw)):,}", ha='center', va='bottom', fontsize=8, rotation=0, zorder=6)

    yticks = ax.get_yticks()
    if 'H3' in modality:
        ax.set_yticklabels([f"{int(max(0, round(10**y - 1))):,}" for y in yticks])
    else:
        ax.set_yticklabels([f"{int(max(0, round(np.expm1(y)))):,}" for y in yticks])

    plt.xticks(rotation=45, ha='right')
    plt.xlabel('')
    plt.ylabel(y_label)
    plt.title(modality)
    plt.tight_layout()
    outpath = f'./figures/Type_VIOLIN_{modality}.pdf'
    plt.savefig(outpath)
    plt.close()
    print(f"Saved: {outpath}")

# Calls for the three subsetted AnnData objects (already subset to shared_barcodes above)
vio(filtered_rna, 'RNA')
vio(filtered_k9me3, 'H3K9me3')
vio(filtered_k27me3, 'H3K27me3')
