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


#################################################
############# RNA ANALYSIS ###############
#################################################

cell_lines = [
    "Sc_MWD_WSU_CD47_Inhibitor",
    "Sc_MWD_WSU_Untreated",
]


# List to hold AnnData objects for each sample
adatas = []
sample_names = []

#For R code
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
#ro.pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython



##Doublet Analysis
%R library(Seurat)
%R library(scater)
%R library(scDblFinder)
%R library(BiocParallel)
%R library(scran)
%R library(scry)


for sample in cell_lines:
    exp = sample
    try:
        adata_sample = sc.read_h5ad(f"{exp}_rna_CountFiltered.h5ad")
        adata_sample.layers["counts"] = adata_sample.X.copy()
        # Normalize and log1p
        ### SCTRANSFORM Normalization
        print('SCTransform Normalization')
        from rpy2.robjects.packages import importr
        from rpy2.robjects import r, pandas2ri
        anndata2ri.activate()
        #pandas2ri.activate()

        sc.pp.filter_genes(adata_sample,min_cells=5,inplace=True)
        mat = adata_sample.X.copy()

        # Set names for the input matrix
        cell_names = adata_sample.obs_names
        gene_names = adata_sample.var_names
        r.assign('mat', mat.T)
        r.assign('cell_names', cell_names)
        r.assign('gene_names', gene_names)
        r('colnames(mat) <- cell_names')
        r('rownames(mat) <- gene_names')

        seurat = importr('Seurat')
        r('seurat_obj <- CreateSeuratObject(mat)')

        # Run
        r(f'seurat_obj <- SCTransform(seurat_obj,vst.flavor="v2")')

        # Extract the SCT data and add it as a new layer in the original anndata object
        sct_data = np.asarray(r['as.matrix'](r('seurat_obj@assays$SCT@data')))
        adata_sample.layers['SCT_data'] = sct_data.T
        sct_data = np.asarray(r['as.matrix'](r('seurat_obj@assays$SCT@counts')))
        adata_sample.layers['SCT_counts'] = sct_data.T

        ### Shifted logarithm Normalization
        print('Shifted logarithm Normalization')
        scales_counts = sc.pp.normalize_total(adata_sample, target_sum=1e4, inplace=False)
        adata_sample.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

        # choose one of the layers
        adata_sample.X = adata_sample.layers['SCT_data']

        adata_sample.write(f"{exp}_rna_CountFiltered_SCT.h5ad")
        
        adatas.append(adata_sample)
        sample_names.append(exp)
    except Exception as e:
        print(f"Could not load {exp}_rna_CountFiltered.h5ad: {e}")


adata = ad.concat(
    adatas,
    join="inner",
    label="sampleMerge",
    keys=sample_names,
    index_unique=None  # Keep original cell barcodes
)


adata.write("Cd47i_Untreated_FULLrna_merged_SCT.h5ad")


adata = sc.read_h5ad("Cd47i_Untreated_FULLrna_merged_SCT.h5ad")
adata.X = adata.layers['log1p_norm']

# Map samples -> condition
cond_map = {
    "Sc_MWD_WSU_CD47_Inhibitor": "CD47_Inhibitor",
    "Sc_MWD_WSU_Untreated": "Untreated"
}
adata.obs["condition"] = adata.obs["sampleMerge"].map(cond_map).astype("category")


print('Deviance')
xTopDeviantGenes = 2000
ro.globalenv["adata"] = adata

# Feature selection 
%R sce = devianceFeatureSelection(adata, assay="X")

binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T
idx = binomial_deviance.argsort()[-xTopDeviantGenes:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


# Highly variable genes (across all cells)
sc.pp.highly_variable_genes(adata,n_top_genes=2000,flavor='seurat_v3',layer='counts')
sc.pl.highly_variable_genes(adata, show=False, save=f'CD47i_Untreated_HVG_ScanpySeuratV3_onSCT.pdf')


# PCA
sc.tl.pca(adata,n_comps=10)
sc.pl.pca_variance_ratio(adata, n_pcs=10, log=True, show=False, save=f'Cd47i_Untreated_pca_variance_ratio_SCT.pdf')

# PCA plots
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts","log1p_n_genes_by_counts", "condition"],
    dimensions=[(0, 1), (0, 1), (0, 1), (0, 1), (0, 1),],
    ncols=3,
    show=False, save=f'Cd47i_Untreated_PCA1-2_QC_SCT.pdf'
)
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts","log1p_n_genes_by_counts", "condition"],
    dimensions=[(2, 3), (2, 3), (2, 3), (2, 3), (2, 3),],
    ncols=3,
    show=False, save=f'Cd47i_Untreated_PCA3-4_QC_SCT.pdf'
)

# Neighbors/UMAP
sc.pp.neighbors(adata, n_pcs=10, n_neighbors=30)
sc.tl.umap(adata)

# UMAPs colored by QC and sample
sc.pl.umap(
    adata,
    color=[ "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_50_genes", "log1p_total_counts","condition",],
    ncols=3,
    show=False, save=f'Cd47i_Untreated_UMAP_QC_SCT.pdf'
)
sc.pl.umap(
    adata,
    color=["log1p_total_counts", "log1p_n_genes_by_counts", "condition"],
    ncols=3,
    show=False, save=f'Cd47i_Untreated_UMAP_QC_Log_SCT.pdf'
)

sc.tl.leiden(adata, resolution=0.01)
sc.pl.umap(
    adata,
    color=['leiden', 'condition'],
    legend_loc=  "on data",
    ncols=2,
    show=False, save=f'Cd47i_Untreated_UMAP_QC_RNA_Leiden_SCT.pdf'
)




# =========================
# 1) Define marker sets
# =========================
bcell_primary = ["MS4A1","CD19","CD79A","CD79B","PAX5","POU2AF1","SPIB","BCL6"]
bcell_ig      = ["IGHG1","IGHG3","IGHG4","IGKC","IGLC2","IGLC3","IGHA1"]  # secondary (ambient-prone)
mac_core      = ["LYZ","TYROBP","FCER1G","LST1","AIF1","CTSS","CTSB","LAPTM5","C1QA","C1QB","C1QC"]
mac_rec       = ["CSF1R","FCGR1A","FCGR3A","MSR1","CD68","CD14","LILRB1"]

bcell_genes = bcell_primary + bcell_ig
mac_genes   = mac_core + mac_rec


# =========================
# 3) Scores
# =========================
sc.tl.score_genes(adata, bcell_primary, score_name="score_bcell_core", use_raw=False,)
sc.tl.score_genes(adata, mac_core, score_name="score_mac_core",   use_raw=False,)
sc.tl.score_genes(adata, bcell_ig, score_name="score_bcell_Ig", use_raw=False,)
sc.tl.score_genes(adata, mac_rec, score_name="score_mac_rec", use_raw=False,)

adata.obs["score_delta_B_minus_M"] = adata.obs["score_bcell_core"] - adata.obs["score_mac_core"]

def my_vmax(values): return np.quantile(np.abs(values), 0.98)
def my_vmin(values): return -np.quantile(np.abs(values), 0.98)

sc.pl.umap(
    adata,
    color=["score_bcell_core","score_mac_core","score_delta_B_minus_M"],
    color_map="coolwarm",
    vmin=my_vmin, vmax=my_vmax,
    frameon=False, show=False,
    save="_delta_B_minus_M.pdf"
)

# UMAP for key markers
markers_to_show = [g for g in ["MS4A1","CD19","CD79A","LYZ","TYROBP","FCER1G"] if g in adata.var_names]
def my_vmax(values): return 2*np.mean(values)

if markers_to_show:
    sc.pl.umap(
        adata,
        color=markers_to_show,
        ncols=3,
        vmax=my_vmax,
        frameon=False, show=False,
        save="_key_markers.pdf"
    )


# =========================
# 4) Auto thresholds and calls
# =========================
q_hi = adata.obs["score_delta_B_minus_M"].quantile(0.75)
q_lo = adata.obs["score_delta_B_minus_M"].quantile(0.25)

def call_cell(delta, hi=q_hi, lo=q_lo):
    if delta >= hi:
        return "WSU-DLCL2"
    elif delta <= lo:
        return "Macrophage"
    else:
        return ""

adata.obs["Individual_Cell_Call"] = adata.obs["score_delta_B_minus_M"].apply(call_cell).astype("category")

# =========================
# 5) SAVE figures (not showing)
# =========================
# UMAP by discrete label
sc.pl.umap(
    adata,
    color=["Individual_Cell_Call"],
    frameon=False, show=False,
    save="_co_culture_IndividualCellCall.pdf"
)

# Violin distributions (scores & delta)
sc.pl.violin(
    adata,
    keys=["score_delta_B_minus_M"],
    groupby="leiden",
    multi_panel=True,
    show=False, save="_Deltascore_by_leiden.pdf"
)

sc.tl.leiden(adata, key_added="leidenMore", resolution=2)

# Violin distributions (scores & delta)
sc.pl.violin(
    adata,
    keys=["score_delta_B_minus_M"],
    groupby="leidenMore",
    multi_panel=True,
    show=False, save="_Deltascore_by_leidenMore.pdf"
)

sc.pl.umap(adata, color="condition", wspace=0.4, show=False, save=f"Cd47i_Untreated_UMAP_condition.pdf")
sc.pl.umap(adata, color="leidenMore", wspace=0.4,legend_loc='on data', show=False, save=f"Cd47i_Untreated_UMAP_leidenMore.pdf")

def my_vmax(values): return np.quantile(np.abs(values), 0.98)
def my_vmin(values): return -np.quantile(np.abs(values), 0.98)
sc.pl.umap(
    adata,
    color=["score_bcell_core","score_mac_core","score_delta_B_minus_M"],
    color_map="coolwarm",
    vmin=my_vmin, vmax=my_vmax,
    frameon=False, show=False,
    save="_delta_B_minus_M.pdf"
)

adata = adata[~adata.obs["leidenMore"].isin(["26", "27","21","25","28"])].copy()

adata.obs["cluster_call"] = adata.obs["leiden"].map({'0': 'WSU-DLCL2','1': 'Macrophage'})
sc.pl.umap(adata, color=["cluster_call"], show=False, frameon=False, save="_FULLRNA_UMAPForFig_ClusterCall.pdf")
sc.pl.umap(adata, color=["condition"], show=False, frameon=False, save="_FULLRNA_UMAPForFig_Sample.pdf")
sc.pl.umap(
    adata,
    color=["score_bcell_core","score_mac_core","score_delta_B_minus_M"],
    color_map="coolwarm",
    vmin=my_vmin, vmax=my_vmax,
    frameon=False, show=False,
    save="_FULLRNA_UMAPForFig_delta_B_minus_M.pdf"
)

adata.write("Cd47i_Untreated_FULLrna_merged_processedSCT_Called_Clean.h5ad")



# =========================
# 7) Split by cluster_call, re-cluster & re-UMAP per subset, show samples
# =========================

adata = sc.read_h5ad("Cd47i_Untreated_FULLrna_merged_processedSCT_Called_Clean.h5ad")

celltypes = {
    "WSU-DLCL2": (adata.obs["cluster_call"] == "WSU-DLCL2"),
    "Macrophage": (adata.obs["cluster_call"] == "Macrophage"),
}
subsets = {k: adata[v].copy() for k, v in celltypes.items() if v.sum() > 0}

def preprocess_and_embed(ad_sub, name, n_top_genes=3000):
    # Basic reprocessing for clean subset embedding
    sc.pp.highly_variable_genes(ad_sub, n_top_genes=n_top_genes, flavor="seurat_v3")
    sc.tl.pca(ad_sub, n_comps=20)
    sc.pl.pca_variance_ratio(adata, n_pcs=20, log=True, show=False, save=f'SUB{name}_pca_variance_ratio_SCT.pdf')
    sc.pp.neighbors(ad_sub, n_pcs=20, n_neighbors=20)
    sc.tl.leiden(ad_sub, key_added="leiden", resolution=0.5)
    sc.tl.umap(ad_sub)
    sc.pl.umap(ad_sub, color=["condition","leiden"], wspace=0.4,
               frameon=False, show=False, save=f"_{name}_UMAP_condition_leiden.pdf")
    return ad_sub

for name, ad_sub in subsets.items():
    subsets[name] = preprocess_and_embed(ad_sub, name)
    ad_sub.write(f"adata_FULLRNA_{name}.h5ad")




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pseudobulk differential expression for single-cell data using:

- edgeR QLF (glmQLFTest)
- DESeq2 (Wald test via DESeq())

with fallback to effect-sizes only (edgeR) if no replication.

- Forced contrast: CD47_Inhibitor vs Untreated (logFC = ALT - REF).
- Labels placed with adjustText (pip install adjustText).
- Point size ∝ detection (cells with counts > 0 for that gene).
- "Smart" label selection: genes with logCPM>4 and |log2FC|>1, ranked equally by both.
- Volcano subtitles show number of significantly up/down genes.

Outputs in: figures_dge_pseudobulk/
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from adjustText import adjust_text

# Optional: make pseudobulk splitting reproducible
# np.random.seed(0)

# ------------------- Configuration -------------------
OUTDIR = "figures_dge_pseudobulk_log10"
os.makedirs(OUTDIR, exist_ok=True)

adatas = {
    "WSU-DLCL2": sc.read_h5ad("adata_FULLRNA_WSU-DLCL2.h5ad"),
    "Macrophage": sc.read_h5ad("adata_FULLRNA_Macrophage.h5ad"),
}

CONDITION_KEY = "condition"    # Must exist in .obs
SAMPLE_KEY    = None           # Set to a column for true replicates; keep None if none
CELLTYPE_KEY  = None           # Optional: per-cell-type DE (e.g., "cluster_call")
REPLICATE_KEY = None           # Optional: blocking factor (e.g., "donor")

# Forced contrast: logFC = ALT - REF
COND_REF = "Untreated"
COND_ALT = "CD47_Inhibitor"

PB_MIN_CELLS_PER_UNIT = 30

# How many pseudobulk replicates to create per condition (when SAMPLE_KEY is None)
REPLICATES_PER_SAMPLE = 4  # <-- change this to the number of pseudobulks you want

# DE / plotting thresholds
FDR_THR      = 0.0000000001
LFC_THR      = 1.0

# Annotation controls
TOP_ANNOTATE          = 10
ANNO_MIN_LOGCPM       = 1.0   # must exceed this mean expression
ANNO_MIN_ABS_LFC      = 1   # must exceed this absolute LFC

# Toggle annotation strategy: "smart" or "simple"
ANNOTATION_MODE = "simple"     # change to "simple" to avoid smart ranking

# Point-size scaling
POINTSIZE_MIN = 6.0
POINTSIZE_MAX = 60.0
POINTSIZE_POWER = 0.5  # sqrt scaling of detection rate


# ------------------- Helper Functions -------------------
def filter_genes_drop_malat1_mt(A):
    v  = pd.Index(A.var_names.astype(str))
    up = v.str.upper()
    keep = (~up.str.contains("MALAT1")) & (~up.str.startswith(("MT-","MT_","MT."))) & (~up.str.contains("ENSG"))
    if keep.sum() == 0:
        raise ValueError("After filtering MALAT1/MT genes, no genes remain.")
    return A[:, keep].copy()

def pick_counts_layer(A):
    for cand in ("counts","SCT_counts","raw_counts"):
        if cand in A.layers:
            return cand
    X = A.X
    if sp.issparse(X):
        ok = np.all(np.equal(np.asarray(X.data), np.asarray(X.data).astype(int)))
    else:
        ok = np.all(np.equal(X, X.astype(int)))
    if not ok:
        raise ValueError("No counts-like layer found (counts/SCT_counts/raw_counts) and .X is not integer.")
    return None

def _sum_over_cells(X):
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel()
    return X.sum(axis=0)

def _gene_detection(A, counts_layer):
    """
    Return DataFrame indexed by gene with:
      n_cells_detected: number of cells with counts > 0
      pct_cells_detected: proportion of cells with counts > 0
    """
    X = A.layers[counts_layer] if counts_layer is not None else A.X
    if sp.issparse(X):
        nnz = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        nnz = (X > 0).sum(axis=0)
    n_cells = A.n_obs
    out = pd.DataFrame({
        "n_cells_detected": nnz.astype(float),
        "pct_cells_detected": nnz.astype(float) / float(n_cells)
    }, index=A.var_names)
    out.index.name = None
    return out, n_cells

def make_pseudobulk(A, counts_layer, cond_key, sample_key=None,
                    celltype_key=None, min_cells=30, reps_per_sample=1):
    """
    Build pseudobulk count matrices.

    - If sample_key is None: split cells within each condition into
      `reps_per_sample` random chunks → pseudobulk "samples" per condition.
    - If sample_key is set: aggregation is per (sample_key, [celltype_key]),
      optionally further split into `reps_per_sample` chunks.
    """
    if cond_key not in A.obs:
        raise KeyError(f"obs['{cond_key}'] missing")
    A.obs[cond_key] = A.obs[cond_key].astype("category")

    # Keep only the two forced conditions, in case extra exist
    if set(A.obs[cond_key].cat.categories) != {COND_REF, COND_ALT}:
        A = A[A.obs[cond_key].isin([COND_REF, COND_ALT])].copy()
        A.obs[cond_key] = A.obs[cond_key].astype("category")
    if A.obs[cond_key].nunique() != 2:
        raise ValueError(f"Need exactly two conditions: {COND_REF} and {COND_ALT}")

    Xmat = A.layers[counts_layer] if counts_layer is not None else A.X

    # Grouping columns for initial units
    group_cols = []
    if sample_key is not None and sample_key in A.obs:
        group_cols.append(sample_key)
    if celltype_key is not None and celltype_key in A.obs:
        group_cols.append(celltype_key)

    gdf = A.obs[[c for c in group_cols if c] + [cond_key]].copy()
    gdf["_idx"] = np.arange(A.n_obs)

    if group_cols:
        # First aggregate at (sample, [cell_type]) level
        keep_units = gdf.groupby(group_cols)["_idx"].transform("size") >= min_cells
        gdf = gdf.loc[keep_units]
        if gdf.empty:
            raise ValueError("No pseudo-bulk units remain after min_cells filter.")
        base_id = gdf[group_cols[0]].astype(str)
        for extra in group_cols[1:]:
            base_id = base_id + "|" + gdf[extra].astype(str)
    else:
        # No sample IDs: base units are just the two conditions
        base_id = gdf[cond_key].astype(str)

    gdf = gdf.assign(_base=base_id.values)

    units, metas = [], []
    for unit, sub in gdf.groupby("_base"):
        idx = sub["_idx"].to_numpy()

        # Split cells of this base unit into pseudobulk replicates
        if reps_per_sample > 1:
            idx = np.random.permutation(idx)
            chunks = np.array_split(idx, reps_per_sample)
        else:
            chunks = [idx]

        for rep_i, chunk in enumerate(chunks, start=1):
            if len(chunk) == 0:
                continue

            vec = _sum_over_cells(Xmat[chunk])
            units.append(vec)

            # "sample" column:
            # - if sample_key exists: keep the biological sample ID
            # - otherwise: treat each pseudobulk as its own sample
            if sample_key is not None and sample_key in sub.columns:
                sample_name = str(sub.iloc[0][sample_key])
            else:
                sample_name = f"{unit}__rep{rep_i}"

            row = {
                "unit":      f"{unit}__rep{rep_i}",
                "condition": str(sub.iloc[0][cond_key]),
                "sample":    sample_name,
            }
            if celltype_key is not None and celltype_key in sub.columns:
                row["cell_type"] = str(sub.iloc[0][celltype_key])
            metas.append(row)

    counts = np.vstack(units)
    counts_df = pd.DataFrame(
        counts.T,
        index=A.var_names,
        columns=[m["unit"] for m in metas],
    )
    meta_df = pd.DataFrame(metas).set_index("unit")
    meta_df["condition"] = meta_df["condition"].astype("category")
    if "cell_type" in meta_df:
        meta_df["cell_type"] = meta_df["cell_type"].astype("category")

    return counts_df, meta_df

def have_replicates(meta_df):
    """
    Decide if we have enough replication to run edgeR/DESeq2 with p-values.

    - If SAMPLE_KEY is defined: require >=2 distinct *samples* per condition
      (true sample-level replication).
    - If SAMPLE_KEY is None: treat each pseudobulk unit as an independent
      sample and require >=2 pseudobulks per condition.
    """
    if SAMPLE_KEY is None:
        n_ref = (meta_df["condition"] == COND_REF).sum()
        n_alt = (meta_df["condition"] == COND_ALT).sum()
        return (n_ref >= 2) and (n_alt >= 2)
    else:
        return (
            meta_df.loc[meta_df["condition"] == COND_REF, "sample"].nunique() >= 2
            and meta_df.loc[meta_df["condition"] == COND_ALT, "sample"].nunique() >= 2
        )

# ---------- Smart annotation ranking ----------
def _smart_annot_rank(df, top_n, min_logcpm=ANNO_MIN_LOGCPM, min_abs_lfc=ANNO_MIN_ABS_LFC):
    """
    df must contain: gene, logCPM, logFC (and optionally mlog10FDR).
    Returns top_n by equal-weight score of expression and |FC| past thresholds.
    """
    d = df.copy()
    d["abs_logFC"] = d["logFC"].abs()

    # Must exceed both minima
    cand = d[(d["logCPM"] > min_logcpm) & (d["abs_logFC"] > min_abs_lfc)].copy()
    if cand.empty:
        return cand

    # Excess over thresholds
    cand["expr_excess"] = cand["logCPM"]   - min_logcpm
    cand["fc_excess"]   = cand["abs_logFC"] - min_abs_lfc

    # Robust scaling by 90th percentile to balance terms
    def _scale(vals):
        vals = np.asarray(vals)
        if vals.size == 0:
            return vals
        p90 = np.percentile(vals, 90)
        denom = p90 if p90 > 0 else (vals.max() if vals.max() > 0 else 1.0)
        return vals / denom

    cand["score"] = 0.5 * _scale(cand["expr_excess"].values) + 0.5 * _scale(cand["fc_excess"].values)
    cand = cand.sort_values("score", ascending=False).head(top_n)
    return cand

def _sizes_from_detection(detected_prop):
    """
    Map detection proportion [0..1] -> point size.
    Use sqrt scaling to reduce dynamic range.
    """
    p = np.clip(np.asarray(detected_prop).astype(float), 0.0, 1.0)
    p = np.power(p, POINTSIZE_POWER)
    return POINTSIZE_MIN + (POINTSIZE_MAX - POINTSIZE_MIN) * p


# ---------- Plotting ----------
def volcano_plot(df, title, outbase, fdr_thr=0.05, lfc_thr=1.0, annotate_top=200):
    """
    df must have: gene, logFC, FDR, logCPM, pct_cells_detected
    """
    d = df[["gene","logFC","FDR","logCPM","pct_cells_detected"]].dropna(subset=["gene","logFC","FDR"]).copy()
    d["mlog10FDR"] = -np.log10(d["FDR"].clip(lower=np.finfo(float).tiny))

    sig_up = (d["FDR"] < fdr_thr) & (d["logFC"] >=  lfc_thr)
    sig_dn = (d["FDR"] < fdr_thr) & (d["logFC"] <= -lfc_thr)
    ns     = ~(sig_up | sig_dn)

    # For subtitle
    n_up = int(sig_up.sum())
    n_dn = int(sig_dn.sum())

    sizes = _sizes_from_detection(d["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6,5.6), dpi=220)
    ax.scatter(d.loc[ns,"logFC"],     d.loc[ns,"mlog10FDR"],     s=sizes[ns], c="#c7c7c7", alpha=0.45, linewidths=0, label="NS")
    ax.scatter(d.loc[sig_dn,"logFC"], d.loc[sig_dn,"mlog10FDR"], s=sizes[sig_dn], c="#3b82f6", alpha=0.85, linewidths=0, label="Down")
    ax.scatter(d.loc[sig_up,"logFC"], d.loc[sig_up,"mlog10FDR"], s=sizes[sig_up], c="#ef4444", alpha=0.85, linewidths=0, label="Up")

    ax.axvline( lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axvline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-np.log10(fdr_thr), color="black", ls="--", lw=1, alpha=0.5)

    ax.set_xlabel("log2 fold-change (CD47_Inhibitor − Untreated)")
    ax.set_ylabel("-log10(FDR)")

    subtitle = f"Up: {n_up}   Down: {n_dn}   (FDR < {fdr_thr}, |log2FC| ≥ {lfc_thr})"
    ax.set_title(f"{title}\n{subtitle}", fontsize=9)

    for spine in ("top","right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # Gene labels
    texts = []
    if ANNOTATION_MODE == "smart":
        # --- original smart ranking ---
        d_sig = d.loc[sig_up | sig_dn, ["gene","logFC","logCPM","mlog10FDR"]].copy()
        top = _smart_annot_rank(d_sig, annotate_top)
    else:
        # --- simple ranking: top by FDR then |logFC| among significant genes ---
        d_sig = d.loc[sig_up | sig_dn, ["gene","logFC","FDR","logCPM","mlog10FDR"]].copy()
        if not d_sig.empty:
            d_sig["abs_logFC"] = d_sig["logFC"].abs()
            top = (
                d_sig
                .sort_values(["FDR", "abs_logFC"], ascending=[True, False])
                .head(annotate_top)
            )
        else:
            top = d_sig

    if not top.empty:
        # build a lookup for y positions
        y_lookup = d.set_index("gene")["mlog10FDR"].to_dict()
        for _, r in top.iterrows():
            x = r["logFC"]
            y = y_lookup.get(r["gene"], r.get("mlog10FDR", np.nan))
            texts.append(ax.text(x, y, r["gene"], fontsize=6, ha="center", va="bottom"))
        adjust_text(
            texts, ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(arrowstyle="-", lw=0.5, color="black", alpha=0.6),
        )

    fig.tight_layout()
    fig.savefig(outbase + ".pdf")
    fig.savefig(outbase + ".png", dpi=220)
    plt.close(fig)

def ma_plot(df, title, outbase, lfc_thr=1.0, annotate_top=200):
    # df must have: gene, logCPM, logFC, pct_cells_detected
    x = df["logCPM"].values   # log2 CPM
    y = df["logFC"].values
    sizes = _sizes_from_detection(df["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6,5.6), dpi=220)
    ax.scatter(x, y, s=sizes, alpha=0.65, linewidths=0, c="#a8a8a8")
    ax.axhline( lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel("mean expression (logCPM)")
    ax.set_ylabel("log2 fold-change (CD47_Inhibitor − Untreated)")
    ax.set_title(title + "\n(no replicates: effect sizes only)", fontsize=9)
    for spine in ("top","right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # Gene labels
    texts = []
    if ANNOTATION_MODE == "smart":
        # --- original smart ranking ---
        pick = df[["gene","logCPM","logFC"]].copy()
        top = _smart_annot_rank(pick, annotate_top)
    else:
        # --- simple ranking: top by |logFC| ---
        pick = df[["gene","logCPM","logFC"]].copy()
        if not pick.empty:
            pick["abs_logFC"] = pick["logFC"].abs()
            top = pick.sort_values("abs_logFC", ascending=False).head(annotate_top)
        else:
            top = pick

    if not top.empty:
        # place texts at their (x,y) coordinates
        lut_x = df.set_index("gene")["logCPM"].to_dict()
        lut_y = df.set_index("gene")["logFC"].to_dict()
        for _, r in top.iterrows():
            gx = lut_x.get(r["gene"], r["logCPM"])
            gy = lut_y.get(r["gene"], r["logFC"])
            texts.append(ax.text(gx, gy, r["gene"], fontsize=6, ha="center", va="bottom"))
        adjust_text(
            texts, ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(arrowstyle="-", lw=0.5, color="black", alpha=0.6),
        )

    fig.tight_layout()
    fig.savefig(outbase + ".pdf")
    fig.savefig(outbase + ".png", dpi=220)
    plt.close(fig)


# ------------------- R-side via rpy2 -------------------
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

edgeR = importr("edgeR")
limma = importr("limma")
deseq2 = importr("DESeq2")

# Full edgeR pipeline: force group levels so contrast is ALT - REF
ro.r(f"""
run_edgeR <- function(counts, group, replicate=NULL, cell_type=NULL) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  if (!is.null(replicate)) replicate <- factor(replicate)
  if (!is.null(cell_type)) cell_type <- factor(cell_type)

  library(edgeR); library(limma)

  if (!is.null(cell_type)) {{
    # Per cell-type design
    groupct <- factor(paste0(group, ".", cell_type))

    if (is.null(replicate)) {{
      design <- model.matrix(~ 0 + groupct)
    }} else {{
      design <- model.matrix(~ 0 + groupct + replicate)
    }}
    colnames(design) <- sub("^groupct", "", colnames(design))

    y <- DGEList(counts=counts)
    keep <- filterByExpr(y, design=design)
    y <- y[keep,, keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)

    res_list <- list()
    for (ct in levels(cell_type)) {{
      cname1 <- paste0("{COND_ALT}", ".", ct)
      cname0 <- paste0("{COND_REF}", ".", ct)
      contrast <- makeContrasts(
        contrasts=paste0("`", cname1, "`-`", cname0, "`"),
        levels=design
      )
      qlf <- glmQLFTest(fit, contrast=contrast)
      tt  <- topTags(qlf, n=Inf)$table
      tt$gene <- rownames(tt)
      res_list[[ct]] <- tt
    }}
    return(res_list)

  }} else {{
    # No cell-type stratification
    if (is.null(replicate)) {{
      design <- model.matrix(~ 0 + group)
    }} else {{
      design <- model.matrix(~ 0 + group + replicate)
    }}
    colnames(design) <- sub("^group", "", colnames(design))

    y <- DGEList(counts=counts)
    keep <- filterByExpr(y, design=design)
    y <- y[keep,, keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)

    contrast <- makeContrasts(
      contrasts="`{COND_ALT}`-`{COND_REF}`",
      levels=design
    )
    qlf <- glmQLFTest(fit, contrast=contrast)
    tt  <- topTags(qlf, n=Inf)$table
    tt$gene <- rownames(tt)
    return(tt)
  }}
}}
""")

# Effect-size-only path (no p-values), forced levels for ALT - REF
ro.r(f"""
edgeR_effect_only <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  if (nlevels(group) != 2) stop("Need exactly 2 groups")
  y <- edgeR::DGEList(counts=counts)
  y <- edgeR::calcNormFactors(y)
  lcpms <- edgeR::cpm(y, log=TRUE, prior.count=0.5)  # log2 CPM
  mu0 <- rowMeans(lcpms[, group=="{COND_REF}", drop=FALSE])
  mu1 <- rowMeans(lcpms[, group=="{COND_ALT}", drop=FALSE])
  out  <- data.frame(
    gene   = rownames(counts),
    logFC  = mu1 - mu0,           # CD47_Inhibitor − Untreated
    logCPM = rowMeans(lcpms),
    PValue = NA_real_,
    FDR    = NA_real_
  )
  rownames(out) <- NULL
  out
}}
""")

# DESeq2 pipeline: same forced contrast ALT - REF
ro.r(f"""
run_DESeq2 <- function(counts, group, replicate=NULL) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))

  suppressPackageStartupMessages(library(DESeq2))

  coldata <- data.frame(group=group)
  rownames(coldata) <- colnames(counts)

  if (!is.null(replicate)) {{
    coldata$replicate <- factor(replicate)
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=coldata,
                                  design=~ replicate + group)
  }} else {{
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=coldata,
                                  design=~ group)
  }}

  dds <- DESeq(dds)

  res <- results(dds, contrast=c("group", "{COND_ALT}", "{COND_REF}"))

  # baseMean is mean normalized count per gene; use log2(baseMean + 0.5) ~ logCPM-like
  logCPM <- log2(res$baseMean + 0.5)

  out <- data.frame(
    gene   = rownames(res),
    logFC  = res$log2FoldChange,
    logCPM = logCPM,
    PValue = res$pvalue,
    FDR    = res$padj
  )
  rownames(out) <- NULL
  out
}}
""")

def run_edgeR(counts_df, meta_df, out_prefix=None, title_prefix=None,
              fdr_thr=FDR_THR, lfc_thr=LFC_THR):
    # Ensure columns of counts_df are in the same order as rows of meta_df
    counts_df = counts_df.loc[:, meta_df.index]

    group = meta_df["condition"].astype(str).values
    replicate = (
        meta_df[REPLICATE_KEY].astype(str).values
        if (REPLICATE_KEY is not None and REPLICATE_KEY in meta_df.columns)
        else None
    )
    cell_type = (
        meta_df["cell_type"].astype(str).values
        if "cell_type" in meta_df.columns
        else None
    )

    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    rfun = ro.globalenv["run_edgeR"]

    if cell_type is not None:
        # Per cell-type analysis
        with localconverter(conv):
            if replicate is None:
                # Let R use default replicate=NULL
                res_list = rfun(counts_df, group, cell_type=cell_type)
            else:
                # Provide both replicate and cell_type
                res_list = rfun(counts_df, group, replicate, cell_type)
        return res_list  # handle per-CT in caller
    else:
        # No cell-type stratification
        with localconverter(conv):
            if replicate is None:
                dfR = rfun(counts_df, group)
            else:
                dfR = rfun(counts_df, group, replicate)
            df = ro.conversion.rpy2py(dfR)

        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)
        return df

def effect_only_table(counts_df, meta_df):
    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    with localconverter(conv):
        dfR = ro.globalenv["edgeR_effect_only"](counts_df, meta_df["condition"].astype(str).values)
        df  = ro.conversion.rpy2py(dfR)
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df[["gene","logFC","logCPM","PValue","FDR"]]

def run_DESeq2(counts_df, meta_df):
    # Ensure columns of counts_df are in the same order as rows of meta_df
    counts_df = counts_df.loc[:, meta_df.index]

    group = meta_df["condition"].astype(str).values
    replicate = (
        meta_df[REPLICATE_KEY].astype(str).values
        if (REPLICATE_KEY is not None and REPLICATE_KEY in meta_df.columns)
        else None
    )

    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    rfun = ro.globalenv["run_DESeq2"]

    with localconverter(conv):
        if replicate is None:
            dfR = rfun(counts_df, group)
        else:
            dfR = rfun(counts_df, group, replicate)
        df = ro.conversion.rpy2py(dfR)

    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df


# ------------------- Main Loop -------------------
def quick_pca_plot(counts_df, meta_df, title, out_png):
    try:
        from sklearn.decomposition import PCA
    except ImportError:
        return
    X = counts_df.T.values
    lib = X.sum(1, keepdims=True)
    Xnorm = np.log1p(1e6 * X / np.maximum(lib, 1))
    pc = PCA(n_components=2).fit_transform(Xnorm)
    fig, ax = plt.subplots(figsize=(4.5,4))
    ax.scatter(pc[:,0], pc[:,1], s=40, c=pd.Categorical(meta_df["condition"]).codes)
    for i, t in enumerate(meta_df.index):
        ax.text(pc[i,0], pc[i,1], t, fontsize=6, alpha=0.6)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


for label, A0 in adatas.items():
    print(f"\n### {label} ###")
    A = filter_genes_drop_malat1_mt(A0)
    if CONDITION_KEY not in A.obs:
        raise KeyError(f"{label}: obs['{CONDITION_KEY}'] missing")
    # Keep only the two target conditions in the forced order
    A = A[A.obs[CONDITION_KEY].isin([COND_REF, COND_ALT])].copy()
    A.obs[CONDITION_KEY] = A.obs[CONDITION_KEY].astype("category")

    counts_layer = pick_counts_layer(A)
    print(f"{label}: using counts from layer = {counts_layer or 'X (integers)'}")

    # Detection stats BEFORE pseudobulk (for point sizes)
    det_df, n_cells_total = _gene_detection(A, counts_layer)

    counts_df, meta_df = make_pseudobulk(
        A, counts_layer,
        cond_key=CONDITION_KEY,
        sample_key=SAMPLE_KEY,
        celltype_key=CELLTYPE_KEY,
        min_cells=PB_MIN_CELLS_PER_UNIT,
        reps_per_sample=REPLICATES_PER_SAMPLE
    )

    quick_pca_plot(counts_df, meta_df,
                   title=f"{label} pseudobulk PCA (log CPM)",
                   out_png=os.path.join(OUTDIR, f"{label}_pseudobulk_PCA.png"))

    base_edgeR  = os.path.join(OUTDIR, f"{label}_pb_edgeR")
    base_DESeq2 = os.path.join(OUTDIR, f"{label}_pb_DESeq2")

    if have_replicates(meta_df):
        # ---------- edgeR QLF ----------
        title_edgeR = f"{label} — {COND_ALT} vs {COND_REF} (edgeR QLF)"
        df_edgeR = run_edgeR(counts_df, meta_df,
                             out_prefix=base_edgeR,
                             title_prefix=title_edgeR,
                             fdr_thr=FDR_THR, lfc_thr=LFC_THR)
        # merge detection
        df_edgeR = df_edgeR.merge(det_df, left_on="gene", right_index=True, how="left")
        keep_cols = [c for c in ["gene","logFC","logCPM","F","PValue","FDR",
                                 "n_cells_detected","pct_cells_detected"] if c in df_edgeR.columns]
        df_edgeR = df_edgeR[keep_cols]
        df_edgeR.to_csv(f"{base_edgeR}.csv", index=False)
        volcano_plot(df_edgeR, title=title_edgeR,
                     outbase=f"{base_edgeR}__volcano",
                     fdr_thr=FDR_THR, lfc_thr=LFC_THR,
                     annotate_top=TOP_ANNOTATE)

        # ---------- DESeq2 ----------
        title_DE = f"{label} — {COND_ALT} vs {COND_REF} (DESeq2)"
        df_DE = run_DESeq2(counts_df, meta_df)
        df_DE = df_DE.merge(det_df, left_on="gene", right_index=True, how="left")
        keep_cols_DE = [c for c in ["gene","logFC","logCPM","PValue","FDR",
                                    "n_cells_detected","pct_cells_detected"] if c in df_DE.columns]
        df_DE = df_DE[keep_cols_DE]
        df_DE.to_csv(f"{base_DESeq2}.csv", index=False)
        volcano_plot(df_DE, title=title_DE,
                     outbase=f"{base_DESeq2}__volcano",
                     fdr_thr=FDR_THR, lfc_thr=LFC_THR,
                     annotate_top=TOP_ANNOTATE)

    else:
        print(f"{label}: no sample-level replication detected → reporting effect sizes only (no p-values/FDR).")
        df = effect_only_table(counts_df.loc[:, meta_df.index], meta_df)
        # merge detection
        df = df.merge(det_df, left_on="gene", right_index=True, how="left")
        base = base_edgeR  # keep naming consistent
        df.to_csv(base + "__effect_only.csv", index=False)
        ma_plot(df,
                title=f"{label} — {COND_ALT} vs {COND_REF}",
                outbase=base + "__effect_only_MA",
                lfc_thr=LFC_THR, annotate_top=TOP_ANNOTATE)

print("\nDone. Outputs in:", OUTDIR)











#################################################
############# DNA ANALYSIS ###############
#################################################

cell_lines = [
    "Sc_MWD_WSU_CD47_Inhibitor",
    "Sc_MWD_WSU_Untreated",
]

BL = '/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz'
GG = snap.genome.GRCh38

# Set parameters for each mark as previously used
mark_params = {
    "ac": {"n": 120000, "bs": 600},
    "me3": {"n": 80000, "bs": 5000}
}

for mark in ['ac','me3']:
    mark_name = "H3K27ac" if mark == "ac" else "H3K27me3"
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
    # Map samples -> condition
    cond_map = {
        "Sc_MWD_WSU_CD47_Inhibitor": "CD47_Inhibitor",
        "Sc_MWD_WSU_Untreated": "Untreated"
    }
    adata.obs["condition"] = adata.obs["sampleMerge"].map(cond_map).astype("category")
    adata.uns["reference_sequences"] = adata_sample.uns.get("reference_sequences", "")

    adata.write(f"Cd47i_Untreated_FULL{mark}_merged_processed.h5ad")




# ------------------- Imports -------------------
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad

# ------------------- Small utilities -------------------
def _first_present(A: ad.AnnData, candidates):
    for c in candidates:
        if c in A.obs.columns:
            return c
    return None

def _ensure_categorical_condition(A: ad.AnnData):
    if "condition" not in A.obs:
        raise ValueError("AnnData is missing obs['condition'].")
    if not pd.api.types.is_categorical_dtype(A.obs["condition"]):
        A.obs["condition"] = A.obs["condition"].astype("category")
    return list(A.obs["condition"].cat.categories)

def _cond_palette(ordered_conditions):
    pal = sns.color_palette("Set2", n_colors=len(ordered_conditions))
    return {c: pal[i] for i, c in enumerate(ordered_conditions)}

def _set_log10p1_axis_round_ticks(ax, raw_values, max_ticks=6, q=(1, 99), pad=0.08):
    """
    Axis is linear in log10(x+1) space, but tick *labels* are raw counts.
    Picks round-ish ticks (1/2/3/4/5/6/8 * 10^k) within the data range.
    Uses quantiles to avoid a single outlier forcing ugly ticks.
    """
    raw = pd.to_numeric(pd.Series(raw_values), errors="coerce").to_numpy()
    raw = raw[np.isfinite(raw)]
    if raw.size == 0:
        return
    raw = raw[raw >= 0]

    lo = float(np.percentile(raw, q[0]))
    hi = float(np.percentile(raw, q[1]))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return
    if hi <= lo:
        lo = float(np.min(raw))
        hi = float(np.max(raw))
        if hi <= lo:
            hi = lo + 1.0

    # Avoid anchoring at 0; keep lower bound above ~0 unless data truly near 0
    lo = max(lo, 1.0)

    multipliers = [1, 2, 3, 4, 5, 6, 8]
    min_k = int(np.floor(np.log10(lo))) - 1
    max_k = int(np.ceil(np.log10(hi))) + 1

    ticks_all = []
    for k in range(min_k, max_k + 1):
        base = 10 ** k
        for m in multipliers:
            t = m * base
            if t > 0:
                ticks_all.append(int(t))
    ticks_all = sorted(set(ticks_all))

    # ticks within the core range
    inrange = [t for t in ticks_all if t >= lo and t <= hi]
    if len(inrange) < 2:
        # fallback: just use rounded endpoints
        inrange = sorted(set([int(round(lo)), int(round(hi))]))

    # choose ~evenly spaced ticks in log space
    log_in = np.log10(np.array(inrange, dtype=float) + 1.0)
    desired = np.linspace(log_in.min(), log_in.max(), min(max_ticks, len(inrange)))
    chosen = []
    for d in desired:
        idx = int(np.argmin(np.abs(log_in - d)))
        chosen.append(inrange[idx])
    chosen = sorted(set(chosen))

    # expand to cover lo/hi with nearest outside ticks (still round)
    if chosen and chosen[0] > lo:
        below = [t for t in ticks_all if t < lo]
        if below:
            chosen = [below[-1]] + chosen
    if chosen and chosen[-1] < hi:
        above = [t for t in ticks_all if t > hi]
        if above:
            chosen = chosen + [above[0]]

    # downsample if expansion made too many
    if len(chosen) > max_ticks:
        logc = np.log10(np.array(chosen, dtype=float) + 1.0)
        desired = np.linspace(logc.min(), logc.max(), max_ticks)
        new = []
        for d in desired:
            idx = int(np.argmin(np.abs(logc - d)))
            new.append(chosen[idx])
        chosen = sorted(set(new))

    # apply ticks/labels in log space
    ytick_pos = np.log10(np.array(chosen, dtype=float) + 1.0)
    ax.set_yticks(ytick_pos)
    ax.set_yticklabels([f"{int(t):,}" for t in chosen])

def _annotate_group_counts_and_raw_median(ax, df, x_col, y_plot_col, raw_col, y_is_log10p1=False):
    """
    Annotate N and RAW median (not median in transformed space).
    This fixes the common bug where median(log(x)) is back-transformed.
    """
    ymin, ymax = ax.get_ylim()
    headroom = 0.12 * (ymax - ymin)
    ax.set_ylim(ymin, ymax + headroom)
    ypos = ymax + 0.04 * (ymax - ymin)

    # preserve x order as drawn
    xticks = [t.get_text() for t in ax.get_xticklabels()]
    for xt, group in enumerate(xticks):
        sub = df[df[x_col] == group]
        if sub.shape[0] == 0:
            continue
        raw_med = float(np.median(sub[raw_col].to_numpy()))
        ax.text(
            xt, ypos,
            f"{sub.shape[0]} cells\n{int(round(raw_med)):,} median",
            ha="center", va="bottom", fontsize=10
        )

# ------------------- Core plotting helpers -------------------
def violin_by_condition(
    A: ad.AnnData,
    value_col: str,
    ylabel: str,
    outpath: str,
    title: str,
    log10_plus1: bool,
    order=None,
    palette=None,
):
    if value_col not in A.obs.columns:
        return  # silently skip missing metric

    raw = pd.to_numeric(A.obs[value_col], errors="coerce").fillna(0).to_numpy()
    raw = np.clip(raw, 0, None)

    if log10_plus1:
        y = np.log10(raw + 1.0)
    else:
        y = raw

    df = pd.DataFrame({
        "condition": A.obs["condition"].astype(str).values,
        "y": y,
        "raw": raw,
    })

    if order is None:
        order = sorted(df["condition"].unique())
    if palette is None:
        palette = _cond_palette(order)

    fig, ax = plt.subplots(figsize=(6.8, 5.0))

    # cut=0 avoids violins extending beyond observed range (often looks "wrong" after transforms)
    sns.violinplot(
        data=df, x="condition", y="y",
        order=order, palette=palette,
        inner=None, cut=2, linewidth=1.2, ax=ax
    )

    # boxplot overlay is in the same plotted space (log if log10_plus1=True)
    sns.boxplot(
        data=df, x="condition", y="y",
        order=order, ax=ax, width=0.22,
        showcaps=True, dodge=False,
        boxprops={'facecolor':'none','linewidth':1.5},
        whiskerprops={'linewidth':1.5},
        medianprops={'color':'black'},
        flierprops={'marker':'','markersize':0}
    )

    ymin = np.nanmin(df["y"].to_numpy())
    ymax = np.nanmax(df["y"].to_numpy())
    pad = 0.04 * (ymax - ymin if ymax > ymin else 1.0)
    ax.set_ylim(ymin - 0.3, ymax + 0.3)

    if log10_plus1:
        _set_log10p1_axis_round_ticks(ax, df["raw"].values, max_ticks=6, q=(0, 100), pad=0.15)

    _annotate_group_counts_and_raw_median(
        ax, df, x_col="condition",
        y_plot_col="y", raw_col="raw",
        y_is_log10p1=log10_plus1
    )

    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath)
    plt.close(fig)

# ------------------- Public entrypoint -------------------
def make_condition_qc_plots(rna_ad: ad.AnnData, ac_ad: ad.AnnData, me3_ad: ad.AnnData,
                           outdir: str, prefix: str = "Showcase"):
    """
    Create condition-split violins for:
      - RNA: UMI per cell (log10(x+1)) and nGenes per cell (log10(x+1))
      - AC: "Fragments per cell" (log10(x+1))
      - ME3: "Fragments per cell" (log10(x+1))

    (No TSSE violins and no scatter plots.)
    """
    os.makedirs(outdir, exist_ok=True)

    # Prefer true fragment counters; only fall back to nCount_ATAC if no other option exists.
    FRAG_CANDIDATES_STRICT = [
        "n_fragment", "n_fragments", "fragments", "nFragments", "passed_filters"
    ]
    FRAG_CANDIDATES_FALLBACK = FRAG_CANDIDATES_STRICT + ["nCount_ATAC"]

    # ----- RNA -----
    if rna_ad is not None:
        rna_order = _ensure_categorical_condition(rna_ad)
        rna_pal = _cond_palette(rna_order)

        umi_col = _first_present(rna_ad, ["total_counts", "n_counts", "nCount_RNA", "umi", "counts", "SCT_counts"])
        genes_col = _first_present(rna_ad, ["n_genes_by_counts", "nFeature_RNA", "genes"])

        if umi_col is not None:
            violin_by_condition(
                rna_ad, umi_col, "UMIs per cell",
                os.path.join(outdir, f"{prefix}_RNA_UMI_byCondition.pdf"),
                title=f"RNA – UMIs per Cell (by condition)\nCells: {rna_ad.n_obs:,}",
                log10_plus1=True, order=rna_order, palette=rna_pal
            )

        if genes_col is not None:
            violin_by_condition(
                rna_ad, genes_col, "Genes per cell",
                os.path.join(outdir, f"{prefix}_RNA_nGenes_byCondition.pdf"),
                title=f"RNA – Genes per Cell (by condition)\nCells: {rna_ad.n_obs:,}",
                log10_plus1=True, order=rna_order, palette=rna_pal
            )

    # ----- H3K27ac -----
    if ac_ad is not None:
        ac_order = _ensure_categorical_condition(ac_ad)
        ac_pal = _cond_palette(ac_order)

        frag_col = _first_present(ac_ad, FRAG_CANDIDATES_FALLBACK)
        if frag_col is not None:
            violin_by_condition(
                ac_ad, frag_col, "Fragments per cell",
                os.path.join(outdir, f"{prefix}_H3K27ac_Fragments_byCondition.pdf"),
                title=f"H3K27ac – Fragments per Cell (by condition)\nCells: {ac_ad.n_obs:,}\nColumn: {frag_col}",
                log10_plus1=True, order=ac_order, palette=ac_pal
            )

    # ----- H3K27me3 -----
    if me3_ad is not None:
        me3_order = _ensure_categorical_condition(me3_ad)
        me3_pal = _cond_palette(me3_order)

        frag_col = _first_present(me3_ad, FRAG_CANDIDATES_FALLBACK)
        if frag_col is not None:
            violin_by_condition(
                me3_ad, frag_col, "Fragments per cell",
                os.path.join(outdir, f"{prefix}_H3K27me3_Fragments_byCondition.pdf"),
                title=f"H3K27me3 – Fragments per Cell (by condition)\nCells: {me3_ad.n_obs:,}\nColumn: {frag_col}",
                log10_plus1=True, order=me3_order, palette=me3_pal
            )


adata1 = sc.read_h5ad("./Cd47i_Untreated_rna_Overlap.h5ad")
adata2 = sc.read_h5ad("./Cd47i_Untreated_ac_Overlap.h5ad")
adata3 = sc.read_h5ad("./Cd47i_Untreated_me3_Overlap.h5ad")

make_condition_qc_plots(
    rna_ad=adata1,     # AnnData for RNA
    ac_ad=adata2,       # AnnData for H3K27ac
    me3_ad=adata3,     # AnnData for H3K27me3
    outdir="QC_Overlap3",
    prefix="Overlap3Data"
)


adataRNA = sc.read_h5ad("Cd47i_Untreated_FULLrna_merged_processedSCT_Called_Clean.h5ad")
adataMe3 = sc.read_h5ad("Cd47i_Untreated_FULLme3_merged_processed.h5ad")
adataAc = sc.read_h5ad("Cd47i_Untreated_FULLac_merged_processed.h5ad")

make_condition_qc_plots(
    rna_ad=adataRNA,     # AnnData for RNA
    ac_ad=adataAc,       # AnnData for H3K27ac
    me3_ad=adataMe3,     # AnnData for H3K27me3
    outdir="Showcase_QC_FULL",
    prefix="FULLData"
)


# Define the set of barcodes for each dataset after filtering
rna_barcodes = set(adataRNA.obs.index)
ac_barcodes = set(adataAc.obs.index)
me3_barcodes = set(adataMe3.obs.index)

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
plt.savefig(f'./figures/{exp}_VennPassingCells_Cd47i_Untreated.pdf')
plt.close()


# Subset adata2 and adata3 to only the cells that are present in adata1
adata1 = adataRNA[adataRNA.obs_names.isin(shared_barcodes)].copy()
adata2 = adataAc[adataAc.obs_names.isin(shared_barcodes)].copy()
adata3 = adataMe3[adataMe3.obs_names.isin(shared_barcodes)].copy()

adata1.write("./Cd47i_Untreated_rna_Overlap.h5ad")
adata2.write("./Cd47i_Untreated_ac_Overlap.h5ad")
adata3.write("./Cd47i_Untreated_me3_Overlap.h5ad")


adata = adata1.copy()
celltypes = {
    "WSU-DLCL2": (adata.obs["cluster_call"] == "WSU-DLCL2"),
    "Macrophage": (adata.obs["cluster_call"] == "Macrophage"),
}
subsets = {k: adata[v].copy() for k, v in celltypes.items() if v.sum() > 0}

def preprocess_and_embed(ad_sub, name, n_top_genes=2000):
    # Basic reprocessing for clean subset embedding
    sc.pp.highly_variable_genes(ad_sub, n_top_genes=n_top_genes, flavor="seurat_v3")
    sc.tl.pca(ad_sub, n_comps=10)
    sc.pl.pca_variance_ratio(adata, n_pcs=10, log=True, show=False, save=f'SUB{name}_pca_variance_ratio_SCT.pdf')
    sc.pp.neighbors(ad_sub, n_pcs=10, n_neighbors=30)
    sc.tl.leiden(ad_sub, key_added="leiden", resolution=0.5)
    sc.tl.umap(ad_sub)
    sc.pl.umap(ad_sub, color=["condition","leiden"], wspace=0.4,
               frameon=False, show=False, save=f"_{name}_RNAOverlap3UMAP_condition_leiden.pdf")
    return ad_sub

for name, ad_sub in subsets.items():
    subsets[name] = preprocess_and_embed(ad_sub, name)
    ad_sub.write(f"Overlap3RNA_{name}_Cd47i_Untreated_Reclust.h5ad")





#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from anndata import AnnData
from typing import Union, Optional, Dict, Any


def plot_linked_umaps_procrustes(
    adata_rna: Union[str, AnnData],
    adata_ac:  Union[str, AnnData],
    adata_me3: Union[str, AnnData],
    color_key: str,
    output_prefix: str,
    umap_key: str = "X_umap",
    fig_size=(18, 6),
    fig_dpi=300,
    lift_rna_frac: float = 0.10,
    lift_me3_frac: float = 0.10,
    gap_frac: float = 0.30,
    point_size: float = 18,
    line_width: float = 0.22,
    line_alpha: float = 0.20,
    point_alpha: float = 0.98,
    legend_max_cols: int = 1,
    rasterize: bool = True,
    palette_override: Optional[Dict[Any, str]] = None,
    max_lines: Optional[int] = None,   # set e.g. 30000 to cap; None = all lines
    line_seed: int = 0,
) -> Dict[str, Any]:
    """
    Linked UMAPs with TRANSLATION ONLY (no Procrustes / no rotation / no scaling):
      - Intersect cells across modalities
      - Take native UMAP coords from each modality
      - Horizontally translate AC and ME3 panels (gap_frac) into ONE shared axis
      - Optionally vertically lift RNA + ME3 (lift_*_frac) to reduce line overlap
      - Draw lines in the SAME coordinate system as the points (so they always land correctly)
      - Color by obs[color_key] (prefer RNA; fallback AC/ME3)
    """

    # ---------- helpers ----------
    def _ensure_adata(x):
        if isinstance(x, AnnData):
            return x
        if isinstance(x, str):
            return sc.read_h5ad(x)
        raise TypeError("adata inputs must be AnnData or file paths to .h5ad")

    def _require_umap(a: AnnData, name: str):
        if umap_key not in a.obsm_keys():
            raise KeyError(f"{name} is missing .obsm['{umap_key}'].")

    def _get_obs_series(a1: AnnData, a2: AnnData, a3: AnnData, key: str, idx) -> pd.Series:
        if key in a1.obs.columns:
            return a1.obs.loc[idx, key]
        if key in a2.obs.columns:
            return a2.obs.loc[idx, key]
        if key in a3.obs.columns:
            return a3.obs.loc[idx, key]
        raise KeyError(f"'{key}' not found in .obs of RNA/AC/ME3.")

    def x_extent(arr2):  # arr2: (n,2)
        return float(arr2[:, 0].max() - arr2[:, 0].min())

    def y_extent(arr2):
        return float(arr2[:, 1].max() - arr2[:, 1].min())

    def _make_discrete_palette(categories) -> Dict[Any, str]:
        cats = list(categories)

        # Start with explicit overrides (exact key match)
        pal = {c: None for c in cats}
        if palette_override:
            for c in cats:
                if c in palette_override:
                    pal[c] = palette_override[c]

        # Defaults (conditions; cell types; condition×celltype combos)
        # These are chosen to mirror the “related shades within group” idea from Fig 3.
        defaults = {
            # Condition-only
            "Untreated":        "#4C72B0",
            "CD47_Inhibitor":   "#DD8452",

            # Cell type-only (robust to singular/plural)
            "Macrophage":       "#8C564B",  # Fig3 Macrophages (brown)
            "Macrophages":      "#8C564B",
            "WSU-DLCL2":        "#2CA02C",  # Fig3 greens

            # Condition × cell type (4-color paired scheme)
            "Untreated__Macrophage":        "#8C564B",
            "CD47_Inhibitor__Macrophage":   "#C49C94",
            "Untreated__Macrophages":       "#8C564B",
            "CD47_Inhibitor__Macrophages":  "#C49C94",
            "Untreated__WSU-DLCL2":         "#2CA02C",
            "CD47_Inhibitor__WSU-DLCL2":    "#98DF8A",
        }

        for c in cats:
            if pal[c] is None and c in defaults:
                pal[c] = defaults[c]

        # Fill remaining categories deterministically
        remaining = [c for c in cats if pal[c] is None]
        if remaining:
            try:
                base_list = list(sc.pl.palettes.default_102)
            except Exception:
                base_list = None

            if base_list and len(remaining) <= len(base_list):
                for c, col in zip(remaining, base_list):
                    pal[c] = col
            else:
                cmap = plt.get_cmap("turbo")
                for i, c in enumerate(remaining):
                    pal[c] = mcolors.to_hex(cmap(i / max(1, (len(remaining) - 1))))

        return pal

    # ---------- I/O ----------
    a1 = _ensure_adata(adata_rna)
    a2 = _ensure_adata(adata_ac)
    a3 = _ensure_adata(adata_me3)

    _require_umap(a1, "RNA")
    _require_umap(a2, "H3K27ac")
    _require_umap(a3, "H3K27me3")

    # ---------- intersect cells ----------
    common = a1.obs_names.intersection(a2.obs_names).intersection(a3.obs_names)
    if len(common) == 0:
        raise ValueError("No overlapping cell barcodes across RNA/AC/ME3. Check obs_names consistency.")
    common = common.sort_values()

    a1 = a1[common].copy()
    a2 = a2[common].copy()
    a3 = a3[common].copy()

    umap1 = pd.DataFrame(a1.obsm[umap_key], index=common, columns=["x1", "y1"])
    umap2 = pd.DataFrame(a2.obsm[umap_key], index=common, columns=["x2", "y2"])
    umap3 = pd.DataFrame(a3.obsm[umap_key], index=common, columns=["x3", "y3"])

    # ---------- TRANSLATE PANELS ONLY (no rotation / no scaling) ----------
    w = max(
        x_extent(umap1[["x1", "y1"]].values),
        x_extent(umap2[["x2", "y2"]].values),
        x_extent(umap3[["x3", "y3"]].values),
    )
    gap = gap_frac * w

    x2min = umap2["x2"].min()
    umap2["x2"] += (umap1["x1"].max() - x2min) + gap

    x3min = umap3["x3"].min()
    umap3["x3"] += (umap2["x2"].max() - x3min) + gap

    H = max(
        y_extent(umap1[["x1", "y1"]].values),
        y_extent(umap2[["x2", "y2"]].values),
        y_extent(umap3[["x3", "y3"]].values),
    )
    umap1["y1"] += lift_rna_frac * H
    umap3["y3"] += lift_me3_frac * H

    # ---------- combine + labels ----------
    df = pd.concat([umap1, umap2, umap3], axis=1)

    lab = _get_obs_series(a1, a2, a3, color_key, common)
    if not isinstance(lab.dtype, pd.CategoricalDtype):
        lab = lab.astype("category")
    df[color_key] = lab

    cats = list(df[color_key].cat.categories)
    palette = _make_discrete_palette(cats)

    rgba_points = np.vstack([
        mcolors.to_rgba(palette.get(v, "#BDBDBD"), alpha=point_alpha)
        for v in df[color_key].astype(object).values
    ])
    rgba_lines_all = np.vstack([
        mcolors.to_rgba(palette.get(v, "#BDBDBD"), alpha=line_alpha)
        for v in df[color_key].astype(object).values
    ])

    # Optional line subsampling (for huge N)
    n = df.shape[0]
    if max_lines is None or max_lines >= n:
        idx = np.arange(n)
    else:
        rng = np.random.default_rng(line_seed)
        idx = rng.choice(n, size=max_lines, replace=False)

    # ---------- plot (single axis so lines land exactly) ----------
    fig, ax = plt.subplots(figsize=fig_size, dpi=fig_dpi)
    ax.set_title(f"RNA → H3K27ac → H3K27me3  (colored by {color_key})", fontsize=18, weight="bold")

    seg12 = np.stack([df.loc[df.index[idx], ["x1", "y1"]].values,
                      df.loc[df.index[idx], ["x2", "y2"]].values], axis=1)
    seg23 = np.stack([df.loc[df.index[idx], ["x2", "y2"]].values,
                      df.loc[df.index[idx], ["x3", "y3"]].values], axis=1)

    rgba_lines = rgba_lines_all[idx]

    lc12 = LineCollection(seg12, colors=rgba_lines, linewidths=line_width, zorder=0)
    lc23 = LineCollection(seg23, colors=rgba_lines, linewidths=line_width, zorder=0)
    lc12.set_rasterized(rasterize)
    lc23.set_rasterized(rasterize)
    ax.add_collection(lc12)
    ax.add_collection(lc23)

    for x, y, panel_name in [("x1", "y1", "RNA"), ("x2", "y2", "H3K27ac"), ("x3", "y3", "H3K27me3")]:
        sca = ax.scatter(
            df[x].values, df[y].values,
            s=point_size,
            c=rgba_points,
            edgecolors="none",
            zorder=1,
        )
        sca.set_rasterized(rasterize)

        ax.text(
            df[x].min(), df[y].max(),
            panel_name,
            fontsize=14, weight="bold",
            ha="left", va="bottom"
        )

    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.03)
    ax.axis("off")

    handles = [
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=palette[c], markeredgecolor="none",
               markersize=7, label=str(c))
        for c in cats
    ]
    ax.legend(
        handles=handles,
        title=color_key,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
        ncol=max(1, legend_max_cols),
        fontsize=9,
        title_fontsize=10,
    )

    plt.tight_layout()

    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    pdf_path = f"{output_prefix}.pdf"
    png_path = f"{output_prefix}.png"
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.savefig(png_path, bbox_inches="tight")
    plt.close(fig)

    return {
        "pdf": pdf_path,
        "png": png_path,
        "coords": df,
        "palette": palette,
        "n_cells": len(common),
        "n_lines": len(idx),
    }


if __name__ == "__main__":
    # ---- inputs ----
    RNA_PATH = "Cd47i_Untreated_FULLrna_merged_processedSCT_Called_Clean.h5ad"
    ME3_PATH = "FULLH3K27me3_Cd47i_Untreated_Clean.h5ad"
    AC_PATH  = "FULLH3K27ac_Cd47i_Untreated_Clean.h5ad"

    CONDITION_KEY = "condition"
    # If your true cell type column is different, change it here (e.g. "cell_type", "leiden_merged_type", etc.)
    CELLTYPE_KEY  = "cluster_call"

    COMBINED_KEY   = "condition__celltype"

    adataRNA = sc.read_h5ad(RNA_PATH)
    adataMe3 = sc.read_h5ad(ME3_PATH)
    adataAc  = sc.read_h5ad(AC_PATH)

    # ---- ensure combined key exists for the “both at the same time” plot ----
    # We build it in ALL three objects for robustness, but the plotting function will prefer RNA.
    for a in (adataRNA, adataAc, adataMe3):
        if CONDITION_KEY not in a.obs.columns:
            raise KeyError(f"Missing '{CONDITION_KEY}' in .obs")
        if CELLTYPE_KEY not in a.obs.columns:
            raise KeyError(f"Missing '{CELLTYPE_KEY}' in .obs (set CELLTYPE_KEY correctly)")

        a.obs[COMBINED_KEY] = (
            a.obs[CONDITION_KEY].astype(str) + "__" + a.obs[CELLTYPE_KEY].astype(str)
        ).astype("category")

    # ---- 1) linked UMAP colored by condition ----
    res_condition = plot_linked_umaps_procrustes(
        adataRNA, adataAc, adataMe3,
        color_key=CONDITION_KEY,
        output_prefix="./figures/LinkedUMAPs_by_condition",
        lift_rna_frac=0.10,
        lift_me3_frac=0.10,
        gap_frac=0.30,
        rasterize=False,
        max_lines=None,
    )
    print("Saved:", res_condition["pdf"], res_condition["png"],
          "cells:", res_condition["n_cells"], "lines:", res_condition["n_lines"])

    # ---- 2) linked UMAP colored by cell type ----
    res_celltype = plot_linked_umaps_procrustes(
        adataRNA, adataAc, adataMe3,
        color_key=CELLTYPE_KEY,
        output_prefix="./figures/LinkedUMAPs_by_celltype",
        lift_rna_frac=0.10,
        lift_me3_frac=0.10,
        gap_frac=0.30,
        rasterize=False,
        max_lines=None,
    )
    print("Saved:", res_celltype["pdf"], res_celltype["png"],
          "cells:", res_celltype["n_cells"], "lines:", res_celltype["n_lines"])

    # ---- 3) linked UMAP colored by condition × cell type (4-color paired scheme for Macrophage/WSU) ----
    # Only the Macrophage/WSU combos are explicitly colored; any other combos (if present) get auto colors.
    res_both = plot_linked_umaps_procrustes(
        adataRNA, adataAc, adataMe3,
        color_key=COMBINED_KEY,
        output_prefix="./figures/LinkedUMAPs_by_condition_and_celltype",
        lift_rna_frac=0.10,
        lift_me3_frac=0.10,
        gap_frac=0.30,
        rasterize=False,
        max_lines=None,
    )
    print("Saved:", res_both["pdf"], res_both["png"],
          "cells:", res_both["n_cells"], "lines:", res_both["n_lines"])

