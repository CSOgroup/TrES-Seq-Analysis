#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DESeq2-only pseudobulk differential analysis on H3K27ac gene activity.

Pipeline:
1) Load ATAC/AC AnnData (.h5ad)
2) Subset to TARGET_CELLTYPE + CONDITIONS
3) Compute gene activity matrix with SnapATAC2 (per-cell gene counts)
   - raw gene activity counts saved in ga.layers["counts"]
4) Build pseudobulk samples:
   - If SAMPLE_KEY provided: aggregate by SAMPLE_KEY within each condition
   - Else: split cells within each condition into REPLICATES_PER_CONDITION chunks
5) Run DESeq2 (Wald), contrast = CD47_Inhibitor vs Untreated
6) Save CSV + Volcano plot + Pseudobulk PCA

NOTE on inference:
- If SAMPLE_KEY is None and you split cells into pseudo-replicates, p-values/FDR can be wildly overconfident.
  Use true replicate IDs (donor/library) in SAMPLE_KEY whenever possible.
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ------------------- User config -------------------
AC_PATH = "FULLH3K27ac_Cd47i_Untreated_Clean.h5ad"

celltype_key  = "cluster_call"
condition_key = "condition"

TARGET_CELLTYPE = "Macrophage"
COND_REF = "Untreated"
COND_ALT = "CD47_Inhibitor"
CONDITIONS = [COND_REF, COND_ALT]

OUTDIR = "figures_deseq2_gene_activity"
os.makedirs(OUTDIR, exist_ok=True)

# Replication config (recommended: set SAMPLE_KEY to a true replicate/library/donor column)
SAMPLE_KEY = None         # e.g. "library_id" or "donor" in ac.obs
REPLICATE_KEY = None      # optional blocking factor in DESeq2 design (e.g. donor) if SAMPLE_KEY is library
REPLICATES_PER_CONDITION = 3  # used ONLY if SAMPLE_KEY is None (pseudo-reps)

PB_MIN_CELLS_PER_UNIT = 30

# Gene activity settings
AC_GENE_ACTIVITY_CACHE = True
AC_GENE_ANNO = "GRCh38"      # snap.genome.GRCh38
AC_UPSTREAM = 2500
AC_DOWNSTREAM = 2500
AC_INCLUDE_GENE_BODY = True

# Plotting thresholds
FDR_THR = 0.01
LFC_THR = 1.0
TOP_ANNOTATE = 80

# Point size from detection in single-cell gene activity
POINTSIZE_MIN = 6.0
POINTSIZE_MAX = 60.0
POINTSIZE_POWER = 0.5  # sqrt scaling

# ------------------- Helpers -------------------
def _require_obs_keys(adata, keys, label):
    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        raise ValueError(f"[{label}] Missing obs keys: {missing}. Available: {list(adata.obs.columns)}")

def _to_float32(X):
    if sp.issparse(X):
        return X.astype(np.float32)
    return np.asarray(X, dtype=np.float32)

def _sum_over_cells(X):
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel()
    return np.asarray(X).sum(axis=0)

def _gene_detection_from_layer(A, layer="counts"):
    X = A.layers[layer] if layer in A.layers else A.X
    if sp.issparse(X):
        nnz = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        nnz = (np.asarray(X) > 0).sum(axis=0)
    n_cells = A.n_obs
    out = pd.DataFrame({
        "n_cells_detected": nnz.astype(float),
        "pct_cells_detected": nnz.astype(float) / float(n_cells)
    }, index=A.var_names)
    out.index.name = None
    return out

def _sizes_from_detection(detected_prop):
    p = np.clip(np.asarray(detected_prop).astype(float), 0.0, 1.0)
    p = np.power(p, POINTSIZE_POWER)
    return POINTSIZE_MIN + (POINTSIZE_MAX - POINTSIZE_MIN) * p

def subset_target(ac):
    _require_obs_keys(ac, [celltype_key, condition_key], "AC")
    # keep only target celltype
    mask_ct = ac.obs[celltype_key].astype(str) == str(TARGET_CELLTYPE)
    ac = ac[mask_ct].copy()
    # keep only the two conditions
    mask_cond = ac.obs[condition_key].astype(str).isin([str(c) for c in CONDITIONS])
    ac = ac[mask_cond].copy()
    ac.obs[condition_key] = ac.obs[condition_key].astype("category")
    return ac

def compute_gene_activity(ac_sub, cache_h5ad):
    """
    Returns AnnData with:
      - X = gene activity counts (raw)
      - layers["counts"] = raw counts (same as X)
    """
    if AC_GENE_ACTIVITY_CACHE and os.path.exists(cache_h5ad):
        ga = sc.read_h5ad(cache_h5ad)
        if "counts" not in ga.layers:
            ga.layers["counts"] = ga.X.copy()
        print(f"[GA] Loaded cached gene activity: {cache_h5ad}")
        return ga

    import snapatac2 as snap

    # Pick gene annotation
    if AC_GENE_ANNO == "GRCh38":
        gene_anno = snap.genome.GRCh38
    else:
        gene_anno = AC_GENE_ANNO

    use_x = ("insertion" not in ac_sub.obsm)

    print(f"[GA] Computing gene activity (use_x={use_x}, up={AC_UPSTREAM}, down={AC_DOWNSTREAM}, "
          f"gene_body={AC_INCLUDE_GENE_BODY})")

    ga = snap.pp.make_gene_matrix(
        ac_sub,
        gene_anno=gene_anno,
        inplace=False,
        file=None,
        use_x=use_x,
        upstream=AC_UPSTREAM,
        downstream=AC_DOWNSTREAM,
        include_gene_body=AC_INCLUDE_GENE_BODY,
    )

    # Ensure we keep RAW counts for DESeq2
    ga.layers["counts"] = ga.X.copy()

    # Optional: keep X as float32 for plotting convenience; DE uses layers["counts"]
    ga.X = _to_float32(ga.X)

    if AC_GENE_ACTIVITY_CACHE:
        ga.write(cache_h5ad)
        print(f"[GA] Saved gene activity cache: {cache_h5ad}")

    return ga

def make_pseudobulk(ga, layer="counts"):
    """
    Returns:
      counts_df: genes x pseudobulk_samples integer
      meta_df  : pseudobulk_samples x metadata (condition, sample[, replicate])
    """
    X = ga.layers[layer] if layer in ga.layers else ga.X

    cond = ga.obs[condition_key].astype(str).values

    # base units: either true samples (SAMPLE_KEY) or per-condition (pseudo reps)
    if SAMPLE_KEY is not None:
        if SAMPLE_KEY not in ga.obs.columns:
            raise KeyError(f"SAMPLE_KEY='{SAMPLE_KEY}' not found in ga.obs")
        sample_ids = ga.obs[SAMPLE_KEY].astype(str).values
        base = pd.Series(sample_ids, index=ga.obs_names)
        # filter units with enough cells
        sizes = base.value_counts()
        keep_samples = sizes[sizes >= PB_MIN_CELLS_PER_UNIT].index
        keep_mask = base.isin(keep_samples).values
        ga2 = ga[keep_mask].copy()
        X2 = ga2.layers[layer] if layer in ga2.layers else ga2.X
        cond2 = ga2.obs[condition_key].astype(str).values
        sample2 = ga2.obs[SAMPLE_KEY].astype(str).values

        units = []
        metas = []
        for sid in pd.unique(sample2):
            idx = np.where(sample2 == sid)[0]
            if idx.size == 0:
                continue
            vec = _sum_over_cells(X2[idx, :])
            units.append(vec)
            row = {
                "unit": sid,
                "condition": pd.unique(cond2[idx])[0],
                "sample": sid,
            }
            if REPLICATE_KEY is not None and REPLICATE_KEY in ga2.obs.columns:
                row["replicate"] = str(pd.unique(ga2.obs[REPLICATE_KEY].astype(str).values[idx])[0])
            metas.append(row)

        counts = np.vstack(units)  # samples x genes
        counts = np.round(counts).astype(np.int64, copy=False)

        counts_df = pd.DataFrame(
            counts.T,
            index=ga2.var_names,
            columns=[m["unit"] for m in metas],
        )
        meta_df = pd.DataFrame(metas).set_index("unit")
        meta_df["condition"] = meta_df["condition"].astype("category")
        return counts_df, meta_df

    # --- pseudo-replicates by splitting cells within each condition ---
    units = []
    metas = []
    for c in CONDITIONS:
        idx_all = np.where(cond == str(c))[0]
        if idx_all.size < PB_MIN_CELLS_PER_UNIT:
            print(f"[PB] Condition {c}: only {idx_all.size} cells (<{PB_MIN_CELLS_PER_UNIT}); skipping.")
            continue

        idx_all = np.random.permutation(idx_all)
        chunks = np.array_split(idx_all, REPLICATES_PER_CONDITION)

        for rep_i, idx in enumerate(chunks, start=1):
            if idx.size < PB_MIN_CELLS_PER_UNIT:
                continue
            vec = _sum_over_cells(X[idx, :])
            units.append(vec)
            metas.append({
                "unit": f"{c}__rep{rep_i}",
                "condition": c,
                "sample": f"{c}__rep{rep_i}",
            })

    if not units:
        raise ValueError("No pseudobulk units created. Check cell counts / thresholds.")

    counts = np.vstack(units)  # samples x genes
    counts = np.round(counts).astype(np.int64, copy=False)

    counts_df = pd.DataFrame(
        counts.T,
        index=ga.var_names,
        columns=[m["unit"] for m in metas],
    )
    meta_df = pd.DataFrame(metas).set_index("unit")
    meta_df["condition"] = meta_df["condition"].astype("category")
    return counts_df, meta_df

def quick_pca_plot(counts_df, meta_df, title, out_png):
    try:
        from sklearn.decomposition import PCA
    except ImportError:
        return
    X = counts_df.T.values.astype(float)
    lib = X.sum(1, keepdims=True)
    Xnorm = np.log1p(1e6 * X / np.maximum(lib, 1.0))
    pc = PCA(n_components=2).fit_transform(Xnorm)
    fig, ax = plt.subplots(figsize=(4.8, 4.2), dpi=180)
    ax.scatter(pc[:, 0], pc[:, 1], s=55, c=pd.Categorical(meta_df["condition"]).codes)
    for i, name in enumerate(meta_df.index):
        ax.text(pc[i, 0], pc[i, 1], name, fontsize=6, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)

def volcano_plot(df, title, outbase, fdr_thr=0.05, lfc_thr=1.0, annotate_top=80):
    d = df[["gene","logFC","FDR","logCPM","pct_cells_detected"]].dropna(subset=["gene","logFC","FDR"]).copy()
    d["mlog10FDR"] = -np.log10(d["FDR"].astype(float).clip(lower=np.finfo(float).tiny))

    sig_up = (d["FDR"] < fdr_thr) & (d["logFC"] >=  lfc_thr)
    sig_dn = (d["FDR"] < fdr_thr) & (d["logFC"] <= -lfc_thr)
    ns     = ~(sig_up | sig_dn)

    sizes = _sizes_from_detection(d["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=220)
    ax.scatter(d.loc[ns,"logFC"],     d.loc[ns,"mlog10FDR"],     s=sizes[ns], c="#c7c7c7", alpha=0.45, linewidths=0)
    ax.scatter(d.loc[sig_dn,"logFC"], d.loc[sig_dn,"mlog10FDR"], s=sizes[sig_dn], c="#3b82f6", alpha=0.85, linewidths=0)
    ax.scatter(d.loc[sig_up,"logFC"], d.loc[sig_up,"mlog10FDR"], s=sizes[sig_up], c="#ef4444", alpha=0.85, linewidths=0)

    ax.axvline( lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axvline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-np.log10(fdr_thr), color="black", ls="--", lw=1, alpha=0.5)

    ax.set_xlabel(f"log2 fold-change ({COND_ALT} − {COND_REF})")
    ax.set_ylabel("-log10(FDR)")

    n_up = int(sig_up.sum())
    n_dn = int(sig_dn.sum())
    subtitle = f"Up: {n_up}   Down: {n_dn}   (FDR < {fdr_thr}, |log2FC| ≥ {lfc_thr})"
    ax.set_title(f"{title}\n{subtitle}", fontsize=9)

    for spine in ("top","right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # annotate top significant by FDR then |logFC|
    texts = []
    d_sig = d.loc[sig_up | sig_dn, ["gene","logFC","FDR","mlog10FDR"]].copy()
    if not d_sig.empty:
        d_sig["abs_logFC"] = d_sig["logFC"].abs()
        top = d_sig.sort_values(["FDR","abs_logFC"], ascending=[True, False]).head(annotate_top)
        y_lookup = d.set_index("gene")["mlog10FDR"].to_dict()
        for _, r in top.iterrows():
            x = r["logFC"]
            y = y_lookup.get(r["gene"], r["mlog10FDR"])
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

# ------------------- DESeq2 via rpy2 -------------------
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

deseq2 = importr("DESeq2")

ro.r(f"""
run_DESeq2 <- function(counts, group, replicate=NULL) {{
  suppressPackageStartupMessages(library(DESeq2))

  counts <- as.matrix(counts)
  # enforce integer counts (DESeq2 expects integers)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  group <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))

  coldata <- data.frame(group=group)
  rownames(coldata) <- colnames(counts)

  if (!is.null(replicate)) {{
    coldata$replicate <- factor(replicate)
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~ replicate + group)
  }} else {{
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~ group)
  }}

  dds <- DESeq(dds)

  res <- results(dds, contrast=c("group", "{COND_ALT}", "{COND_REF}"))

  # logCPM-ish for plotting
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

def run_DESeq2(counts_df, meta_df):
    counts_df = counts_df.loc[:, meta_df.index]

    group = meta_df["condition"].astype(str).values
    replicate = None
    if REPLICATE_KEY is not None and "replicate" in meta_df.columns:
        replicate = meta_df["replicate"].astype(str).values

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

# ------------------- Main -------------------
def main():
    print("[Load] AC:", AC_PATH)
    ac = sc.read_h5ad(AC_PATH)

    print("[Subset] Target celltype + conditions")
    ac = subset_target(ac)
    print(f"[AC] {TARGET_CELLTYPE} cells retained: {ac.n_obs}")
    print("[AC] Condition counts:")
    print(ac.obs[condition_key].astype(str).value_counts())

    cache_h5ad = os.path.join(
        OUTDIR, f"{TARGET_CELLTYPE}.H3K27ac_gene_activity.{COND_ALT}_vs_{COND_REF}.h5ad"
    )
    ga = compute_gene_activity(ac, cache_h5ad)

    # Detection stats at single-cell gene activity level (raw counts)
    det_df = _gene_detection_from_layer(ga, layer="counts")

    print("[Pseudobulk] Building pseudobulk count matrix")
    counts_df, meta_df = make_pseudobulk(ga, layer="counts")
    print("[Pseudobulk] Samples:", meta_df.shape[0], "Genes:", counts_df.shape[0])
    print(meta_df["condition"].value_counts())

    quick_pca_plot(
        counts_df, meta_df,
        title=f"{TARGET_CELLTYPE} H3K27ac gene activity pseudobulk PCA",
        out_png=os.path.join(OUTDIR, f"{TARGET_CELLTYPE}_GA_pseudobulk_PCA.png")
    )

    print("[DESeq2] Running DESeq2")
    df = run_DESeq2(counts_df, meta_df)

    # merge detection stats
    df = df.merge(det_df, left_on="gene", right_index=True, how="left")

    out_csv = os.path.join(
        OUTDIR, f"{TARGET_CELLTYPE}.H3K27ac_GA.DESeq2.{COND_ALT}_vs_{COND_REF}.csv"
    )
    df.to_csv(out_csv, index=False)
    print("[Save] CSV:", out_csv)

    volcano_plot(
        df.rename(columns={"pct_cells_detected": "pct_cells_detected"}),
        title=f"{TARGET_CELLTYPE} — H3K27ac gene activity (DESeq2)\n{COND_ALT} vs {COND_REF}",
        outbase=os.path.join(OUTDIR, f"{TARGET_CELLTYPE}.H3K27ac_GA.DESeq2_volcano"),
        fdr_thr=FDR_THR,
        lfc_thr=LFC_THR,
        annotate_top=TOP_ANNOTATE
    )
    print("[Done] Outputs in:", OUTDIR)

if __name__ == "__main__":
    main()
