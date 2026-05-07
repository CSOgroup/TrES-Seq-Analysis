
#!/usr/bin/env python3
# dge_markers_leiden_pseudobulk.py
# Pseudobulk cluster-marker analysis:
#   For each cell type in obs["leiden_merged_type"]:
#     cluster vs rest (logFC = cluster − rest)
# using edgeR QLF + DESeq2 with pseudobulk replicates.

import os
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ------------------- Configuration -------------------
OUTDIR = "figures_cluster_markers_pseudobulk"
os.makedirs(OUTDIR, exist_ok=True)

# <<< CHANGE THIS TO YOUR H5AD FILE >>>
H5AD_PATH = "./MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"

# Cell-type / cluster annotation
CELLTYPE_KEY = "leiden_merged_type"

# Minimum cells in a cluster to consider it for marker calling
PB_MIN_CELLS_PER_GROUP = 50

# Minimum cells per condition (cluster / rest) for pseudobulk building
PB_MIN_CELLS_PER_CONDITION = 25

# How many pseudobulk replicates to create per condition
REPLICATES_PER_CONDITION = 3

# DE / plotting thresholds
FDR_THR = 0.01
LFC_THR = 1.0

# Annotation controls
TOP_ANNOTATE     = 80
ANNO_MIN_LOGCPM  = 1.0
ANNO_MIN_ABS_LFC = 1.0

# Toggle annotation strategy: "smart" or "simple"
ANNOTATION_MODE = "simple"

# Point-size scaling
POINTSIZE_MIN   = 6.0
POINTSIZE_MAX   = 60.0
POINTSIZE_POWER = 0.5  # sqrt-like scaling of detection rate

# Binary contrast labels: logFC = ALT - REF = cluster - rest
COND_REF = "rest"
COND_ALT = "cluster"


# ------------------- Helper Functions -------------------
def filter_genes_drop_malat1_mt(A):
    v  = pd.Index(A.var_names.astype(str))
    up = v.str.upper()
    keep = (
        (~up.str.contains("MALAT1"))
        & (~up.str.startswith(("MT-", "MT_", "MT.")))
        & (~up.str.contains("ENSG"))
    )
    if keep.sum() == 0:
        raise ValueError("After filtering MALAT1/MT/ENSG genes, no genes remain.")
    return A[:, keep].copy()


def pick_counts_layer(A):
    for cand in ("counts", "SCT_counts", "raw_counts"):
        if cand in A.layers:
            return cand
    X = A.X
    if sp.issparse(X):
        ok = np.all(np.equal(np.asarray(X.data), np.asarray(X.data).astype(int)))
    else:
        ok = np.all(np.equal(X, X.astype(int)))
    if not ok:
        raise ValueError(
            "No counts-like layer found (counts/SCT_counts/raw_counts) and .X is not integer."
        )
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
    out = pd.DataFrame(
        {
            "n_cells_detected": nnz.astype(float),
            "pct_cells_detected": nnz.astype(float) / float(n_cells),
        },
        index=A.var_names,
    )
    out.index.name = None
    return out, n_cells


def make_pseudobulk_binary(A, counts_layer, cond_key,
                           min_cells=30, reps_per_cond=4):
    """
    Build pseudobulk count matrices for a *binary* contrast:

    - cond_key must have exactly two levels: COND_REF and COND_ALT.
    - Cells from each condition are shuffled and split into `reps_per_cond`
      approximately equal-sized chunks.
    - Each chunk becomes one pseudobulk "sample".
    """
    if cond_key not in A.obs:
        raise KeyError(f"obs['{cond_key}'] missing")

    A = A.copy()
    A.obs[cond_key] = A.obs[cond_key].astype("category")

    if set(A.obs[cond_key].cat.categories) != {COND_REF, COND_ALT}:
        A = A[A.obs[cond_key].isin([COND_REF, COND_ALT])].copy()
        A.obs[cond_key] = A.obs[cond_key].astype("category")

    if A.obs[cond_key].nunique() != 2:
        raise ValueError(f"Need exactly two conditions: {COND_REF} and {COND_ALT}")

    Xmat = A.layers[counts_layer] if counts_layer is not None else A.X

    units = []
    metas = []

    for cond in (COND_REF, COND_ALT):
        mask = (A.obs[cond_key] == cond).to_numpy()
        idx = np.where(mask)[0]
        if idx.size < min_cells:
            print(
                f"  [WARN] condition '{cond}' has only {idx.size} cells "
                f"(< {min_cells}); skipping this contrast."
            )
            return None, None

        idx = np.random.permutation(idx)
        chunks = np.array_split(idx, reps_per_cond)

        for rep_i, chunk in enumerate(chunks, start=1):
            if chunk.size == 0:
                continue
            vec = _sum_over_cells(Xmat[chunk])
            units.append(vec)
            unit_name = f"{cond}__rep{rep_i}"
            metas.append(
                {
                    "unit": unit_name,
                    "condition": cond,
                    "sample": unit_name,  # each pseudobulk is its own "sample"
                }
            )

    counts = np.vstack(units)
    counts_df = pd.DataFrame(
        counts.T,
        index=A.var_names,
        columns=[m["unit"] for m in metas],
    )
    meta_df = pd.DataFrame(metas).set_index("unit")
    meta_df["condition"] = meta_df["condition"].astype("category")
    return counts_df, meta_df


def have_replicates(meta_df, min_units_per_condition=2):
    """
    Decide if we have enough pseudobulk samples to run edgeR/DESeq2 with p-values:

    Require >= min_units_per_condition pseudobulks per condition.
    """
    n_ref = (meta_df["condition"] == COND_REF).sum()
    n_alt = (meta_df["condition"] == COND_ALT).sum()
    return (n_ref >= min_units_per_condition) and (n_alt >= min_units_per_condition)


# ---------- Smart annotation ranking ----------
def _smart_annot_rank(df, top_n,
                      min_logcpm=ANNO_MIN_LOGCPM,
                      min_abs_lfc=ANNO_MIN_ABS_LFC):
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
    cand["expr_excess"] = cand["logCPM"] - min_logcpm
    cand["fc_excess"] = cand["abs_logFC"] - min_abs_lfc

    # Robust scaling by 90th percentile to balance terms
    def _scale(vals):
        vals = np.asarray(vals)
        if vals.size == 0:
            return vals
        p90 = np.percentile(vals, 90)
        denom = p90 if p90 > 0 else (vals.max() if vals.max() > 0 else 1.0)
        return vals / denom

    cand["score"] = 0.5 * _scale(cand["expr_excess"].values) + 0.5 * _scale(
        cand["fc_excess"].values
    )
    cand = cand.sort_values("score", ascending=False).head(top_n)
    return cand


def _sizes_from_detection(detected_prop):
    """
    Map detection proportion [0..1] -> point size.
    Use sqrt-like scaling to reduce dynamic range.
    """
    p = np.clip(np.asarray(detected_prop).astype(float), 0.0, 1.0)
    p = np.power(p, POINTSIZE_POWER)
    return POINTSIZE_MIN + (POINTSIZE_MAX - POINTSIZE_MIN) * p


# ---------- Plotting ----------
def volcano_plot(df, title, outbase,
                 fdr_thr=0.05, lfc_thr=1.0, annotate_top=200):
    """
    df must have: gene, logFC, FDR, logCPM, pct_cells_detected
    """
    d = (
        df[["gene", "logFC", "FDR", "logCPM", "pct_cells_detected"]]
        .dropna(subset=["gene", "logFC", "FDR"])
        .copy()
    )
    d["mlog10FDR"] = -np.log10(d["FDR"].clip(lower=np.finfo(float).tiny))

    sig_up = (d["FDR"] < fdr_thr) & (d["logFC"] >= lfc_thr)
    sig_dn = (d["FDR"] < fdr_thr) & (d["logFC"] <= -lfc_thr)
    ns = ~(sig_up | sig_dn)

    # For subtitle
    n_up = int(sig_up.sum())
    n_dn = int(sig_dn.sum())

    sizes = _sizes_from_detection(d["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=220)
    ax.scatter(
        d.loc[ns, "logFC"],
        d.loc[ns, "mlog10FDR"],
        s=sizes[ns],
        c="#c7c7c7",
        alpha=0.45,
        linewidths=0,
        label="NS",
    )
    ax.scatter(
        d.loc[sig_dn, "logFC"],
        d.loc[sig_dn, "mlog10FDR"],
        s=sizes[sig_dn],
        c="#3b82f6",
        alpha=0.85,
        linewidths=0,
        label="Down",
    )
    ax.scatter(
        d.loc[sig_up, "logFC"],
        d.loc[sig_up, "mlog10FDR"],
        s=sizes[sig_up],
        c="#ef4444",
        alpha=0.85,
        linewidths=0,
        label="Up",
    )

    ax.axvline(lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axvline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-np.log10(fdr_thr), color="black", ls="--", lw=1, alpha=0.5)

    ax.set_xlabel(f"log2 fold-change ({COND_ALT} − {COND_REF})")
    ax.set_ylabel("-log10(FDR)")

    subtitle = f"Up: {n_up}   Down: {n_dn}   (FDR < {fdr_thr}, |log2FC| ≥ {lfc_thr})"
    ax.set_title(f"{title}\n{subtitle}", fontsize=9)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # Gene labels
    texts = []
    if ANNOTATION_MODE == "smart":
        d_sig = d.loc[sig_up | sig_dn, ["gene", "logFC", "logCPM", "mlog10FDR"]].copy()
        top = _smart_annot_rank(d_sig, annotate_top)
    else:
        d_sig = d.loc[
            sig_up | sig_dn, ["gene", "logFC", "FDR", "logCPM", "mlog10FDR"]
        ].copy()
        if not d_sig.empty:
            d_sig["abs_logFC"] = d_sig["logFC"].abs()
            top = (
                d_sig.sort_values(["FDR", "abs_logFC"], ascending=[True, False])
                .head(annotate_top)
            )
        else:
            top = d_sig

    if not top.empty:
        y_lookup = d.set_index("gene")["mlog10FDR"].to_dict()
        for _, r in top.iterrows():
            x = r["logFC"]
            y = y_lookup.get(r["gene"], r.get("mlog10FDR", np.nan))
            texts.append(
                ax.text(x, y, r["gene"], fontsize=6, ha="center", va="bottom")
            )
        adjust_text(
            texts,
            ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(
                arrowstyle="-", lw=0.5, color="black", alpha=0.6
            ),
        )

    fig.tight_layout()
    fig.savefig(outbase + ".pdf")
    fig.savefig(outbase + ".png", dpi=220)
    plt.close(fig)


def ma_plot(df, title, outbase, lfc_thr=1.0, annotate_top=200):
    # df must have: gene, logCPM, logFC, pct_cells_detected
    x = df["logCPM"].values  # log2 CPM
    y = df["logFC"].values
    sizes = _sizes_from_detection(df["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=220)
    ax.scatter(x, y, s=sizes, alpha=0.65, linewidths=0, c="#a8a8a8")
    ax.axhline(lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel("mean expression (logCPM)")
    ax.set_ylabel(f"log2 fold-change ({COND_ALT} − {COND_REF})")
    ax.set_title(title, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # Gene labels
    texts = []
    if ANNOTATION_MODE == "smart":
        pick = df[["gene", "logCPM", "logFC"]].copy()
        top = _smart_annot_rank(pick, annotate_top)
    else:
        pick = df[["gene", "logCPM", "logFC"]].copy()
        if not pick.empty:
            pick["abs_logFC"] = pick["logFC"].abs()
            top = pick.sort_values("abs_logFC", ascending=False).head(annotate_top)
        else:
            top = pick

    if not top.empty:
        lut_x = df.set_index("gene")["logCPM"].to_dict()
        lut_y = df.set_index("gene")["logFC"].to_dict()
        for _, r in top.iterrows():
            gx = lut_x.get(r["gene"], r["logCPM"])
            gy = lut_y.get(r["gene"], r["logFC"])
            texts.append(
                ax.text(gx, gy, r["gene"], fontsize=6, ha="center", va="bottom")
            )
        adjust_text(
            texts,
            ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(
                arrowstyle="-", lw=0.5, color="black", alpha=0.6
            ),
        )

    fig.tight_layout()
    fig.savefig(outbase + ".pdf")
    fig.savefig(outbase + ".png", dpi=220)
    plt.close(fig)


def quick_pca_plot(counts_df, meta_df, title, out_png):
    try:
        from sklearn.decomposition import PCA
    except ImportError:
        return
    X = counts_df.T.values
    lib = X.sum(1, keepdims=True)
    Xnorm = np.log1p(1e6 * X / np.maximum(lib, 1))
    pc = PCA(n_components=2).fit_transform(Xnorm)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    codes = pd.Categorical(meta_df["condition"]).codes
    ax.scatter(pc[:, 0], pc[:, 1], s=40, c=codes)
    for i, t in enumerate(meta_df.index):
        ax.text(pc[i, 0], pc[i, 1], t, fontsize=6, alpha=0.6)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


# ------------------- R-side via rpy2 -------------------
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

edgeR = importr("edgeR")
limma = importr("limma")
deseq2 = importr("DESeq2")

# Full edgeR pipeline for binary contrast: ALT - REF = cluster - rest
ro.r(f"""
run_edgeR_binary <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  if (nlevels(group) != 2) stop("Need exactly 2 groups: {COND_REF} and {COND_ALT}")
  library(edgeR); library(limma)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  y <- DGEList(counts=counts)
  keep <- filterByExpr(y, group=group)
  y <- y[keep,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  contrast <- makeContrasts(
    contrasts="{COND_ALT}-{COND_REF}",
    levels=design
  )
  qlf <- glmQLFTest(fit, contrast=contrast)
  tt  <- topTags(qlf, n=Inf)$table
  tt$gene <- rownames(tt)
  tt
}}
""")

# Effect-size-only path (no p-values), forced levels ALT - REF
ro.r(f"""
edgeR_effect_only_binary <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  if (nlevels(group) != 2) stop("Need exactly 2 groups: {COND_REF}, {COND_ALT}")
  y <- edgeR::DGEList(counts=counts)
  y <- edgeR::calcNormFactors(y)
  lcpms <- edgeR::cpm(y, log=TRUE, prior.count=0.5)
  mu_ref <- rowMeans(lcpms[, group=="{COND_REF}", drop=FALSE])
  mu_alt <- rowMeans(lcpms[, group=="{COND_ALT}", drop=FALSE])
  out  <- data.frame(
    gene   = rownames(counts),
    logFC  = mu_alt - mu_ref,           # cluster − rest
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
run_DESeq2_binary <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  suppressPackageStartupMessages(library(DESeq2))
  coldata <- data.frame(group=group)
  rownames(coldata) <- colnames(counts)
  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData=coldata,
                                design=~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("group", "{COND_ALT}", "{COND_REF}"))
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


def run_edgeR(counts_df, meta_df):
    # Ensure columns of counts_df are in the same order as rows of meta_df
    counts_df = counts_df.loc[:, meta_df.index]
    group = meta_df["condition"].astype(str).values

    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    rfun = ro.globalenv["run_edgeR_binary"]

    with localconverter(conv):
        dfR = rfun(counts_df, group)
        df = ro.conversion.rpy2py(dfR)

    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df


def effect_only_table(counts_df, meta_df):
    counts_df = counts_df.loc[:, meta_df.index]
    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    with localconverter(conv):
        dfR = ro.globalenv["edgeR_effect_only_binary"](
            counts_df, meta_df["condition"].astype(str).values
        )
        df = ro.conversion.rpy2py(dfR)
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df[["gene", "logFC", "logCPM", "PValue", "FDR"]]


def run_DESeq2(counts_df, meta_df):
    # Ensure columns of counts_df are in the same order as rows of meta_df
    counts_df = counts_df.loc[:, meta_df.index]
    group = meta_df["condition"].astype(str).values

    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    rfun = ro.globalenv["run_DESeq2_binary"]

    with localconverter(conv):
        dfR = rfun(counts_df, group)
        df = ro.conversion.rpy2py(dfR)

    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df


# ------------------- Utility -------------------
def slugify(name):
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return s or "cluster"


# ------------------- Main Loop -------------------
def main():
    print(f"Reading AnnData from: {H5AD_PATH}")
    A0 = sc.read_h5ad(H5AD_PATH)
    if CELLTYPE_KEY not in A0.obs:
        raise KeyError(f"obs['{CELLTYPE_KEY}'] missing in AnnData.")

    print(f"Original shape: {A0.n_obs} cells × {A0.n_vars} genes")
    print("Cell-type counts (raw):")
    print(A0.obs[CELLTYPE_KEY].value_counts())

    A = filter_genes_drop_malat1_mt(A0)
    print(f"After gene filtering: {A.n_obs} cells × {A.n_vars} genes")

    counts_layer = pick_counts_layer(A)
    print(f"Using counts from layer = {counts_layer or 'X (integers)'}")

    # Detection stats over all cells (for point sizes)
    det_df, n_cells_total = _gene_detection(A, counts_layer)
    print(f"Computed detection stats over {n_cells_total} cells.")

    vc = A.obs[CELLTYPE_KEY].value_counts()
    print("\nWill process clusters with ≥ "
          f"{PB_MIN_CELLS_PER_GROUP} cells:\n{vc}")

    for ct, n_ct in vc.items():
        if n_ct < PB_MIN_CELLS_PER_GROUP:
            print(f"\n[SKIP] {ct}: only {n_ct} cells "
                  f"(< {PB_MIN_CELLS_PER_GROUP}).")
            continue

        print(f"\n=== {ct} ===  (n = {n_ct} cells)")

        cond_key = "__marker_cond__"
        cond = np.where(A.obs[CELLTYPE_KEY] == ct, COND_ALT, COND_REF)

        A_ct = A.copy()
        A_ct.obs[cond_key] = cond

        counts_df, meta_df = make_pseudobulk_binary(
            A_ct,
            counts_layer,
            cond_key=cond_key,
            min_cells=PB_MIN_CELLS_PER_CONDITION,
            reps_per_cond=REPLICATES_PER_CONDITION,
        )
        if counts_df is None:
            print(f"  -> Not enough cells for robust pseudobulk in {ct}; skipping.")
            continue

        ct_tag = slugify(ct)

        quick_pca_plot(
            counts_df,
            meta_df,
            title=f"{ct} vs rest pseudobulk PCA (log CPM)",
            out_png=os.path.join(OUTDIR, f"{ct_tag}_pseudobulk_PCA.png"),
        )

        base_edgeR = os.path.join(OUTDIR, f"{ct_tag}_pb_edgeR")
        base_DESeq2 = os.path.join(OUTDIR, f"{ct_tag}_pb_DESeq2")

        if have_replicates(meta_df):
            # ---------- edgeR QLF ----------
            title_edgeR = f"{ct} — {COND_ALT} vs {COND_REF} (edgeR QLF)"
            df_edgeR = run_edgeR(counts_df, meta_df)
            df_edgeR = df_edgeR.merge(
                det_df, left_on="gene", right_index=True, how="left"
            )
            keep_cols = [
                c
                for c in [
                    "gene",
                    "logFC",
                    "logCPM",
                    "F",
                    "PValue",
                    "FDR",
                    "n_cells_detected",
                    "pct_cells_detected",
                ]
                if c in df_edgeR.columns
            ]
            df_edgeR = df_edgeR[keep_cols]
            df_edgeR.to_csv(f"{base_edgeR}.csv", index=False)
            volcano_plot(
                df_edgeR,
                title=title_edgeR,
                outbase=f"{base_edgeR}__volcano",
                fdr_thr=FDR_THR,
                lfc_thr=LFC_THR,
                annotate_top=TOP_ANNOTATE,
            )

            # ---------- DESeq2 ----------
            title_DE = f"{ct} — {COND_ALT} vs {COND_REF} (DESeq2)"
            df_DE = run_DESeq2(counts_df, meta_df)
            df_DE = df_DE.merge(
                det_df, left_on="gene", right_index=True, how="left"
            )
            keep_cols_DE = [
                c
                for c in [
                    "gene",
                    "logFC",
                    "logCPM",
                    "PValue",
                    "FDR",
                    "n_cells_detected",
                    "pct_cells_detected",
                ]
                if c in df_DE.columns
            ]
            df_DE = df_DE[keep_cols_DE]
            df_DE.to_csv(f"{base_DESeq2}.csv", index=False)
            volcano_plot(
                df_DE,
                title=title_DE,
                outbase=f"{base_DESeq2}__volcano",
                fdr_thr=FDR_THR,
                lfc_thr=LFC_THR,
                annotate_top=TOP_ANNOTATE,
            )

        else:
            print(
                f"  -> {ct}: not enough pseudobulk replicates per group "
                "for edgeR/DESeq2; reporting effect sizes only."
            )
            df = effect_only_table(counts_df, meta_df)
            df = df.merge(det_df, left_on="gene", right_index=True, how="left")
            base = base_edgeR
            df.to_csv(base + "__effect_only.csv", index=False)
            ma_plot(
                df,
                title=f"{ct} — {COND_ALT} vs {COND_REF}",
                outbase=base + "__effect_only_MA",
                lfc_thr=LFC_THR,
                annotate_top=TOP_ANNOTATE,
            )

    print("\nDone. Outputs in:", OUTDIR)


if __name__ == "__main__":
    main()











import os
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ------------------- Configuration -------------------
OUTDIR = "figures_cluster_markers_pseudobulk"
os.makedirs(OUTDIR, exist_ok=True)

# <<< CHANGE THIS TO YOUR H5AD FILE >>>
H5AD_PATH = "./MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"

# Cell-type / cluster annotation
CELLTYPE_KEY = "leiden_merged_type"

adata = sc.read(H5AD_PATH)

# ----------------- CONFIG -----------------
celltype_key       = "leiden_merged_type"
dge_dir            = "figures_cluster_markers_pseudobulk"  # where <slug>_pb_DESeq2.csv live
top_n              = 2          # top genes per cell type BEFORE dedup
max_fdr            = 0.001       # FDR cutoff
min_abs_lfc        = 0.5        # |log2FC| cutoff
min_logCPM         = 0.5        # mean expression cutoff (DESeq2 logCPM)
min_pct_detected   = 0.02       # fraction of cells (0–1) with >0 counts

def slugify(name: str) -> str:
    """Must match how you named <cluster>_pb_DESeq2.csv."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return s or "cluster"


# ensure categorical order is defined
adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
clusters = list(adata.obs[celltype_key].cat.categories)

var_names = {}       # dict for scanpy.pl.dotplot: {cluster_label: [genes]}
used_genes = set()   # global set to de-duplicate genes across clusters

for ct in clusters:
    slug = slugify(ct)
    csv_path = os.path.join(dge_dir, f"{slug}_pb_DESeq2.csv")
    if not os.path.exists(csv_path):
        print(f"[skip] {ct}: {csv_path} not found")
        continue

    df = pd.read_csv(csv_path)
    if "gene" not in df.columns:
        print(f"[skip] {ct}: no 'gene' column in {csv_path}")
        continue

    # numeric columns + convenience abs_logFC
    for c in ["logFC", "logCPM", "PValue", "FDR", "pct_cells_detected"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["logFC"])
    df["abs_logFC"] = df["logFC"].abs()

    # ----------------- FILTERS (toggles) -----------------
    if max_fdr is not None and "FDR" in df.columns:
        df = df[df["FDR"] <= max_fdr]

    if min_abs_lfc is not None:
        df = df[df["abs_logFC"] >= min_abs_lfc]

    if min_logCPM is not None and "logCPM" in df.columns:
        df = df[df["logCPM"] >= min_logCPM]

    if min_pct_detected is not None and "pct_cells_detected" in df.columns:
        df = df[df["pct_cells_detected"] >= min_pct_detected]

    if df.empty:
        print(f"[skip] {ct}: all genes filtered out")
        continue

    # ----------------- RANKING (per cluster) -----------------
    # smaller FDR is better, larger abs_logFC is better
    # ranking within cluster: FDR then |logFC|
    if "FDR" not in df.columns:
        df["FDR"] = np.nan
    
    df = df.sort_values(
        ["FDR", "abs_logFC"],
        ascending=[True, False],
    )

 # pick TOP_N genes per cell type, scanning deeper to avoid global duplicates
    genes_for_ct = []
    dropped_dups = []

    for g in df["gene"].astype(str).tolist():
        if g not in adata.var_names:
            continue
        if g in used_genes:
            dropped_dups.append(g)
            continue
        genes_for_ct.append(g)
        used_genes.add(g)
        if len(genes_for_ct) == top_n:
            break

    if dropped_dups:
        print(f"[dedup] {ct}: skipped already-used genes (showing up to 10): {dropped_dups[:10]}")

    if len(genes_for_ct) < top_n:
        raise RuntimeError(
            f"{ct}: only found {len(genes_for_ct)}/{top_n} unique genes after scanning DESeq2 list. "
            f"Consider relaxing filters (max_fdr/min_abs_lfc/min_logCPM/min_pct_detected) or disabling global dedup."
        )

    var_names[ct] = genes_for_ct
    print(f"[OK] {ct}: using {len(genes_for_ct)} genes -> {genes_for_ct}")


print("\nFinal var_names dict going into scanpy.pl.dotplot:")
for ct, genes in var_names.items():
    print(f"  {ct}: {genes}")

sc.tl.dendrogram(
    adata,
    groupby=celltype_key,
    use_rep="X_pca",   # or None to use .X, or "X_umap" etc. if you prefer
)

# ---- the actual dotplot ----
sc.pl.dotplot(
    adata,
    var_names,               # dict: {cluster: [genes]}
    groupby=celltype_key,    # rows = your cell types
    standard_scale="var",    # scale per gene across groups (nice for markers)
    swap_axes=False,         # genes on x, clusters on y; flip if you prefer
    figsize=(8, 4),
    cmap="Reds",             # or e.g. "viridis"
    dendrogram=True,         # <---- this turns the tree on
    save="AllCellTypes_PseudobulkDESeq2.pdf",               # set to ".pdf" / ".png" to auto-save
)


































#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse as sp
import snapatac2 as snap
from matplotlib.gridspec import GridSpec

# ----------------- CONFIG -----------------

RNA_PATH = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"

AC_PATH  = "DNA_merged/MergedTonsils_H3K27ac_merged_processed_barcodeRewritten_withRNAlabels.h5ad"
ME3_PATH = "DNA_merged/MergedTonsils_H3K27me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad"
K9_PATH  = "DNA_merged/MergedTonsils_H3K9me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad"

DGE_DIR  = "figures_cluster_markers_pseudobulk"



# ----------------- CONFIG -----------------
celltype_key       = "leiden_merged_type"
dge_dir            = "figures_cluster_markers_pseudobulk"  # where <slug>_pb_DESeq2.csv live
top_n              = 2          # top genes per cell type BEFORE dedup
max_fdr            = 0.001       # FDR cutoff
min_abs_lfc        = 0.5        # |log2FC| cutoff
min_logCPM         = 0.5        # mean expression cutoff (DESeq2 logCPM)
min_pct_detected   = 0.02       # fraction of cells (0–1) with >0 counts

PROMOTER_WINDOW = 2000

AGG_FUNC   = "mean"          # "mean" or "median"
NORM_MODE  = "zscore_gene"   # "none", "minmax_gene", "zscore_gene"

OUTDIR = f"heatmaps_top{top_n}_markers_RNA_H3K27ac_H3K27me3_H3K9me3_groupmeans"
os.makedirs(OUTDIR, exist_ok=True)

try:
    GENE_ANNO = snap.genome.hg38
except Exception:
    GENE_ANNO = None


# ----------------- HELPERS -----------------

def slugify(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return s or "cluster"


def _resolve_anno_for_snap(gene_anno):
    if gene_anno is None:
        raise ValueError("GENE_ANNO is None. Set it to snap.genome.hg38 or a GTF/GFF path.")
    return getattr(gene_anno, "annotation", gene_anno)


def make_gene_activity(adata_mark, gene_anno, window_bp, layer_name="log1p"):
    ann = _resolve_anno_for_snap(gene_anno)
    ga = snap.pp.make_gene_matrix(
        adata_mark,
        ann,
        upstream=window_bp,
        downstream=window_bp,
        include_gene_body=True,
        id_type="gene",
    )

    for col in ["sample", "tonsil_id", "leiden_merged_type", "leiden"]:
        if col in adata_mark.obs:
            ga.obs[col] = adata_mark.obs[col].copy()

    X = ga.X
    if sp.issparse(X):
        Xlog = X.copy()
        Xlog.data = np.log1p(Xlog.data)
    else:
        Xlog = np.log1p(X)
    ga.layers[layer_name] = Xlog
    return ga


def compute_group_means(ad, genes, groupby, layer=None, group_order=None, agg_func="mean"):
    df = sc.get.obs_df(ad, keys=genes, layer=layer)
    df[groupby] = ad.obs[groupby].astype(str).values

    if agg_func == "median":
        g = df.groupby(groupby).median(numeric_only=True)
    else:
        g = df.groupby(groupby).mean(numeric_only=True)

    if group_order is not None:
        g = g.reindex(group_order)

    g = g[genes]
    g = g.fillna(0.0)
    return g


def _normalize_matrix(M, mode="none", z_clip=1.5):
    if mode == "none":
        return M, "mean expression"

    if mode == "minmax_gene":
        M_scaled = np.zeros_like(M, dtype=float)
        for i in range(M.shape[0]):
            row = M[i, :]
            rmin = float(np.min(row))
            rmax = float(np.max(row))
            if rmax > rmin:
                M_scaled[i, :] = (row - rmin) / (rmax - rmin)
            else:
                M_scaled[i, :] = 0.0
        return M_scaled, "per-gene min–max scaled mean"

    if mode == "zscore_gene":
        M_scaled = np.zeros_like(M, dtype=float)
        for i in range(M.shape[0]):
            row = M[i, :]
            mu = float(np.mean(row))
            sd = float(np.std(row))
            if sd > 0:
                z = (row - mu) / sd
                if z_clip is not None:
                    z = np.clip(z, -z_clip, z_clip)
                M_scaled[i, :] = z
            else:
                M_scaled[i, :] = 0.0
        return M_scaled, f"per-gene z-score (clip ±{z_clip})"

    return M, "mean expression"


def plot_group_mean_heatmap(
    mean_df,
    title,
    fname,
    cmap="viridis",
    norm_mode="none",
    var_group_positions=None,
    var_group_labels=None,
    x_label=None,
    agg_label="mean",
):
    """
    mean_df: DataFrame
      index   = groups / cell types
      columns = genes

    Plots:
      x-axis = genes
      y-axis = groups / cell types
    """
    groups = list(mean_df.index)
    genes = list(mean_df.columns)

    # matrix for plotting: groups × genes
    M = mean_df.to_numpy()

    # normalize per gene
    M_scaled, norm_label = _normalize_matrix(M.T, mode=norm_mode)
    M_scaled = M_scaled.T

    n_groups, n_genes = M_scaled.shape

    fig_w = max(6, 0.35 * n_genes + 2.5)
    fig_h = max(4.5, 0.45 * n_groups + 2.5)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
    gs = GridSpec(nrows=1, ncols=2, width_ratios=[20, 1], wspace=0.15, figure=fig)

    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(M_scaled, aspect="auto", interpolation="nearest", cmap=cmap)

    # axes / labels
    ax.set_xticks(np.arange(n_genes))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(n_groups))
    ax.set_yticklabels(groups, fontsize=8)

    ax.set_xlabel("marker genes", fontsize=9)
    if x_label is not None:
        ax.set_ylabel(x_label, fontsize=9)
    else:
        ax.set_ylabel("cell type", fontsize=9)

    ax.set_title(title, fontsize=10)

    # vertical lines between gene blocks
    if var_group_positions is not None and len(var_group_positions) > 1:
        for (_, end) in var_group_positions[:-1]:
            ax.axvline(end + 0.5, color="white", linewidth=0.5)

    # optional top labels for gene blocks
    if var_group_positions is not None and var_group_labels is not None:
        centers = [0.5 * (s + e) for (s, e) in var_group_positions]
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(centers)
        ax2.set_xticklabels(var_group_labels, rotation=45, ha="left", fontsize=8)
        ax2.set_xlabel("cluster of origin", fontsize=9)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(f"{norm_label} (agg: {agg_label})", fontsize=8)

    pdf = os.path.join(OUTDIR, f"{fname}.pdf")
    png = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[saved] {pdf}")

# ----------------- LOAD DATA -----------------

print(f"[INFO] Loading RNA AnnData from: {RNA_PATH}")
adata_rna = sc.read_h5ad(RNA_PATH)

print(f"[INFO] Loading H3K27ac from:  {AC_PATH}")
print(f"[INFO] Loading H3K27me3 from: {ME3_PATH}")
print(f"[INFO] Loading H3K9me3 from:  {K9_PATH}")

adata_ac_raw  = sc.read_h5ad(AC_PATH)
adata_me3_raw = sc.read_h5ad(ME3_PATH)
adata_k9_raw  = sc.read_h5ad(K9_PATH)

if celltype_key not in adata_rna.obs:
    raise KeyError(f"obs['{celltype_key}'] is missing from RNA AnnData.")

# Fix / enforce categorical ordering from RNA
adata_rna.obs[celltype_key] = adata_rna.obs[celltype_key].astype("category")
clusters = list(adata_rna.obs[celltype_key].cat.categories)
print("[INFO] cell types (order):", clusters)

for mod_name, mod in [
    ("H3K27ac", adata_ac_raw),
    ("H3K27me3", adata_me3_raw),
    ("H3K9me3", adata_k9_raw),
]:
    if celltype_key not in mod.obs:
        raise KeyError(f"{celltype_key} missing from {mod_name} object")
    mod.obs[celltype_key] = mod.obs[celltype_key].astype("category")
    mod.obs[celltype_key] = mod.obs[celltype_key].cat.set_categories(clusters)


# ----------------- TOP MARKERS FROM RNA DGE -----------------

var_names = {}
used_genes = set()

for ct in clusters:
    slug = slugify(ct)
    csv_path = os.path.join(DGE_DIR, f"{slug}_pb_DESeq2.csv")

    if not os.path.exists(csv_path):
        print(f"[skip] {ct}: {csv_path} not found")
        continue

    df = pd.read_csv(csv_path)
    if "gene" not in df.columns:
        print(f"[skip] {ct}: no 'gene' column in {csv_path}")
        continue

    for c in ["logFC", "logCPM", "PValue", "FDR", "pct_cells_detected"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["logFC"])
    df["abs_logFC"] = df["logFC"].abs()

    if max_fdr is not None and "FDR" in df.columns:
        df = df[df["FDR"] <= max_fdr]
    if min_abs_lfc is not None:
        df = df[df["abs_logFC"] >= min_abs_lfc]
    if min_logCPM is not None and "logCPM" in df.columns:
        df = df[df["logCPM"] >= min_logCPM]
    if min_pct_detected is not None and "pct_cells_detected" in df.columns:
        df = df[df["pct_cells_detected"] >= min_pct_detected]

    if df.empty:
        print(f"[skip] {ct}: all genes filtered out after thresholds")
        continue

    if "FDR" not in df.columns:
        df["FDR"] = np.nan

    df = df.sort_values(["FDR", "abs_logFC"], ascending=[True, False])

    genes_for_ct = []
    dropped_dups = []

    for g in df["gene"].astype(str).tolist():
        if g not in adata_rna.var_names:
            continue
        if g in used_genes:
            dropped_dups.append(g)
            continue
        genes_for_ct.append(g)
        used_genes.add(g)
        if len(genes_for_ct) == top_n:
            break

    if dropped_dups:
        print(f"[dedup] {ct}: skipped already-used genes (showing up to 10): {dropped_dups[:10]}")

    if len(genes_for_ct) < top_n:
        raise RuntimeError(
            f"{ct}: only found {len(genes_for_ct)}/{top_n} unique genes after scanning DESeq2 list."
        )

    var_names[ct] = genes_for_ct
    print(f"[OK] {ct}: using {len(genes_for_ct)} genes -> {genes_for_ct}")

print("\nFinal marker list:")
for ct, genes in var_names.items():
    print(f"  {ct}: {genes}")


# ----------------- BUILD GENE-ACTIVITY MATRICES -----------------

print("\n[INFO] Building gene-activity matrices for H3K27ac / H3K27me3 / H3K9me3 ...")

ga_ac  = make_gene_activity(adata_ac_raw,  GENE_ANNO, PROMOTER_WINDOW)
ga_me3 = make_gene_activity(adata_me3_raw, GENE_ANNO, PROMOTER_WINDOW)
ga_k9  = make_gene_activity(adata_k9_raw,  GENE_ANNO, PROMOTER_WINDOW)

for ga in (ga_ac, ga_me3, ga_k9):
    if celltype_key not in ga.obs:
        raise KeyError(f"{celltype_key} missing after gene activity creation")
    ga.obs[celltype_key] = ga.obs[celltype_key].astype("category")
    ga.obs[celltype_key] = ga.obs[celltype_key].cat.set_categories(clusters)


# ----------------- SHARED GENE LIST ACROSS ALL 4 MATRICES -----------------

genes_flat = []
for ct in clusters:
    if ct not in var_names:
        continue
    genes_flat.extend(var_names[ct])

seen = set()
genes_flat_unique = []
for g in genes_flat:
    if g not in seen:
        seen.add(g)
        genes_flat_unique.append(g)

print("\nTop-marker genes (before multi-modality intersection):")
print(genes_flat_unique)

var_rna = set(adata_rna.var_names)
var_ac  = set(ga_ac.var_names)
var_me3 = set(ga_me3.var_names)
var_k9  = set(ga_k9.var_names)

genes_shared = []
missing_ac = []
missing_me3 = []
missing_k9 = []

for g in genes_flat_unique:
    in_rna = g in var_rna
    in_ac  = g in var_ac
    in_me3 = g in var_me3
    in_k9  = g in var_k9

    if in_rna and in_ac and in_me3 and in_k9:
        genes_shared.append(g)
    else:
        if not in_ac:
            missing_ac.append(g)
        if not in_me3:
            missing_me3.append(g)
        if not in_k9:
            missing_k9.append(g)

if missing_ac:
    print(f"[WARN] genes missing in H3K27ac GA: {missing_ac}")
if missing_me3:
    print(f"[WARN] genes missing in H3K27me3 GA: {missing_me3}")
if missing_k9:
    print(f"[WARN] genes missing in H3K9me3 GA: {missing_k9}")

if not genes_shared:
    raise RuntimeError("No marker genes are present in all of RNA / H3K27ac / H3K27me3 / H3K9me3.")

print("\nGenes present in ALL modalities:")
print(genes_shared)

filtered_var_names = {}
for ct in clusters:
    if ct not in var_names:
        continue
    kept = [g for g in var_names[ct] if g in genes_shared]
    if kept:
        filtered_var_names[ct] = kept

genes_ordered = []
var_group_positions = []
var_group_labels = []

start = 0
for ct in clusters:
    genes_ct = filtered_var_names.get(ct, [])
    if not genes_ct:
        continue
    genes_ordered.extend(genes_ct)
    end = start + len(genes_ct) - 1
    var_group_positions.append((start, end))
    var_group_labels.append(ct)
    start = end + 1

if not genes_ordered:
    raise RuntimeError("No genes left after intersecting with all modality matrices.")

print("\nFinal gene order for heatmaps:")
for ct, (g_start, g_end) in zip(var_group_labels, var_group_positions):
    print(f"  {ct}: {genes_ordered[g_start:g_end+1]}")


# ----------------- LAYERS -----------------

RNA_LAYER = None
AC_LAYER  = "log1p"
ME3_LAYER = "log1p"
K9_LAYER  = "log1p"


# ----------------- GROUP-MEAN MATRICES -----------------

rna_mean = compute_group_means(
    ad=adata_rna,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=RNA_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)

ac_mean = compute_group_means(
    ad=ga_ac,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=AC_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)

me3_mean = compute_group_means(
    ad=ga_me3,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=ME3_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)

k9_mean = compute_group_means(
    ad=ga_k9,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=K9_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)

print("[INFO] rna_mean shape:", rna_mean.shape)
print("[INFO] ac_mean shape:", ac_mean.shape)
print("[INFO] me3_mean shape:", me3_mean.shape)
print("[INFO] k9_mean shape:", k9_mean.shape)


# ----------------- HEATMAPS -----------------

groupby_label = celltype_key

plot_group_mean_heatmap(
    mean_df=rna_mean,
    title=f"RNA: {AGG_FUNC} expression per cell type (top{top_n} DESeq2 markers)",
    fname=f"groupmean_RNA_top{top_n}_markers_by_celltype",
    cmap="viridis",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

plot_group_mean_heatmap(
    mean_df=ac_mean,
    title=f"H3K27ac gene activity (±{PROMOTER_WINDOW//1000}kb): {AGG_FUNC} per cell type",
    fname=f"groupmean_H3K27ac_top{top_n}_markers_by_celltype",
    cmap="magma",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

plot_group_mean_heatmap(
    mean_df=me3_mean,
    title=f"H3K27me3 gene activity (±{PROMOTER_WINDOW//1000}kb): {AGG_FUNC} per cell type",
    fname=f"groupmean_H3K27me3_top{top_n}_markers_by_celltype",
    cmap="magma",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

plot_group_mean_heatmap(
    mean_df=k9_mean,
    title=f"H3K9me3 gene activity (±{PROMOTER_WINDOW//1000}kb): {AGG_FUNC} per cell type",
    fname=f"groupmean_H3K9me3_top{top_n}_markers_by_celltype",
    cmap="magma",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

print(f"\nDone. Group-mean heatmaps written to: {OUTDIR}")

sc.pl.dotplot(
    adata_rna,
    genes_ordered,                 # flat list of genes
    groupby=celltype_key,
    standard_scale="var",
    swap_axes=False,
    figsize=(10, 4),
    cmap="Reds",
    dendrogram=False,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    var_group_rotation=45,
    save=f"_RNA_AllCellTypes_PseudobulkDESeq2_top{top_n}.pdf",
)

