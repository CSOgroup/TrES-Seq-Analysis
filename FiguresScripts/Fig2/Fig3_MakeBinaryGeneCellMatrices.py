#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_binary_gene_cell_matrices_COMMON3.py

Build 0/1 matrices (genes x cells) for:
  - RNA expression (Expr)
  - H3K27ac gene-activity (±5 kb around gene body)
  - H3K27me3 gene-activity (±5 kb around gene body)

This version RESTRICTS to genes present in ALL THREE modalities (COMMON3),
so the three outputs (RNA / AC / ME3) have identical shapes and gene order.

No external GTF: RNA gene list from adata_rna.var; GA genes from SnapATAC2 using snap.genome.hg38.

Outputs:
  - binary_Expr_*__COMMON3.tsv.gz
  - binary_H3K27ac_pm1kb_*__COMMON3.tsv.gz
  - binary_H3K27me3_pm1kb_*__COMMON3.tsv.gz
  - (and MARKERS_COMMON3 if rGG markers exist)
  - COMMON3_genes_order.txt
  - COMMON3_cells_order.txt
  - MANIFEST_binary_matrices.tsv

Requirements: scanpy, anndata, snapatac2, pandas, numpy, scipy
"""

import os, gzip
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import scipy.sparse as sp
from typing import Dict, List

# ------------------ CONFIG ------------------
RNA_PATH = "AllCellLines_rna_merged_processed_Filtered.h5ad"
AC_PATH  = "AllCellLines_H3K27ac_merged_processed.h5ad"
ME3_PATH = "AllCellLines_H3K27me3_merged_processed.h5ad"

OUTDIR = "CellLine_BinaryMatrices"
os.makedirs(OUTDIR, exist_ok=True)

# rank_genes_groups keys to pull markers (as in your RNA marker script)
RGG_KEY       = "CL_rgg"
RGG_GROUPBY   = "CL"
RGG_TOP_N     = 50
RGG_ALPHA     = 0.05
RGG_MIN_LFC   = 0.25
RGG_UNIQUE    = True
RGG_USE_SCORE = True

# Gene-activity window
GA_FLANK = 5_000

# For output naming
GA_FLANK_TAG = f"pm{GA_FLANK // 1000}kb" if (GA_FLANK % 1000 == 0) else f"pm{GA_FLANK}bp"
AC_MODALITY  = f"H3K27ac_{GA_FLANK_TAG}"
ME3_MODALITY = f"H3K27me3_{GA_FLANK_TAG}"

# CPM scaling for GA (counts per 1e5)
GA_CPM_SCALE = 1e5

# ------------------ LOAD ------------------
adata_rna = sc.read_h5ad(RNA_PATH)
adata_ac  = sc.read_h5ad(AC_PATH)
adata_me3 = sc.read_h5ad(ME3_PATH)


def _strip_suffix_after_last_underscore(adata, col_name="barcode_orig"):
    """
    Remove everything after the last '_' in obs_names.
    E.g. 'AAAC-1_GM12878_ac' -> 'AAAC-1'
    """
    adata = adata.copy()
    idx = pd.Index(adata.obs_names)
    adata.obs[col_name] = idx  # keep original if you ever need it
    new_idx = idx.str.rsplit("_", n=1).str[0]
    adata.obs_names = new_idx
    return adata

adata_ac = _strip_suffix_after_last_underscore(adata_ac)
adata_me3 = _strip_suffix_after_last_underscore(adata_me3)

# Keep shared cells; preserve RNA order
common_cells = adata_rna.obs_names.intersection(adata_ac.obs_names).intersection(adata_me3.obs_names)
adata_rna = adata_rna[common_cells].copy()
adata_ac  = adata_ac[common_cells].copy()
adata_me3 = adata_me3[common_cells].copy()

# ------------------ HELPERS ------------------
def _normalize_symbols(s: pd.Series) -> pd.Series:
    """Turn 'A,B' -> 'A', strip spaces."""
    return s.astype(str).str.split(',').str[0].str.strip()

def _infer_symbol_col(adata) -> str:
    for c in ['gene_symbol','gene_symbols','symbol','SYMBOL','feature_name','GeneSymbol']:
        if c in adata.var.columns:
            return c
    return None

def _cell_labels_with_sample(adata) -> List[str]:
    """Return 'sample|barcode' when possible, else barcode."""
    if "sample" in adata.obs:
        return [f"{str(adata.obs.loc[c,'sample'])}|{c}" for c in adata.obs_names]
    return list(adata.obs_names)

def _ensure_log1p_layer(adata, prefer='counts', layer_name='log1p'):
    if layer_name in adata.layers:
        return layer_name
    Xsrc = adata.layers[prefer] if prefer in adata.layers else adata.X
    if sp.issparse(Xsrc):
        Xlog = Xsrc.copy(); Xlog.data = np.log1p(Xlog.data)
    else:
        Xlog = np.log1p(Xsrc)
    adata.layers[layer_name] = Xlog
    return layer_name

def _dedup_index(df: pd.DataFrame, how='sum') -> pd.DataFrame:
    """Aggregate duplicate gene symbols (RNA: sum; GA: max)."""
    if df.index.is_unique:
        return df
    dup_count = int(df.index.size - df.index.unique().size)
    print(f"[warn] deduplicating {dup_count} duplicate gene symbol rows via '{how}'")
    return df.groupby(level=0, sort=False).agg(how)

def _z_per_gene_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Row-wise z-score (per gene across cells), NaN-safe -> 0."""
    mu = df.mean(axis=1)
    sd = df.std(axis=1).replace(0, np.nan)
    z  = df.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)
    return z

def _cpm_df(df: pd.DataFrame, scale=GA_CPM_SCALE) -> pd.DataFrame:
    """Column-wise CPM-like normalization (genes x cells)."""
    col_sum = df.sum(axis=0).replace(0.0, np.nan)
    return df.div(col_sum, axis=1).fillna(0.0) * float(scale)

def write_tsv_gz(df: pd.DataFrame, path: str):
    with gzip.open(path, 'wt') as f:
        df.to_csv(f, sep='\t', header=True, index=True, lineterminator='\n')

# ------------------ MAKE GENE ACTIVITY (±1 kb, NO GTF) ------------------
try:
    GENE_ANNO = snap.genome.hg38
except Exception as e:
    raise RuntimeError("snap.genome.hg38 is required (no external GTF used).") from e

def make_gene_activity_tables(adata_mark, gene_anno, flank_bp):
    """
    SnapATAC2 GA (cell×gene) → transpose to genes×cells; also build log1p DF.
    """
    ga = snap.pp.make_gene_matrix(
        adata_mark,
        gene_anno,                # genome annotation object
        include_gene_body=True,
        upstream=flank_bp,
        downstream=flank_bp,
        id_type="gene",
    )
    ga = ga[adata_rna.obs_names, :]  # align cells

    X = ga.X.toarray() if sp.issparse(ga.X) else np.asarray(ga.X)
    counts_gxc = X.T if X.shape == (ga.n_obs, ga.n_vars) else X  # genes×cells

    Xlog = ga.layers.get("log1p", None)
    if Xlog is None:
        if sp.issparse(ga.X):
            Xlog = ga.X.copy(); Xlog.data = np.log1p(Xlog.data)
        else:
            Xlog = np.log1p(ga.X)
    Xlog = Xlog.toarray() if sp.issparse(Xlog) else np.asarray(Xlog)
    log1p_gxc = Xlog.T if Xlog.shape == (ga.n_obs, ga.n_vars) else Xlog

    gene_idx = _normalize_symbols(pd.Series(ga.var_names, index=ga.var_names)).values
    cell_cols = list(ga.obs_names)

    df_counts = pd.DataFrame(counts_gxc, index=gene_idx, columns=cell_cols)
    df_log1p  = pd.DataFrame(log1p_gxc,  index=gene_idx, columns=cell_cols)

    # GA duplicates: presence-like => MAX
    df_counts = _dedup_index(df_counts, how='max')
    df_log1p  = _dedup_index(df_log1p,  how='max')
    return df_counts, df_log1p

ga_ac_counts,  ga_ac_log1p  = make_gene_activity_tables(adata_ac,  GENE_ANNO, GA_FLANK)
ga_me3_counts, ga_me3_log1p = make_gene_activity_tables(adata_me3, GENE_ANNO, GA_FLANK)

# ------------------ RNA matrices (genes×cells) ------------------
rna_symcol = _infer_symbol_col(adata_rna)
rna_symbols = _normalize_symbols(adata_rna.var[rna_symcol]) if rna_symcol else _normalize_symbols(pd.Series(adata_rna.var_names, index=adata_rna.var_names))

def _matrix_to_genesxcells(adata, layer) -> pd.DataFrame:
    X = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    df = pd.DataFrame(X, index=adata.obs_names, columns=rna_symbols.values)  # cells×genes
    return df.T  # genes×cells

rna_counts_layer = "counts" if "counts" in adata_rna.layers else None
rna_log1p_layer  = _ensure_log1p_layer(adata_rna, prefer="counts", layer_name="log1p")

rna_counts_gxc = _matrix_to_genesxcells(adata_rna, rna_counts_layer) if rna_counts_layer else None
rna_log1p_gxc  = _matrix_to_genesxcells(adata_rna, rna_log1p_layer)

# RNA duplicates: SUM
if rna_counts_gxc is not None:
    rna_counts_gxc = _dedup_index(rna_counts_gxc, how='sum')
rna_log1p_gxc = _dedup_index(rna_log1p_gxc, how='sum')

# ------------------ COMMON3 gene set (same for all outputs) ------------------
rna_all_genes  = list(rna_log1p_gxc.index)
ac_all_genes   = list(ga_ac_counts.index)
me3_all_genes  = list(ga_me3_counts.index)

ac_set   = set(ac_all_genes)
me3_set  = set(me3_all_genes)

# Preserve RNA order
COMMON3_GENES = [g for g in rna_all_genes if (g in ac_set and g in me3_set)]
if len(COMMON3_GENES) == 0:
    raise RuntimeError("COMMON3 gene set is empty. Gene symbols may be mismatched across modalities.")
print(f"[info] COMMON3 genes: {len(COMMON3_GENES)}")

# --- Hard order anchors (genes & cells) + export for verification ---
GENE_ORDER_COMMON3 = list(COMMON3_GENES)       # preserve RNA-derived order
CELL_ORDER         = list(adata_rna.obs_names) # RNA cell order
CELL_LABELS        = _cell_labels_with_sample(adata_rna)

with open(os.path.join(OUTDIR, "COMMON3_genes_order.txt"), "w") as f:
    f.write("\n".join(GENE_ORDER_COMMON3) + "\n")
with open(os.path.join(OUTDIR, "COMMON3_cells_order.txt"), "w") as f:
    f.write("\n".join(CELL_LABELS) + "\n")
print(f"[info] wrote COMMON3_genes_order.txt ({len(GENE_ORDER_COMMON3)} genes) and COMMON3_cells_order.txt ({len(CELL_LABELS)} cells)")

# ------------------ Markers from RNA (intersect with COMMON3) ------------------
def top_markers_from_rgg(
    adata, rgg_key='CL_rgg', groupby='CL', top_n=10,
    alpha=0.05, min_logfc=0.25, ensure_unique=True, use_score=True
) -> Dict[str, List[str]]:
    if rgg_key not in adata.uns or groupby not in adata.obs:
        return {}
    rgg = adata.uns[rgg_key]
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    groups = [g for g in adata.obs[groupby].cat.categories if g in list(rgg['names'].dtype.names)]

    valid_symbols = set(rna_symbols.values)
    used = set(); out = {}
    for g in groups:
        names  = np.array(rgg['names'][g], dtype=str)
        pvals  = np.array(rgg['pvals_adj'][g], dtype=float)
        lfc    = np.array(rgg['logfoldchanges'][g], dtype=float)
        scores = np.array(rgg['scores'][g], dtype=float)
        names  = np.array(_normalize_symbols(pd.Series(names)), dtype=str)
        df = pd.DataFrame({'gene':names, 'p_adj':pvals, 'logFC':lfc, 'score':scores})
        df = df[df['gene'].isin(valid_symbols)]
        df = df[(df['p_adj'] < alpha) & (df['logFC'] > min_logfc) & (df['score'] > 0)]
        df = df.sort_values('score' if use_score else 'logFC', ascending=False)
        if ensure_unique and used:
            df = df[~df['gene'].isin(used)]
        picks = df['gene'].head(top_n).tolist()
        out[g] = picks
        used.update(picks)
    return out

markers_dict  = top_markers_from_rgg(
    adata_rna, rgg_key=RGG_KEY, groupby=RGG_GROUPBY,
    top_n=RGG_TOP_N, alpha=RGG_ALPHA, min_logfc=RGG_MIN_LFC,
    ensure_unique=RGG_UNIQUE, use_score=RGG_USE_SCORE
)
markers_union   = []
if markers_dict:
    for g in adata_rna.obs[RGG_GROUPBY].cat.categories:
        if g in markers_dict:
            markers_union.extend(markers_dict[g])
markers_union      = list(dict.fromkeys(markers_union))
MARKERS_COMMON3    = [g for g in markers_union if g in set(GENE_ORDER_COMMON3)]
print(f"[info] MARKERS_COMMON3 genes: {len(MARKERS_COMMON3)}")

# ------------------ Strategies ------------------
# RNA strategies: map name -> (DataFrame genes×cells, apply_fn or flag)
rna_strategies = {}
if rna_counts_gxc is not None:
    rna_strategies["expr_umi1"] = (rna_counts_gxc,  lambda M: (M >= 1).astype(np.uint8))
    rna_strategies["expr_umi2"] = (rna_counts_gxc,  lambda M: (M >= 2).astype(np.uint8))
rna_strategies["expr_log1p_gt0"]   = (rna_log1p_gxc, lambda M: (M > 0).astype(np.uint8))
rna_strategies["expr_log1p_z_gt1"] = (rna_log1p_gxc, "PER_GENE_Z>1")  # special flag

# GA (AC/ME3) base matrices restricted to COMMON3
ac_counts = ga_ac_counts.loc[GENE_ORDER_COMMON3]
ac_log1p  = ga_ac_log1p.loc[GENE_ORDER_COMMON3]
ac_cpm    = _cpm_df(ac_counts, GA_CPM_SCALE)

me3_counts = ga_me3_counts.loc[GENE_ORDER_COMMON3]
me3_log1p  = ga_me3_log1p.loc[GENE_ORDER_COMMON3]
me3_cpm    = _cpm_df(me3_counts, GA_CPM_SCALE)

ga_strategies_ac = {
    "ga_any1":        (ac_counts, lambda M: (M >= 1).astype(np.uint8)),
    "ga_min2":        (ac_counts, lambda M: (M >= 2).astype(np.uint8)),
    "ga_cpm_ge1":     (ac_cpm,    lambda M: (M >= 1.0).astype(np.uint8)),
    "ga_log1p_z_gt1": (ac_log1p,  "PER_GENE_Z>1"),
}
ga_strategies_me3 = {
    "ga_any1":        (me3_counts, lambda M: (M >= 1).astype(np.uint8)),
    "ga_min2":        (me3_counts, lambda M: (M >= 2).astype(np.uint8)),
    "ga_cpm_ge1":     (me3_cpm,    lambda M: (M >= 1.0).astype(np.uint8)),
    "ga_log1p_z_gt1": (me3_log1p,  "PER_GENE_Z>1"),
}

# ------------------ Writer ------------------
manifest_rows = []

def _process_and_write(modality: str, strategy: str, df_raw_gxc: pd.DataFrame,
                       rule, subset_genes: List[str], tag: str):
    """
    df_raw_gxc: genes×cells (floats)
    rule: lambda M->binary or "PER_GENE_Z>1"
    subset_genes: desired output gene list (we'll intersect with df_raw_gxc.index)
    """
    # Intersect while preserving requested order
    present = set(df_raw_gxc.index)
    genes   = [g for g in subset_genes if g in present]
    if not genes:
        print(f"[skip] {modality}/{strategy}/{tag}: 0 genes after intersection.")
        return

    # --- Enforce canonical orders (rows & columns) ---
    if tag == "COMMON3":
        df = df_raw_gxc.reindex(GENE_ORDER_COMMON3).fillna(0.0).astype(float)
        df = df.reindex(columns=CELL_ORDER)
    else:
        df = df_raw_gxc.reindex(genes).fillna(0.0).astype(float)
        df = df.reindex(columns=CELL_ORDER)

    # Binarize
    if isinstance(rule, str) and rule == "PER_GENE_Z>1":
        z = _z_per_gene_matrix(df)
        bin_df = (z > 1.0).astype(np.uint8)
    else:
        bin_df = rule(df)

    # Rename columns to the *same* sample|barcode labels everywhere
    bin_df.columns = CELL_LABELS

    # Write
    fname = f"binary_{modality}_{strategy}__{tag}.tsv.gz"
    path  = os.path.join(OUTDIR, fname)
    write_tsv_gz(bin_df, path)

    notes = {
        "expr_umi1": "RNA counts >= 1",
        "expr_umi2": "RNA counts >= 2",
        "expr_log1p_gt0": "RNA log1p > 0",
        "expr_log1p_z_gt1": "RNA log1p per-gene z-score > 1",
        "ga_any1": "GA raw count >= 1",
        "ga_min2": "GA raw count >= 2",
        "ga_cpm_ge1": f"GA CPM per {int(GA_CPM_SCALE):,} >= 1",
        "ga_log1p_z_gt1": "GA log1p per-gene z-score > 1",
    }.get(strategy, "")

    manifest_rows.append({
        "file": fname, "modality": modality, "strategy": strategy,
        "genes": tag, "notes": notes
    })
    print(f"[ok] {fname}  ({bin_df.shape[0]} genes × {bin_df.shape[1]} cells)")

# ------------------ RUN (COMMON3 + MARKERS_COMMON3) ------------------
# Subset RNA matrices to COMMON3 once, with enforced row order
rna_counts_common = rna_counts_gxc.loc[GENE_ORDER_COMMON3] if rna_counts_gxc is not None else None
rna_log1p_common  = rna_log1p_gxc.loc[GENE_ORDER_COMMON3]

rna_strategies_common = {}
if rna_counts_common is not None:
    rna_strategies_common["expr_umi1"] = (rna_counts_common,  rna_strategies["expr_umi1"][1])
    rna_strategies_common["expr_umi2"] = (rna_counts_common,  rna_strategies["expr_umi2"][1])
rna_strategies_common["expr_log1p_gt0"]   = (rna_log1p_common, lambda M: (M > 0).astype(np.uint8))
rna_strategies_common["expr_log1p_z_gt1"] = (rna_log1p_common, "PER_GENE_Z>1")

# RNA
for strat, (df_raw, rule) in rna_strategies_common.items():
    _process_and_write("Expr", strat, df_raw, rule, GENE_ORDER_COMMON3, "COMMON3")
    if len(MARKERS_COMMON3):
        _process_and_write("Expr", strat, df_raw, rule, MARKERS_COMMON3, "MARKERS_COMMON3")

# AC
for strat, (df_raw, rule) in ga_strategies_ac.items():
    _process_and_write(AC_MODALITY, strat, df_raw, rule, GENE_ORDER_COMMON3, "COMMON3")
    if len(MARKERS_COMMON3):
        _process_and_write(AC_MODALITY, strat, df_raw, rule, MARKERS_COMMON3, "MARKERS_COMMON3")

# ME3
for strat, (df_raw, rule) in ga_strategies_me3.items():
    _process_and_write(ME3_MODALITY, strat, df_raw, rule, GENE_ORDER_COMMON3, "COMMON3")
    if len(MARKERS_COMMON3):
        _process_and_write(ME3_MODALITY, strat, df_raw, rule, MARKERS_COMMON3, "MARKERS_COMMON3")

# ------------------ Manifest ------------------
manifest = pd.DataFrame(manifest_rows, columns=["file","modality","strategy","genes","notes"])
manifest_path = os.path.join(OUTDIR, "MANIFEST_binary_matrices.tsv")
manifest.to_csv(manifest_path, sep='\t', index=False)
print(f"[done] Manifest: {manifest_path}")
