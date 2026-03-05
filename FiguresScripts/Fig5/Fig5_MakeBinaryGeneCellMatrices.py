#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_binary_gene_cell_matrices_COMMON3.py  (UPDATED: ALL + byClusterCond in ONE folder)

Build 0/1 matrices (genes x cells) for:
  - RNA expression (Expr)
  - H3K27ac gene-activity (±1 kb around gene body)
  - H3K27me3 gene-activity (±1 kb around gene body)

This version RESTRICTS to genes present in ALL THREE modalities (COMMON3),
so the three outputs (RNA / AC / ME3) have identical gene rows and ordering.

OUTPUT LAYOUT
  OUTROOT/
    ALL/                      # all cells together, natural order
      binary_*.tsv.gz
      COMMON3_genes_order.txt
      COMMON3_cells_order.txt
      MANIFEST_binary_matrices.tsv
    byClusterCond/            # all cells together, but columns ordered by cluster_call×condition
      binary_*.tsv.gz
      COMMON3_genes_order.txt
      COMMON3_cells_order.txt
      group_sizes.tsv
      cell_to_group.tsv
      MANIFEST_binary_matrices.tsv

Key behavior change vs the previous adaptation:
  - byClusterCond is NOT split into per-group subfolders anymore.
  - Instead, byClusterCond writes ONE set of matrices (like “cell lines” in the original script),
    with columns labeled as "<cluster_call>_<condition>|<barcode>" and ordered by group.

No external GTF: RNA gene list from adata_rna.var; GA genes from SnapATAC2 using snap.genome.hg38.

Requirements: scanpy, anndata, snapatac2, pandas, numpy, scipy
"""

import os, gzip, re
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import scipy.sparse as sp
from typing import Dict, List, Tuple, Optional

# ------------------ CONFIG ------------------
RNA_PATH = "RNA_Cd47i_Untreated_FiltShared.h5ad"
AC_PATH  = "H3K27ac_Cd47i_Untreated_FiltShared.h5ad"
ME3_PATH = "H3K27me3_Cd47i_Untreated_FiltShared.h5ad"

OUTROOT = "BinaryMatrices_OUT"   # will create ALL/ and byClusterCond/
os.makedirs(OUTROOT, exist_ok=True)

# Split keys (must exist in at least one modality's .obs)
CONDITION_COL = "condition"
CLUSTER_COL   = "cluster_call"

# rank_genes_groups keys to pull markers (optional; if not present, markers are skipped)
RGG_KEY       = "CL_rgg"
RGG_GROUPBY   = "CL"
RGG_TOP_N     = 50
RGG_ALPHA     = 0.05
RGG_MIN_LFC   = 0.25
RGG_UNIQUE    = True
RGG_USE_SCORE = True

# Gene-activity window (±1 kb around gene body)
GA_FLANK = 5_000

# For output naming
GA_FLANK_TAG = f"pm{GA_FLANK // 1000}kb" if (GA_FLANK % 1000 == 0) else f"pm{GA_FLANK}bp"
AC_MODALITY  = f"H3K27ac_{GA_FLANK_TAG}"
ME3_MODALITY = f"H3K27me3_{GA_FLANK_TAG}"

# CPM scaling for GA (counts per 1e5)
GA_CPM_SCALE = 1e5

# ------------------ HELPERS ------------------
def _sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))

def _normalize_symbols(s: pd.Series) -> pd.Series:
    """Turn 'A,B' -> 'A', strip spaces."""
    return s.astype(str).str.split(',').str[0].str.strip()

def _infer_symbol_col(adata) -> Optional[str]:
    for c in ['gene_symbol','gene_symbols','symbol','SYMBOL','feature_name','GeneSymbol']:
        if c in adata.var.columns:
            return c
    return None

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

def _maybe_strip_suffix_after_last_underscore(adata, ref_index: pd.Index, col_name="barcode_orig"):
    """
    If stripping after last '_' increases overlap with ref_index, apply it; else keep as-is.
    """
    adata = adata.copy()
    idx0 = pd.Index(adata.obs_names)
    overlap0 = idx0.intersection(ref_index).size

    if not idx0.str.contains("_").any():
        return adata

    new_idx = idx0.str.rsplit("_", n=1).str[0]
    overlap1 = pd.Index(new_idx).intersection(ref_index).size

    if overlap1 > overlap0:
        adata.obs[col_name] = idx0.astype(str)
        adata.obs_names = pd.Index(new_idx)
        # guard against duplicates
        if not adata.obs_names.is_unique:
            raise ValueError("[ERR] Stripping suffix created duplicate obs_names; aborting.")
        print(f"[info] stripped suffix after last '_' for better barcode match: {overlap0} -> {overlap1}")
    return adata

def _pull_obs_col(colname: str, adatas: Tuple[sc.AnnData, sc.AnnData, sc.AnnData]) -> pd.Series:
    """Find obs col in (RNA, AC, ME3) in that order; return aligned to RNA obs_names."""
    a1, a2, a3 = adatas
    for ad in (a1, a2, a3):
        if colname in ad.obs.columns:
            return ad.obs[colname].copy().reindex(a1.obs_names)
    raise KeyError(f"Required column '{colname}' not found in .obs of RNA/AC/ME3.")

def _cell_labels_all(adata_rna: sc.AnnData) -> List[str]:
    """ALL folder labels: sample|barcode if available, else barcode."""
    if "sample" in adata_rna.obs.columns:
        return [f"{str(adata_rna.obs.loc[c,'sample'])}|{c}" for c in adata_rna.obs_names]
    return list(adata_rna.obs_names)

# ------------------ LOAD ------------------
adata_rna = sc.read_h5ad(RNA_PATH)
adata_ac  = sc.read_h5ad(AC_PATH)
adata_me3 = sc.read_h5ad(ME3_PATH)

# Improve barcode matching only if needed
adata_ac  = _maybe_strip_suffix_after_last_underscore(adata_ac,  adata_rna.obs_names)
adata_me3 = _maybe_strip_suffix_after_last_underscore(adata_me3, adata_rna.obs_names)

# Keep shared cells; preserve RNA order
common_cells = adata_rna.obs_names.intersection(adata_ac.obs_names).intersection(adata_me3.obs_names)
if common_cells.size == 0:
    raise RuntimeError("[ERR] No overlapping cells across RNA/AC/ME3 after optional stripping.")
adata_rna = adata_rna[common_cells].copy()
adata_ac  = adata_ac[common_cells].copy()
adata_me3 = adata_me3[common_cells].copy()

# Ensure split metadata present in all modalities (copy into each .obs)
condition = _pull_obs_col(CONDITION_COL, (adata_rna, adata_ac, adata_me3))
cluster   = _pull_obs_col(CLUSTER_COL,   (adata_rna, adata_ac, adata_me3))
for ad in (adata_rna, adata_ac, adata_me3):
    ad.obs[CONDITION_COL] = condition.values
    ad.obs[CLUSTER_COL]   = cluster.values

# One canonical group label per cell
group_cc_cond = (cluster.astype(str) + "_" + condition.astype(str)).map(_sanitize)
for ad in (adata_rna, adata_ac, adata_me3):
    ad.obs["group_cc_cond"] = group_cc_cond.values

print("[info] Cells per cluster_call×condition:")
print(pd.Series(group_cc_cond.values).value_counts().sort_index())

# ------------------ MAKE GENE ACTIVITY (±1 kb, NO GTF) ------------------
try:
    GENE_ANNO = snap.genome.hg38
except Exception as e:
    raise RuntimeError("snap.genome.hg38 is required (no external GTF used).") from e

def make_gene_activity_tables(
    adata_mark: sc.AnnData,
    gene_anno,
    flank_bp: int,
    cell_order: List[str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    SnapATAC2 GA (cell×gene) → transpose to genes×cells; also build log1p DF.
    """
    ga = snap.pp.make_gene_matrix(
        adata_mark,
        gene_anno,
        include_gene_body=True,
        upstream=flank_bp,
        downstream=flank_bp,
        id_type="gene",
    )
    ga = ga[cell_order, :]  # align cells to RNA

    X = ga.X.toarray() if sp.issparse(ga.X) else np.asarray(ga.X)
    counts_gxc = X.T  # genes×cells

    Xlog = ga.layers.get("log1p", None)
    if Xlog is None:
        if sp.issparse(ga.X):
            Xlog = ga.X.copy(); Xlog.data = np.log1p(Xlog.data)
        else:
            Xlog = np.log1p(ga.X)
    Xlog = Xlog.toarray() if sp.issparse(Xlog) else np.asarray(Xlog)
    log1p_gxc = Xlog.T  # genes×cells

    gene_idx = _normalize_symbols(pd.Series(ga.var_names, index=ga.var_names)).values
    cell_cols = list(ga.obs_names)

    df_counts = pd.DataFrame(counts_gxc, index=gene_idx, columns=cell_cols)
    df_log1p  = pd.DataFrame(log1p_gxc,  index=gene_idx, columns=cell_cols)

    # GA duplicates: presence-like => MAX
    df_counts = _dedup_index(df_counts, how='max')
    df_log1p  = _dedup_index(df_log1p,  how='max')
    return df_counts, df_log1p

# ------------------ RNA matrices (genes×cells) ------------------
rna_symcol = _infer_symbol_col(adata_rna)
rna_symbols = (
    _normalize_symbols(adata_rna.var[rna_symcol])
    if rna_symcol else
    _normalize_symbols(pd.Series(adata_rna.var_names, index=adata_rna.var_names))
)

def _matrix_to_genesxcells(adata: sc.AnnData, layer: Optional[str], gene_symbols: pd.Series) -> pd.DataFrame:
    X = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    df = pd.DataFrame(X, index=adata.obs_names, columns=gene_symbols.values)  # cells×genes
    return df.T  # genes×cells

rna_counts_layer = "counts" if "counts" in adata_rna.layers else None
rna_log1p_layer  = _ensure_log1p_layer(adata_rna, prefer="counts", layer_name="log1p")

rna_counts_gxc = _matrix_to_genesxcells(adata_rna, rna_counts_layer, rna_symbols) if rna_counts_layer else None
rna_log1p_gxc  = _matrix_to_genesxcells(adata_rna, rna_log1p_layer,  rna_symbols)

# RNA duplicates: SUM
if rna_counts_gxc is not None:
    rna_counts_gxc = _dedup_index(rna_counts_gxc, how='sum')
rna_log1p_gxc = _dedup_index(rna_log1p_gxc, how='sum')

# ------------------ GA matrices (genes×cells) for ALL cells ------------------
CELL_ORDER_ALL = list(adata_rna.obs_names)

ga_ac_counts_all,  ga_ac_log1p_all  = make_gene_activity_tables(adata_ac,  GENE_ANNO, GA_FLANK, CELL_ORDER_ALL)
ga_me3_counts_all, ga_me3_log1p_all = make_gene_activity_tables(adata_me3, GENE_ANNO, GA_FLANK, CELL_ORDER_ALL)

# ------------------ COMMON3 gene set (ALL) ------------------
rna_all_genes = list(rna_log1p_gxc.index)
ac_set  = set(ga_ac_counts_all.index)
me3_set = set(ga_me3_counts_all.index)

COMMON3_GENES_ALL = [g for g in rna_all_genes if (g in ac_set and g in me3_set)]
if len(COMMON3_GENES_ALL) == 0:
    raise RuntimeError("COMMON3 gene set is empty. Gene symbols may be mismatched across modalities.")
print(f"[info] COMMON3 genes (ALL): {len(COMMON3_GENES_ALL)}")

GENE_ORDER_ALL = list(COMMON3_GENES_ALL)

# ------------------ Markers from RNA (optional; intersect with COMMON3) ------------------
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

markers_dict = top_markers_from_rgg(
    adata_rna, rgg_key=RGG_KEY, groupby=RGG_GROUPBY,
    top_n=RGG_TOP_N, alpha=RGG_ALPHA, min_logfc=RGG_MIN_LFC,
    ensure_unique=RGG_UNIQUE, use_score=RGG_USE_SCORE
)

markers_union = []
if markers_dict:
    for g in adata_rna.obs[RGG_GROUPBY].cat.categories:
        if g in markers_dict:
            markers_union.extend(markers_dict[g])
markers_union = list(dict.fromkeys(markers_union))
MARKERS_COMMON3_ALL = [g for g in markers_union if g in set(GENE_ORDER_ALL)]
print(f"[info] MARKERS_COMMON3 genes (ALL): {len(MARKERS_COMMON3_ALL)}")

# ------------------ STRATEGIES (defined once) ------------------
def build_strategies(
    rna_counts: Optional[pd.DataFrame],
    rna_log1p: pd.DataFrame,
    ac_counts: pd.DataFrame,
    ac_log1p: pd.DataFrame,
    me3_counts: pd.DataFrame,
    me3_log1p: pd.DataFrame,
) -> Tuple[Dict, Dict, Dict]:
    rna_strats = {}
    if rna_counts is not None:
        rna_strats["expr_umi1"] = (rna_counts,  lambda M: (M >= 1).astype(np.uint8))
        rna_strats["expr_umi2"] = (rna_counts,  lambda M: (M >= 2).astype(np.uint8))
    rna_strats["expr_log1p_gt0"]   = (rna_log1p, lambda M: (M > 0).astype(np.uint8))
    rna_strats["expr_log1p_z_gt1"] = (rna_log1p, "PER_GENE_Z>1")

    ac_cpm  = _cpm_df(ac_counts, GA_CPM_SCALE)
    me3_cpm = _cpm_df(me3_counts, GA_CPM_SCALE)

    ac_strats = {
        "ga_any1":        (ac_counts, lambda M: (M >= 1).astype(np.uint8)),
        "ga_min2":        (ac_counts, lambda M: (M >= 2).astype(np.uint8)),
        "ga_cpm_ge1":     (ac_cpm,    lambda M: (M >= 1.0).astype(np.uint8)),
        "ga_log1p_z_gt1": (ac_log1p,  "PER_GENE_Z>1"),
    }
    me3_strats = {
        "ga_any1":        (me3_counts, lambda M: (M >= 1).astype(np.uint8)),
        "ga_min2":        (me3_counts, lambda M: (M >= 2).astype(np.uint8)),
        "ga_cpm_ge1":     (me3_cpm,    lambda M: (M >= 1.0).astype(np.uint8)),
        "ga_log1p_z_gt1": (me3_log1p,  "PER_GENE_Z>1"),
    }
    return rna_strats, ac_strats, me3_strats

# Restrict “base” matrices to ALL COMMON3 gene order once
rna_counts_all = rna_counts_gxc.loc[GENE_ORDER_ALL] if rna_counts_gxc is not None else None
rna_log1p_all  = rna_log1p_gxc.loc[GENE_ORDER_ALL]

ac_counts_all  = ga_ac_counts_all.loc[GENE_ORDER_ALL]
ac_log1p_all   = ga_ac_log1p_all.loc[GENE_ORDER_ALL]
me3_counts_all = ga_me3_counts_all.loc[GENE_ORDER_ALL]
me3_log1p_all  = ga_me3_log1p_all.loc[GENE_ORDER_ALL]

RNA_STRATS_ALL, AC_STRATS_ALL, ME3_STRATS_ALL = build_strategies(
    rna_counts_all, rna_log1p_all,
    ac_counts_all, ac_log1p_all,
    me3_counts_all, me3_log1p_all
)

# ------------------ EXPORT CORE ------------------
def export_folder(
    outdir: str,
    cell_names: List[str],
    cell_labels: List[str],
    gene_order: List[str],
    markers_common3: List[str],
    rna_strats: Dict,
    ac_strats: Dict,
    me3_strats: Dict,
):
    os.makedirs(outdir, exist_ok=True)
    manifest_rows = []

    # write order anchors
    with open(os.path.join(outdir, "COMMON3_genes_order.txt"), "w") as f:
        f.write("\n".join(gene_order) + "\n")
    with open(os.path.join(outdir, "COMMON3_cells_order.txt"), "w") as f:
        f.write("\n".join(cell_labels) + "\n")

    def _process_and_write(modality: str, strategy: str, df_raw_gxc: pd.DataFrame,
                           rule, subset_genes: List[str], tag: str):
        present = set(df_raw_gxc.index)
        genes = [g for g in subset_genes if g in present]
        if not genes:
            print(f"[skip] {outdir}: {modality}/{strategy}/{tag}: 0 genes after intersection.")
            return

        # subset to requested genes and requested cells (preserve order)
        if tag == "COMMON3":
            df = df_raw_gxc.reindex(gene_order).fillna(0.0).astype(float)
        else:
            df = df_raw_gxc.reindex(genes).fillna(0.0).astype(float)

        df = df.reindex(columns=cell_names)

        # binarize
        if isinstance(rule, str) and rule == "PER_GENE_Z>1":
            z = _z_per_gene_matrix(df)
            bin_df = (z > 1.0).astype(np.uint8)
        else:
            bin_df = rule(df)

        # rename columns to stable labels
        bin_df.columns = cell_labels

        fname = f"binary_{modality}_{strategy}__{tag}.tsv.gz"
        path  = os.path.join(outdir, fname)
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
        print(f"[ok] {os.path.relpath(path)}  ({bin_df.shape[0]} genes × {bin_df.shape[1]} cells)")

    # RNA
    for strat, (df_raw, rule) in rna_strats.items():
        _process_and_write("Expr", strat, df_raw, rule, gene_order, "COMMON3")
        if len(markers_common3):
            _process_and_write("Expr", strat, df_raw, rule, markers_common3, "MARKERS_COMMON3")

    # AC
    for strat, (df_raw, rule) in ac_strats.items():
        _process_and_write(AC_MODALITY, strat, df_raw, rule, gene_order, "COMMON3")
        if len(markers_common3):
            _process_and_write(AC_MODALITY, strat, df_raw, rule, markers_common3, "MARKERS_COMMON3")

    # ME3
    for strat, (df_raw, rule) in me3_strats.items():
        _process_and_write(ME3_MODALITY, strat, df_raw, rule, gene_order, "COMMON3")
        if len(markers_common3):
            _process_and_write(ME3_MODALITY, strat, df_raw, rule, markers_common3, "MARKERS_COMMON3")

    manifest = pd.DataFrame(manifest_rows, columns=["file","modality","strategy","genes","notes"])
    manifest_path = os.path.join(outdir, "MANIFEST_binary_matrices.tsv")
    manifest.to_csv(manifest_path, sep='\t', index=False)
    print(f"[done] Manifest: {manifest_path}")

# ------------------ EXPORT: ALL (single folder) ------------------
out_all = os.path.join(OUTROOT, "ALL")
cell_names_all  = CELL_ORDER_ALL
cell_labels_all = _cell_labels_all(adata_rna)

export_folder(
    outdir=out_all,
    cell_names=cell_names_all,
    cell_labels=cell_labels_all,
    gene_order=GENE_ORDER_ALL,
    markers_common3=MARKERS_COMMON3_ALL,
    rna_strats=RNA_STRATS_ALL,
    ac_strats=AC_STRATS_ALL,
    me3_strats=ME3_STRATS_ALL,
)

# ------------------ EXPORT: byClusterCond (single folder; grouped columns) ------------------
out_bcc = os.path.join(OUTROOT, "byClusterCond")
os.makedirs(out_bcc, exist_ok=True)

groups = pd.Series(adata_rna.obs["group_cc_cond"].values, index=adata_rna.obs_names)

# define grouped cell order: group-sorted, within-group preserve RNA order
group_levels = sorted(groups.unique().tolist())
cell_names_bcc: List[str] = []
cell_labels_bcc: List[str] = []

for g in group_levels:
    cells_g = [c for c in CELL_ORDER_ALL if groups.loc[c] == g]
    cell_names_bcc.extend(cells_g)
    cell_labels_bcc.extend([f"{g}|{c}" for c in cells_g])

# write convenience mapping tables
group_sizes = pd.Series(groups.values).value_counts().reindex(group_levels).fillna(0).astype(int)
pd.DataFrame({"group": group_sizes.index, "n_cells": group_sizes.values}).to_csv(
    os.path.join(out_bcc, "group_sizes.tsv"), sep="\t", index=False
)
pd.DataFrame({
    "barcode": CELL_ORDER_ALL,
    "group_cc_cond": groups.reindex(CELL_ORDER_ALL).values,
    "cluster_call": adata_rna.obs[CLUSTER_COL].reindex(CELL_ORDER_ALL).astype(str).values,
    "condition": adata_rna.obs[CONDITION_COL].reindex(CELL_ORDER_ALL).astype(str).values,
}).to_csv(os.path.join(out_bcc, "cell_to_group.tsv"), sep="\t", index=False)

export_folder(
    outdir=out_bcc,
    cell_names=cell_names_bcc,
    cell_labels=cell_labels_bcc,
    gene_order=GENE_ORDER_ALL,
    markers_common3=MARKERS_COMMON3_ALL,
    rna_strats=RNA_STRATS_ALL,
    ac_strats=AC_STRATS_ALL,
    me3_strats=ME3_STRATS_ALL,
)

print(f"[DONE] Wrote outputs to: {OUTROOT}")
