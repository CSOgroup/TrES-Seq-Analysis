#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gzip
import shutil
import tempfile
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import scipy.sparse as sp
from scipy.io import mmwrite

# --------------------------------------------------
# paths / settings
# --------------------------------------------------
RNA_PATH = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"
AC_PATH  = "DNA_merged/MergedTonsils_H3K27ac_merged_processed_barcodeRewritten_withRNAlabels.h5ad"
ME3_PATH = "DNA_merged/MergedTonsils_H3K27me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad"

OUTDIR = "overlap3mods_allgenes_exports"
PREFIX = "tonsil_overlap3mods_Bcells_with_plasma"
os.makedirs(OUTDIR, exist_ok=True)

GROUP_KEY = "leiden_merged_type"
B_PREFIX = "B_"   # includes B_plasma if label starts with B_

# RNA layer used for expression export
# change to "SCT_counts" if that is what you want
RNA_COUNTS_LAYER = "SCT_counts"
RNA_CPM_TARGET_SUM = 1e6

# SnapATAC2 gene activity: gene body + ±1 kb = "pm2kb_genebody" in filenames
UPSTREAM = 2000
DOWNSTREAM = 2000
INCLUDE_GENE_BODY = True
GENE_ANNO = snap.genome.GRCh38

REUSE_GA_IF_EXISTS = True


# --------------------------------------------------
# helpers
# --------------------------------------------------
def _to_dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def _normalize_symbols(symbols):
    return pd.Series(symbols, dtype=object).astype(str).str.split(",").str[0].str.strip()


def _infer_rna_gene_symbols(adata):
    for col in ["gene_symbol", "gene_symbols", "symbol", "SYMBOL", "feature_name", "GeneSymbol"]:
        if col in adata.var.columns:
            return _normalize_symbols(adata.var[col].values)
    return _normalize_symbols(adata.var_names)


def _maybe_strip_suffix_after_last_underscore(adata, ref_index, col_name="barcode_orig"):
    """
    If obs names like BARCODE_suffix match RNA better after stripping the final suffix, do that.
    """
    adata = adata.copy()
    idx0 = pd.Index(adata.obs_names.astype(str))
    overlap0 = idx0.intersection(ref_index).size

    if not idx0.str.contains("_").any():
        return adata

    idx1 = idx0.str.rsplit("_", n=1).str[0]
    overlap1 = pd.Index(idx1).intersection(ref_index).size

    if overlap1 > overlap0:
        adata.obs[col_name] = idx0.values
        adata.obs_names = pd.Index(idx1)
        if not adata.obs_names.is_unique:
            raise ValueError("Stripping final suffix created duplicate obs_names.")
        print(f"[info] stripped final suffix for better overlap: {overlap0} -> {overlap1}")

    return adata


def infer_use_x(adata_mod):
    return ("insertion" not in adata_mod.obsm)


def compute_rna_log1p_cpm_from_counts_layer(rna_adata, counts_layer="counts", target_sum=1e6):
    if counts_layer not in rna_adata.layers:
        raise ValueError(
            f"RNA layer '{counts_layer}' not found. Available layers: {list(rna_adata.layers.keys())}"
        )
    tmp = rna_adata.copy()
    tmp.X = tmp.layers[counts_layer].copy()
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    sc.pp.log1p(tmp)
    return tmp


def make_gene_activity_log1p_cpm(adata_mod, label, out_h5ad):
    if REUSE_GA_IF_EXISTS and os.path.exists(out_h5ad):
        ga = sc.read_h5ad(out_h5ad)
        print(f"[{label}] loaded cached gene activity: {out_h5ad}")
        return ga

    use_x = infer_use_x(adata_mod)
    print(
        f"[{label}] computing gene activity "
        f"(use_x={use_x}, gene_body={INCLUDE_GENE_BODY}, up={UPSTREAM}, down={DOWNSTREAM})"
    )

    ga = snap.pp.make_gene_matrix(
        adata_mod,
        gene_anno=GENE_ANNO,
        inplace=False,
        file=None,
        use_x=use_x,
        upstream=UPSTREAM,
        downstream=DOWNSTREAM,
        include_gene_body=INCLUDE_GENE_BODY,
    )

    sc.pp.normalize_total(ga, target_sum=1e6)
    sc.pp.log1p(ga)

    ga.write(out_h5ad)
    print(f"[{label}] saved gene activity: {out_h5ad}")
    return ga


def genes_by_cells_df(adata, gene_symbols=None):
    """
    Return DataFrame with genes as rows and cells as columns.
    """
    X = _to_dense(adata.X)  # cells x genes
    if gene_symbols is None:
        genes = _normalize_symbols(adata.var_names).values
    else:
        genes = _normalize_symbols(gene_symbols).values

    df = pd.DataFrame(X.T, index=genes, columns=adata.obs_names)

    # collapse duplicates by max for GA-like matrices; sum would also be reasonable for RNA counts,
    # but for already normalized/logged values max is safer than summing
    if not df.index.is_unique:
        df = df.groupby(level=0, sort=False).max()

    return df


def write_mtx_gz(df_gxc, out_prefix):
    """
    df_gxc: genes x cells dataframe
    Writes:
      - genes_by_cells.mtx.gz
      - genes.txt
      - cells.txt
    """
    genes_txt = out_prefix + ".genes.txt"
    cells_txt = out_prefix + ".cells.txt"
    mtx_gz = out_prefix + ".genes_by_cells.mtx.gz"

    pd.Series(df_gxc.index).to_csv(genes_txt, index=False, header=False)
    pd.Series(df_gxc.columns).to_csv(cells_txt, index=False, header=False)

    mat = sp.csr_matrix(df_gxc.values)

    with tempfile.NamedTemporaryFile(suffix=".mtx", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        mmwrite(tmp_path, mat)
        with open(tmp_path, "rb") as f_in, gzip.open(mtx_gz, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def write_cells_by_genes_csv_gz(df_gxc, out_prefix):
    """
    Save cells x genes CSV.gz
    """
    df_cxg = df_gxc.T
    out_csv = out_prefix + ".cells_by_genes.csv.gz"
    df_cxg.to_csv(out_csv, compression="gzip")


# --------------------------------------------------
# load
# --------------------------------------------------
print("[load] reading objects...")
rna = sc.read_h5ad(RNA_PATH)
ac = sc.read_h5ad(AC_PATH)
me3 = sc.read_h5ad(ME3_PATH)

# improve matching if AC/ME3 have suffixed barcodes
ac = _maybe_strip_suffix_after_last_underscore(ac, rna.obs_names)
me3 = _maybe_strip_suffix_after_last_underscore(me3, rna.obs_names)

if GROUP_KEY not in rna.obs.columns:
    raise ValueError(f"{GROUP_KEY} not found in RNA obs")

# --------------------------------------------------
# define B-lineage cells from RNA
# --------------------------------------------------
bmask = rna.obs[GROUP_KEY].astype(str).str.startswith(B_PREFIX)
rna_b = rna[bmask].copy()
print(f"[RNA] B-lineage cells including plasma: {rna_b.n_obs}")

rna_b_cells = list(rna_b.obs_names)
ac_set = set(ac.obs_names)
me3_set = set(me3.obs_names)

common_cells = [c for c in rna_b_cells if (c in ac_set) and (c in me3_set)]
print(f"[overlap3] RNA(B) ∩ AC ∩ ME3: {len(common_cells)}")

if len(common_cells) == 0:
    raise ValueError("No overlap3 B cells found.")

with open(os.path.join(OUTDIR, f"{PREFIX}.overlap_cells.txt"), "w") as f:
    for c in common_cells:
        f.write(str(c) + "\n")

# metadata table
celltypes_df = rna_b.obs.loc[common_cells, [GROUP_KEY]].copy()
if "sample" in rna_b.obs.columns:
    celltypes_df["sample"] = rna_b.obs.loc[common_cells, "sample"].values
if "tonsil_id" in rna_b.obs.columns:
    celltypes_df["tonsil_id"] = rna_b.obs.loc[common_cells, "tonsil_id"].values
celltypes_df.index.name = "cell"
celltypes_df.to_csv(os.path.join(OUTDIR, f"{PREFIX}.overlap3mods.CELLTYPES.csv"))

# subset overlap3 cells
rna_ol = rna_b[common_cells].copy()
ac_ol = ac[common_cells].copy()
me3_ol = me3[common_cells].copy()

# --------------------------------------------------
# RNA log1p(CPM)
# --------------------------------------------------
rna_ol_logcpm = compute_rna_log1p_cpm_from_counts_layer(
    rna_ol,
    counts_layer=RNA_COUNTS_LAYER,
    target_sum=RNA_CPM_TARGET_SUM
)

rna_gene_symbols = _infer_rna_gene_symbols(rna_ol_logcpm)
rna_gxc = genes_by_cells_df(rna_ol_logcpm, gene_symbols=rna_gene_symbols)

# --------------------------------------------------
# gene activity for overlap3 AC / ME3
# --------------------------------------------------
ac_ga_path = os.path.join(
    OUTDIR,
    f"{PREFIX}.H3K27ac.gene_activity_pm2kb_genebody.overlap3mods.h5ad"
)
me3_ga_path = os.path.join(
    OUTDIR,
    f"{PREFIX}.H3K27me3.gene_activity_pm2kb_genebody.overlap3mods.h5ad"
)

ac_ga = make_gene_activity_log1p_cpm(ac_ol, "H3K27ac", ac_ga_path)
me3_ga = make_gene_activity_log1p_cpm(me3_ol, "H3K27me3", me3_ga_path)

ac_gxc = genes_by_cells_df(ac_ga)
me3_gxc = genes_by_cells_df(me3_ga)

# --------------------------------------------------
# restrict to genes present in all 3
# --------------------------------------------------
common_genes = [g for g in rna_gxc.index if (g in ac_gxc.index) and (g in me3_gxc.index)]
print(f"[genes] common across RNA / AC / ME3: {len(common_genes)}")

if len(common_genes) == 0:
    raise ValueError("No common genes across the three modalities.")

rna_gxc = rna_gxc.loc[common_genes, common_cells]
ac_gxc = ac_gxc.loc[common_genes, common_cells]
me3_gxc = me3_gxc.loc[common_genes, common_cells]

# --------------------------------------------------
# export
# --------------------------------------------------
rna_prefix = os.path.join(
    OUTDIR,
    f"{PREFIX}.RNA.log1pCPM_from_counts.overlap3mods"
)
ac_prefix = os.path.join(
    OUTDIR,
    f"{PREFIX}.H3K27ac.gene_activity.log1pCPM.overlap3mods"
)
me3_prefix = os.path.join(
    OUTDIR,
    f"{PREFIX}.H3K27me3.gene_activity.log1pCPM.overlap3mods"
)

# RNA
write_cells_by_genes_csv_gz(rna_gxc, rna_prefix)
write_mtx_gz(rna_gxc, rna_prefix)

# AC
write_cells_by_genes_csv_gz(ac_gxc, ac_prefix)
write_mtx_gz(ac_gxc, ac_prefix)

# ME3
write_cells_by_genes_csv_gz(me3_gxc, me3_prefix)
write_mtx_gz(me3_gxc, me3_prefix)

print("Done.")
print("Outputs in:", OUTDIR)