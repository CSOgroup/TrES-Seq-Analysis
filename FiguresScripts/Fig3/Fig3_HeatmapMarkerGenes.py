#!/usr/bin/env python3
# rna_marker_heatmap_and_tss.py
#
# - Reads an AnnData (.h5ad)
# - Filters genes (drop ENS*, miRNA, unknown/LOC)
# - Selects per-group markers from rGG (adj p < alpha, logFC > min_logfc, score > 0)
# - Sorts by score (or logFC) and enforces uniqueness across groups
# - Plots a heatmap (z-scored per gene) of those markers
# - Writes per-group TSS BED6 files for the *same* markers, in the *same* order
#
# Example:
#   python rna_marker_heatmap_and_tss.py \
#       --input AllCellLines_rna_merged_processed_Filtered.h5ad \
#       --output ./figures/CL_marker_heatmap.pdf \
#       --bed-file ../../human_epdnew_mjdVA.bed \
#       --tss-outdir tss_beds_topmarkers \
#       --top-n 50 \
#       --rgg-key CL_rgg \
#       --min-logfc 1 \
#       --groupby CL \
#       --layer log1p_Norm_counts

import os
import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
import re

# Desired display order for CL groups (cell lines)
CL_PREFERRED_ORDER = ["GM12878", "JJN2", "Karpas", "Karpas422", "WSU"]

# ------------------------------ helpers ------------------------------
def _to_dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)

def _infer_symbol_col(adata):
    """Try to find the column in .var that holds gene symbols."""
    candidates = ['gene_symbol','gene_symbols','symbol','SYMBOL','feature_name','GeneSymbol']
    for c in candidates:
        if c in adata.var.columns:
            return c
    return None

def _normalize_symbols_series(s):
    """Turn 'NAALADL2,NAALADL2-AS2' -> 'NAALADL2' (strip to first token)."""
    return s.astype(str).str.split(',').str[0].str.strip()

def _reorder_categories(adata, groupby, preferred_order):
    """
    Reorder categorical groups so that preferred_order (if present) appears first
    in that order, and any other categories follow.
    """
    if groupby not in adata.obs:
        return
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    cats = list(adata.obs[groupby].cat.categories)
    new_order = [g for g in preferred_order if g in cats] + [g for g in cats if g not in preferred_order]
    adata.obs[groupby] = adata.obs[groupby].cat.reorder_categories(new_order, ordered=True)

def _ensure_group_colors(adata, groupby, palette=None):
    if groupby not in adata.obs:
        return
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')

    # If this is CL, enforce preferred ordering before assigning colors
    if groupby == "CL":
        _reorder_categories(adata, groupby, CL_PREFERRED_ORDER)

    cats = list(adata.obs[groupby].cat.categories)
    if palette is None:
        palette = {
            'GM12878': '#1b9e77',
            'JJN2': '#d95f02',
            'Karpas': '#7570b3',
            'Karpas422': '#7570b3',
            'WSU': '#e7298a',
        }
    colors = [palette.get(c, '#808080') for c in cats]
    key = f'{groupby}_colors'
    if key not in adata.uns or len(adata.uns[key]) != len(cats):
        adata.uns[key] = colors

def build_gene_filter(adata, var_symbol_col=None):
    """
    Returns a boolean mask over variables keeping only annotated, non-ENS, non-miRNA genes.
    """
    symcol = var_symbol_col or _infer_symbol_col(adata)
    if symcol is not None:
        symbols = _normalize_symbols_series(adata.var[symcol])
    else:
        symbols = _normalize_symbols_series(pd.Series(adata.var_names, index=adata.var_names))

    # Try common biotype columns (optional)
    biotype_col = next((c for c in ['gene_biotype','biotype','feature_biotype','gene_type']
                        if c in adata.var.columns), None)
    biotypes = adata.var[biotype_col].astype(str).str.lower() if biotype_col else pd.Series("", index=symbols.index)

    mask = symbols.notna() & (symbols.str.len() > 0)
    mask &= ~symbols.str.match(r'^ENS[A-Z]*\d+', case=False)  # Ensembl-like IDs
    mask &= ~symbols.str.match(r'^(MIR|MIRLET)\w*', case=False)  # miRNA name patterns
    mask &= ~symbols.str.match(r'^(MT-)\w*', case=False)  # MT-RNA name patterns
    mask &= ~symbols.str.contains(r'\bmirna\b|\bmicroRNA\b', case=False, regex=True)
    mask &= ~symbols.str.contains(r'unknown|uncharacteri[sz]ed', case=False, regex=True)
    mask &= ~symbols.str.match(r'^LOC\d+$', case=False)
    mask &= ~biotypes.str.contains(r'\bmirna\b|\bmicroRNA\b', case=False, regex=True)
    return mask

def _valid_gene_universe(adata, var_symbol_col=None):
    """
    Returns the set of *symbols* that are valid after filtering, regardless of what var_names are.
    Uses normalized (comma-stripped) symbols.
    """
    mask = build_gene_filter(adata, var_symbol_col=var_symbol_col)
    symcol = var_symbol_col or _infer_symbol_col(adata)
    if symcol is None:
        return set(_normalize_symbols_series(pd.Index(adata.var_names)[mask]))
    else:
        return set(_normalize_symbols_series(adata.var.loc[mask, symcol]).tolist())

# ------------------------------ marker selection ------------------------------
def top_markers_from_rgg(
    adata, rgg_key='CL_rgg', groupby='CL', top_n=10,
    alpha=0.05, min_logfc=0.25, ensure_unique=True, use_score=True,
    var_symbol_col=None
):
    """
    Return dict[group] -> list of marker symbols that pass:
      - adjusted p < alpha
      - logFC > min_logfc
      - score > 0
    Sorted by 'score' (default) or 'logFC'. Enforces uniqueness across groups if requested.
    Only returns genes in the filtered, annotated universe.
    """
    if rgg_key not in adata.uns:
        raise KeyError(f"{rgg_key!r} not found in adata.uns")
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    # Enforce CL order if appropriate
    if groupby == "CL":
        _reorder_categories(adata, groupby, CL_PREFERRED_ORDER)

    rgg = adata.uns[rgg_key]
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    groups = list(adata.obs[groupby].cat.categories)
    fields = list(rgg['names'].dtype.names)
    missing = [g for g in groups if g not in fields]
    if missing:
        raise ValueError(f"Groups {missing} not present in uns[{rgg_key!r}] fields {fields}")

    valid_symbols = _valid_gene_universe(adata, var_symbol_col=var_symbol_col)

    used = set()
    markers = {}

    for g in groups:
        names_g  = np.array(rgg['names'][g], dtype=str)
        pvals_g  = np.array(rgg['pvals_adj'][g], dtype=float)
        lfc_g    = np.array(rgg['logfoldchanges'][g], dtype=float)
        scores_g = np.array(rgg['scores'][g], dtype=float)

        # normalize rGG names to first token before comma
        names_g = np.array(_normalize_symbols_series(pd.Series(names_g)), dtype=str)

        df = pd.DataFrame({'gene': names_g, 'p_adj': pvals_g, 'logFC': lfc_g, 'score': scores_g})
        df = df[df['gene'].isin(valid_symbols)]
        df = df[(df['p_adj'] < alpha) & (df['logFC'] > min_logfc) & (df['score'] > 0)]
        df = df.sort_values('score' if use_score else 'logFC', ascending=False)

        if ensure_unique and used:
            df = df[~df['gene'].isin(used)]

        picks = df['gene'].head(top_n).tolist()
        markers[g] = picks
        used.update(picks)

    return markers

# ------------------------------ printing stats ------------------------------
def print_selected_marker_stats(
    adata,
    markers_dict,
    rgg_key='CL_rgg',
    groupby='CL',
):
    """
    Print the marker genes actually used (after all filtering) + their stats (adj p, logFC, score).
    markers_dict is dict[group] -> list of genes in the final order used by heatmap/TSS.
    """
    rgg = adata.uns[rgg_key]
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    groups = list(adata.obs[groupby].cat.categories)

    print("\n=== Selected marker genes per group (after BED filtering) ===")
    out_rows = []
    for g in groups:
        genes = markers_dict.get(g, [])
        if not genes:
            continue
        print(f"\n{g} ({len(genes)} genes):")

        names_g  = np.array(rgg['names'][g], dtype=str)
        logfc_g  = np.array(rgg['logfoldchanges'][g], dtype=float)
        padj_g   = np.array(rgg['pvals_adj'][g], dtype=float)
        score_g  = np.array(rgg['scores'][g], dtype=float)

        # normalize name list for lookup consistency (handles "A,B" -> "A")
        names_g = np.array(_normalize_symbols_series(pd.Series(names_g)), dtype=str)

        df_g = (
            pd.DataFrame({'gene': names_g, 'logFC': logfc_g, 'p_adj': padj_g, 'score': score_g})
            .dropna(subset=['gene'])
            .astype({'logFC': float, 'p_adj': float, 'score': float})
            .set_index('gene', drop=False)
        )

        for gene in genes:
            if gene not in df_g.index:
                continue
            rows = df_g.loc[gene]
            if isinstance(rows, pd.DataFrame):
                rows = rows.sort_values(
                    by=['score', 'logFC', 'p_adj'],
                    ascending=[False, False, True],
                ).iloc[0]
            lgfc  = float(rows['logFC'])
            padj  = float(rows['p_adj'])
            score = float(rows['score'])

            print(f"  {gene:<15} logFC={lgfc:6.2f}  p_adj={padj:.2e}  score={score:6.2f}")
            out_rows.append({'group': g, 'gene': gene, 'logFC': lgfc, 'p_adj': padj, 'score': score})

    return pd.DataFrame(out_rows)

# ------------------------------ heatmap ------------------------------
def plot_marker_heatmap_pretty(
    adata,
    markers_dict,
    outfile,
    groupby='CL',
    layer='log1p_Norm_counts',
    percentile_clip=90.0,
    cmap='RdBu_r',
    base_fontsize=6,
    var_symbol_col=None,
):
    """
    Saves a marker heatmap to `outfile` (PDF).

    markers_dict: dict[group] -> list of genes in final order (already filtered).
    """
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')

    # Enforce CL order if appropriate
    if groupby == "CL":
        _reorder_categories(adata, groupby, CL_PREFERRED_ORDER)

    _ensure_group_colors(adata, groupby)

    markers_dict = {k: v for k, v in markers_dict.items() if v}
    if not markers_dict:
        raise ValueError("No markers passed filters; relax thresholds, increase top_n, or check BED coverage.")

    groups = [g for g in adata.obs[groupby].cat.categories if g in markers_dict]
    ordered_markers = {g: markers_dict[g] for g in groups}
    union_genes = [gene for grp in groups for gene in ordered_markers[grp]]  # symbols in display order

    # 2) subset columns by symbol (or by var_names if already symbols)
    symcol = var_symbol_col or _infer_symbol_col(adata)
    if symcol is not None:
        # normalize the symbol column for matching (strip anything after comma)
        symbol_series_norm = _normalize_symbols_series(adata.var[symcol])
        take = symbol_series_norm.isin(union_genes).values
        a_sub = adata[:, take].copy()

        # order columns to match union_genes (deduplicate if multiple map to same normalized name)
        sym_vals_norm = _normalize_symbols_series(a_sub.var[symcol]).values
        seen = set()
        order = []
        for g in union_genes:
            idxs = np.where(sym_vals_norm == g)[0]
            if len(idxs):
                idx = int(idxs[0])
                if idx not in seen:
                    order.append(idx)
                    seen.add(idx)

        a_sub = a_sub[:, order]

        # Rename var_names to normalized symbols (no commas)
        a_sub.var_names = _normalize_symbols_series(a_sub.var[symcol]).values
    else:
        # assume var_names are symbols; normalize them to drop commas too
        vn = _normalize_symbols_series(pd.Series(adata.var_names, index=adata.var_names))
        take = vn.isin(union_genes).values
        a_sub = adata[:, take].copy()
        # order columns to match union_genes
        vn_sub = _normalize_symbols_series(pd.Series(a_sub.var_names))
        order = []
        seen = set()
        for g in union_genes:
            idxs = np.where(vn_sub.values == g)[0]
            if len(idxs):
                idx = int(idxs[0])
                if idx not in seen:
                    order.append(idx)
                    seen.add(idx)
        a_sub = a_sub[:, order]
        a_sub.var_names = _normalize_symbols_series(pd.Series(a_sub.var_names)).values

    # 3) matrix and z-score per gene (columns)
    if layer is None:
        M = _to_dense(a_sub.X)
    else:
        if layer not in a_sub.layers:
            raise KeyError(f"Layer {layer!r} not found in adata.layers")
        M = _to_dense(a_sub.layers[layer])

    mu = M.mean(axis=0, dtype=float)
    sd = M.std(axis=0, dtype=float); sd[sd == 0] = 1.0
    Mz = (M - mu) / sd
    a_sub.layers['__scaled__'] = Mz

    # 4) vmin/vmax from percentile of |Z|
    vmax = float(np.percentile(np.abs(Mz.ravel()), percentile_clip))
    vmin = -vmax

    # 5) draw with Scanpy (no autosave; we save manually)
    sc.set_figure_params(dpi=200, facecolor='white', fontsize=base_fontsize)

    n_genes = len(union_genes)
    #fig_width = max(8, 0.1 * n_genes + 4)  # each gene takes less horizontal space
    fig_height = max(8, 0.02 * n_genes + 4)  # each gene takes vertical space now

    axgrid = sc.pl.heatmap(
        a_sub,
        var_names=ordered_markers,   # dict: group -> list of symbols
        groupby=groupby,
        layer='__scaled__',
        use_raw=False,
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        show=False,
        show_gene_labels=True,
        swap_axes=True,              # <<< FLIP AXES
        figsize=(9, fig_height),     # <<< width fixed, height scales with #genes
        dendrogram=False,
        var_group_rotation=0
    )

    # 6) tidy labels + colorbar
    fig = plt.gcf()
    for ax in fig.axes:
        if ax.get_ylabel() == groupby:
            ax.set_ylabel("")

    cbar = None
    for ax in fig.axes[::-1]:
        if ax.get_label().startswith("colorbar") or "Colorbar" in str(type(ax)):
            cbar = ax
            break
    if cbar is not None:
        cbar.set_ylabel("RNA Expression Z-score", rotation=90, va='center', fontsize=base_fontsize)
        cbar.tick_params(labelsize=base_fontsize-1, length=3)

    fig.tight_layout(pad=0.6)
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    return outfile

# ------------------------------ TSS BED creation ------------------------------
def make_tss_beds_for_top_markers(
    adata,
    markers_dict,
    bed_df,
    outdir="tss_beds_topmarkers",
    groupby="CL",
):
    """
    Create per-group TSS BED6 files for the given markers_dict.

    markers_dict: dict[group] -> list of genes (already filtered to genes with TSS).
    bed_df: BED DataFrame with columns [chrom, start, end, bed_name, score, strand, gene_name].
    """
    os.makedirs(outdir, exist_ok=True)

    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')

    # Same group order as in heatmap / markers selection
    groups_in_order = [
        g for g in adata.obs[groupby].cat.categories
        if g in markers_dict and markers_dict[g]
    ]

    for group in groups_in_order:
        genes_in_order = markers_dict[group]  # SAME order as heatmap

        # Subset BED to these genes only
        bed_sub = bed_df[bed_df["gene_name"].isin(genes_in_order)].copy()

        # Enforce gene order as in markers_dict (heatmap)
        bed_sub["gene_name"] = pd.Categorical(
            bed_sub["gene_name"],
            categories=genes_in_order,
            ordered=True,
        )
        bed_sub = bed_sub.sort_values("gene_name")

        # Ensure integer coordinates
        bed_sub["start"] = bed_sub["start"].astype(int)
        bed_sub["end"] = bed_sub["end"].astype(int)

        # BED6 output: chrom, start, end, gene_name, score, strand
        bed6_df = bed_sub[["chrom", "start", "end", "gene_name", "strand"]].copy()
        # Insert dummy score column (0); adjust if you want rGG scores instead
        bed6_df.insert(4, "score", 0)

        out_path = os.path.join(
            outdir,
            f"tss_top{len(genes_in_order)}_{group.replace(' ', '_')}.bed",
        )
        bed6_df.to_csv(out_path, sep="\t", header=False, index=False)
        print(f"✔️  Saved {len(bed6_df)} TSS entries → {out_path}")

    print(f"\n✅ All TSS BEDs written to: {outdir}")

# ------------------------------ main ------------------------------
def main():
    p = argparse.ArgumentParser(
        description="RNA marker heatmap + matching TSS BEDs (using the same marker sets and order)."
    )
    p.add_argument('--input', required=True, help='Input .h5ad file')
    p.add_argument('--output', required=True, help='Output PDF file for the heatmap')
    p.add_argument('--bed-file', required=True,
                   help='EPDnew-style BED file: chrom, start, end, bed_name, score, strand')
    p.add_argument('--tss-outdir', default='tss_beds_topmarkers',
                   help='Output directory for TSS BED files (default: tss_beds_topmarkers)')
    p.add_argument('--rgg-key', default='CL_rgg',
                   help='Key in adata.uns with rank_genes_groups-like struct (default: CL_rgg)')
    p.add_argument('--groupby', default='CL',
                   help='obs column for groups (default: CL)')
    p.add_argument('--layer', default='log1p_Norm_counts',
                   help='Layer to plot (default: log1p_Norm_counts). Use "None" for X.')
    p.add_argument('--top-n', type=int, default=40,
                   help='Top N markers per group after filtering (default: 40)')
    p.add_argument('--alpha', type=float, default=0.05,
                   help='Adjusted p-value threshold (default: 0.05)')
    p.add_argument('--min-logfc', type=float, default=0.25,
                   help='Minimum logFC (default: 0.25)')
    p.add_argument('--unique', action='store_true',
                   help='Enforce unique markers across groups')
    p.add_argument('--no-unique', dest='unique', action='store_false',
                   help='Allow overlap across groups')
    p.set_defaults(unique=True)
    p.add_argument('--use-score', action='store_true',
                   help='Sort by score (default)')
    p.add_argument('--use-logfc', dest='use_score', action='store_false',
                   help='Sort by logFC instead of score')
    p.set_defaults(use_score=True)
    p.add_argument('--percentile-clip', type=float, default=90.0,
                   help='Percentile for |Z| clipping (default: 90)')
    p.add_argument('--cmap', default='RdBu_r',
                   help='Matplotlib colormap (default: RdBu_r)')
    p.add_argument('--base-fontsize', type=float, default=6,
                   help='Base font size (default: 11)')
    args = p.parse_args()

    # --- Load data ---
    adata = sc.read_h5ad(args.input)

    # Enforce desired CL order globally if groupby is CL
    if args.groupby == "CL":
        _reorder_categories(adata, args.groupby, CL_PREFERRED_ORDER)

    symcol = _infer_symbol_col(adata)
    mask = build_gene_filter(adata, var_symbol_col=symcol)
    kept = int(mask.sum())
    total = int(adata.n_vars)
    print(f"[info] Gene filter: keeping {kept:,} / {total:,} variables")

    # --- Load BED and compute allowed genes (for both heatmap + TSS) ---
    bed_df = pd.read_csv(
        args.bed_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "bed_name", "score", "strand"],
    )
    # Strip suffixes like _1, _2 to get gene_name
    bed_df["gene_name"] = bed_df["bed_name"].str.replace(r"_\d+$", "", regex=True)
    allowed_genes = set(bed_df["gene_name"].astype(str))
    print(f"[info] BED: {len(allowed_genes):,} unique gene_name entries with TSS.")

    # --- Select markers (full set) ---
    raw_markers = top_markers_from_rgg(
        adata,
        rgg_key=args.rgg_key,
        groupby=args.groupby,
        top_n=args.top_n,
        alpha=args.alpha,
        min_logfc=args.min_logfc,
        ensure_unique=args.unique,
        use_score=args.use_score,
        var_symbol_col=symcol,
    )

    # --- Filter markers to genes that have TSS in BED (so heatmap == TSS) ---
    markers_dict = {}
    for g, genes in raw_markers.items():
        if not genes:
            markers_dict[g] = []
            continue
        filtered = [gene for gene in genes if gene in allowed_genes]
        missing = [gene for gene in genes if gene not in allowed_genes]
        if missing:
            print(
                f"⚠️ Group {g}: {len(missing)} marker genes have no TSS in BED "
                f"(will be dropped from heatmap + TSS): {', '.join(missing[:10])}"
                + ("..." if len(missing) > 10 else "")
            )
        markers_dict[g] = filtered

    # Drop groups with no remaining markers
    markers_dict = {g: genes for g, genes in markers_dict.items() if genes}

    if not markers_dict:
        raise ValueError("After BED filtering, no markers remain. "
                         "Check BED coverage, relax thresholds, or reduce top_n.")

    # --- Stats printout (uses final markers_dict) ---
    _ = print_selected_marker_stats(
        adata,
        markers_dict=markers_dict,
        rgg_key=args.rgg_key,
        groupby=args.groupby,
    )

    # --- Heatmap ---
    layer = None if (args.layer is None or args.layer.strip().lower() == 'none') else args.layer
    out = plot_marker_heatmap_pretty(
        adata,
        markers_dict=markers_dict,
        outfile=args.output,
        groupby=args.groupby,
        layer=layer,
        percentile_clip=args.percentile_clip,
        cmap=args.cmap,
        base_fontsize=args.base_fontsize,
        var_symbol_col=symcol,
    )
    print(f"[done] Heatmap written to {out}")

    # --- TSS BEDs (same markers, same order) ---
    make_tss_beds_for_top_markers(
        adata,
        markers_dict=markers_dict,
        bed_df=bed_df,
        outdir=args.tss_outdir,
        groupby=args.groupby,
    )

if __name__ == "__main__":
    sys.exit(main())



python rna_marker_heatmap_and_tss.py \
  --input AllCellLines_rna_merged_processed_Filtered.h5ad \
  --output ./figures/CL_marker_heatmap.pdf \
  --bed-file ../../human_epdnew_mjdVA.bed \
  --tss-outdir tss_beds_topmarkers \
  --top-n 50 \
  --rgg-key CL_rgg \
  --min-logfc 1 \
  --groupby CL \
  --layer log1p_Norm_counts
