#!/usr/bin/env python3
"""Plot the RNA markers and export their TSS intervals in the same order."""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get("DATA_ROOT", SCRIPT_DIR / "inputs")).expanduser()
WORK_DIR = Path(os.environ.get("WORK_DIR", SCRIPT_DIR / "results")).expanduser()
RNA_H5AD = WORK_DIR / "AllCellLines_rna_merged_processed_Filtered.h5ad"
EPD_BED = Path(os.environ.get("EPD_BED", DATA_ROOT / "human_epdnew_mjdVA.bed")).expanduser()
HEATMAP = WORK_DIR / "figures/CL_marker_heatmap.pdf"
TSS_BED_DIR = WORK_DIR / "tss_beds_topmarkers"
CL_PREFERRED_ORDER = ["GM12878", "JJN3", "Karpas", "Karpas422", "WSU"]


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


def _gene_symbols(adata, var_symbol_col=None):
    column = var_symbol_col or _infer_symbol_col(adata)
    values = adata.var[column] if column is not None else pd.Series(adata.var_names, index=adata.var_names)
    return _normalize_symbols_series(values)

def _reorder_categories(adata, groupby, preferred_order):
    """
    Reorder categorical groups so that preferred_order (if present) appears first
    in that order, and any other categories follow.
    """
    if groupby not in adata.obs:
        return
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    cats = list(adata.obs[groupby].cat.categories)
    new_order = [g for g in preferred_order if g in cats] + [g for g in cats if g not in preferred_order]
    adata.obs[groupby] = adata.obs[groupby].cat.reorder_categories(new_order, ordered=True)

def _ensure_group_colors(adata, groupby, palette=None):
    if groupby not in adata.obs:
        return
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        adata.obs[groupby] = adata.obs[groupby].astype('category')

    # If this is CL, enforce preferred ordering before assigning colors
    if groupby == "CL":
        _reorder_categories(adata, groupby, CL_PREFERRED_ORDER)

    cats = list(adata.obs[groupby].cat.categories)
    if palette is None:
        palette = {
            'GM12878': '#1b9e77',
            'JJN3': '#d95f02',
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
    symbols = _gene_symbols(adata, var_symbol_col)

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
    return set(_gene_symbols(adata, var_symbol_col)[mask])


def _ranked_genes(rgg, group):
    return pd.DataFrame({
        "gene": _normalize_symbols_series(pd.Series(rgg["names"][group], dtype=str)),
        "p_adj": np.asarray(rgg["pvals_adj"][group], dtype=float),
        "logFC": np.asarray(rgg["logfoldchanges"][group], dtype=float),
        "score": np.asarray(rgg["scores"][group], dtype=float),
    })

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
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
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
        df = _ranked_genes(rgg, g)
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
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        adata.obs[groupby] = adata.obs[groupby].astype('category')
    groups = list(adata.obs[groupby].cat.categories)

    print("\n=== Selected marker genes per group (after BED filtering) ===")
    for g in groups:
        genes = markers_dict.get(g, [])
        if not genes:
            continue
        print(f"\n{g} ({len(genes)} genes):")

        df_g = _ranked_genes(rgg, g).set_index('gene', drop=False)

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
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        adata.obs[groupby] = adata.obs[groupby].astype('category')

    _ensure_group_colors(adata, groupby)

    markers_dict = {k: v for k, v in markers_dict.items() if v}
    if not markers_dict:
        raise ValueError("No markers passed filters; relax thresholds, increase top_n, or check BED coverage.")

    groups = [g for g in adata.obs[groupby].cat.categories if g in markers_dict]
    ordered_markers = {g: markers_dict[g] for g in groups}
    union_genes = [gene for grp in groups for gene in ordered_markers[grp]]  # symbols in display order

    symbols = _gene_symbols(adata, var_symbol_col)
    first_column = {}
    for index, symbol in enumerate(symbols):
        first_column.setdefault(symbol, index)
    order = [first_column[g] for g in dict.fromkeys(union_genes) if g in first_column]
    a_sub = adata[:, order].copy()
    a_sub.var_names = symbols.iloc[order].to_numpy()

    # Scale each marker across cells.
    if layer is None:
        M = _to_dense(a_sub.X)
    else:
        if layer not in a_sub.layers:
            raise KeyError(f"Layer {layer!r} not found in adata.layers")
        M = _to_dense(a_sub.layers[layer])

    mu = M.mean(axis=0, dtype=float)
    sd = M.std(axis=0, dtype=float)
    sd[sd == 0] = 1.0
    Mz = (M - mu) / sd
    a_sub.layers['__scaled__'] = Mz

    vmax = float(np.percentile(np.abs(Mz.ravel()), percentile_clip))
    vmin = -vmax

    sc.set_figure_params(dpi=200, facecolor='white', fontsize=base_fontsize)

    n_genes = len(union_genes)
    fig_height = max(8, 0.02 * n_genes + 4)

    sc.pl.heatmap(
        a_sub,
        var_names=ordered_markers,
        groupby=groupby,
        layer='__scaled__',
        use_raw=False,
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        show=False,
        show_gene_labels=True,
        swap_axes=True,
        figsize=(9, fig_height),
        dendrogram=False,
        var_group_rotation=0
    )

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
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
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
    p.add_argument('--input', type=Path, default=RNA_H5AD, help='Input .h5ad file')
    p.add_argument('--output', type=Path, default=HEATMAP, help='Output PDF file for the heatmap')
    p.add_argument('--bed-file', type=Path, default=EPD_BED,
                   help='EPDnew-style BED file: chrom, start, end, bed_name, score, strand')
    p.add_argument('--tss-outdir', type=Path, default=TSS_BED_DIR,
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
                   help='Base font size (default: 6)')
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
    print_selected_marker_stats(
        adata,
        markers_dict=markers_dict,
        rgg_key=args.rgg_key,
        groupby=args.groupby,
    )

    # --- Heatmap ---
    args.output.parent.mkdir(parents=True, exist_ok=True)
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
    main()
