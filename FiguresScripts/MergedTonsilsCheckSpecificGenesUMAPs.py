#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as snap
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

# ------------------- Config -------------------
RNA_PATH = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"
AC_PATH  = "DNA_merged/MergedTonsils_H3K27ac_merged_processed_barcodeRewritten_withRNAlabels.h5ad"
ME3_PATH = "DNA_merged/MergedTonsils_H3K27me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad"

OUTDIR = "Tonsil_SpecificGenes_RNA_AC_ME3_overlapOnly"
os.makedirs(OUTDIR, exist_ok=True)

# Genes of interest
GENES = ["BCL6", "IRF4", "XBP1", "SDC1", "MYC"]

# Grouping for plots/heatmaps
GROUPBY = "tonsil_id"
FALLBACK_GROUPBY = "leiden_merged_type"

# Gene-activity window ±bp
WINDOW = 2_000

# UMAP / plotting parameters
N_HVG = 3000
N_PCS = 20
N_NEIGHBORS = 30
METRIC = "cosine"

# UMAP/heatmap scaling helpers
MEDIAN_FACTOR  = 3.0
MEAN_FACTOR    = 2.0
CAP_PERCENTILE = 99.0
ZERO_EPS       = 1e-12

# RNA layer choice
RNA_LAYER_FOR_HEATMAP = "auto"
RNA_LAYER_PRIORITY = ["log1p_norm", "SCT_data", "SCT_counts", "counts", None]

# Binary rules for overlap heatmap ordering
RNA_BIN_RULE = "auto"
AC_BIN_RULE  = "ge1"
ME3_BIN_RULE = "ge1"

PATTERN_ORDER = [
    0b111,  # RNA+AC+ME3
    0b110,  # RNA+AC
    0b100,  # RNA only
    0b010,  # AC only
    0b001,  # ME3 only
    0b011,  # AC+ME3
    0b101,  # RNA+ME3
    0b000,  # none
]

try:
    GENE_ANNO = snap.genome.hg38
except Exception:
    GENE_ANNO = None


# ------------------- I/O -------------------
adata_rna = sc.read_h5ad(RNA_PATH)
adata_ac  = sc.read_h5ad(AC_PATH)
adata_me3 = sc.read_h5ad(ME3_PATH)

# exact overlap of RNA / H3K27ac / H3K27me3
common = (
    adata_rna.obs_names
    .intersection(adata_ac.obs_names)
    .intersection(adata_me3.obs_names)
)

print(f"[INFO] Common cells across RNA / H3K27ac / H3K27me3: {len(common)}")

adata_rna = adata_rna[common].copy()
adata_ac  = adata_ac[common].copy()
adata_me3 = adata_me3[common].copy()


# ------------------- Helpers -------------------
def ensure_counts_layer_for_rna(adata):
    if "counts" in adata.layers:
        return "counts", None
    if adata.raw is not None:
        return None, True
    return None, False


def ensure_log1p_layer_from(adata, prefer_layer="counts", layer_name="log1p"):
    if layer_name in adata.layers:
        return layer_name
    Xsrc = adata.layers[prefer_layer] if prefer_layer in adata.layers else adata.X
    if sp.issparse(Xsrc):
        Xlog = Xsrc.copy()
        Xlog.data = np.log1p(Xlog.data)
    else:
        Xlog = np.log1p(Xsrc)
    adata.layers[layer_name] = Xlog
    return layer_name


def resolve_rna_layer(adata, choice="auto", priority=None):
    if choice != "auto":
        if choice in adata.layers:
            return choice
        print(f"[WARN] RNA layer '{choice}' not found; falling back to auto.")
    priority = (priority or []) + [None]
    for lay in priority:
        if lay is None:
            return None
        if lay in adata.layers:
            return lay
    return None


def prepare_rna_for_embedding(adata):
    adata = adata.copy()
    if "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata


def recompute_rna_umap(adata):
    adata = prepare_rna_for_embedding(adata)
    hvg_kwargs = dict(flavor="seurat_v3", n_top_genes=N_HVG, subset=False)
    if "counts" in adata.layers:
        hvg_kwargs["layer"] = "counts"
    if "sample" in adata.obs.columns:
        hvg_kwargs["batch_key"] = "sample"

    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    sc.pp.scale(adata, zero_center=False, max_value=4)
    sc.pp.pca(adata, n_comps=N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=N_NEIGHBORS, metric=METRIC)
    sc.tl.umap(adata)
    return adata


def copy_umap_from_reference(ref, target, key="X_umap"):
    common_local = ref.obs_names.intersection(target.obs_names)
    if len(common_local) == 0:
        raise ValueError("No shared barcodes between reference and target.")
    target = target[common_local].copy()
    ref_umap = pd.DataFrame(ref.obsm[key], index=ref.obs_names, columns=["UMAP1", "UMAP2"])
    target.obsm[key] = ref_umap.loc[target.obs_names].to_numpy()
    return target


def add_obs_from_rna(reference_rna, target, columns=("sample", "tonsil_id", "leiden_merged_type", "leiden")):
    for col in columns:
        if col in reference_rna.obs:
            target.obs[col] = reference_rna.obs[col].reindex(target.obs_names)


def pick_groupby(adata):
    if GROUPBY in adata.obs:
        return GROUPBY
    if FALLBACK_GROUPBY in adata.obs:
        return FALLBACK_GROUPBY
    if "leiden" in adata.obs:
        return "leiden"
    adata.obs["_all"] = "all"
    return "_all"


def upper_intersection(list_like, var_names):
    want = {g.upper() for g in list_like}
    mapping = {v.upper(): v for v in var_names}
    return [mapping[g] for g in want if g in mapping]


def _resolve_anno_for_snap(gene_anno):
    if gene_anno is None:
        raise ValueError("Set GENE_ANNO to a hg38 GTF/GFF path or ensure snap.genome.hg38 is available.")
    return getattr(gene_anno, "annotation", gene_anno)


def make_gene_activity(adata_mark, gene_anno, window_bp):
    ann = _resolve_anno_for_snap(gene_anno)
    ga = snap.pp.make_gene_matrix(
        adata_mark,
        gene_anno=ann,
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
    ga.layers["log1p"] = Xlog
    return ga


def make_zero_gray_viridis():
    base = mpl.cm.get_cmap("viridis", 256)
    colors = base(np.linspace(0, 1, 256))
    cmap = ListedColormap(colors)
    cmap.set_under("lightgray")
    return cmap


ZERO_GRAY_VIRIDIS = make_zero_gray_viridis()


def _get_gene_values(adata, gene, layer=None):
    try:
        j = adata.var_names.get_loc(gene)
    except KeyError:
        return np.array([])
    Xsrc = adata.layers[layer] if layer is not None else adata.X
    return Xsrc.getcol(j).toarray().ravel() if sp.issparse(Xsrc) else np.asarray(Xsrc[:, j]).ravel()


def _smart_vmax(values, med_factor=MEDIAN_FACTOR, mean_factor=MEAN_FACTOR, cap_pct=CAP_PERCENTILE):
    if values.size == 0:
        return None
    nz = values[values > 0]
    if nz.size == 0:
        return None
    med = np.median(nz)
    vmax_candidate = med_factor * med if med > 0 else (mean_factor * np.mean(nz) if np.mean(nz) > 0 else None)
    if vmax_candidate is None:
        return None
    cap = np.percentile(nz, cap_pct)
    return float(max(1e-9, min(vmax_candidate, cap))) if np.isfinite(cap) else float(max(1e-9, vmax_candidate))


def save_umap(adata, color, fname, layer=None, vmin=0, vmax=None, color_map=None):
    sc.pl.umap(
        adata,
        color=color,
        layer=layer,
        frameon=False,
        title=fname.replace(".pdf", "").replace(".png", ""),
        vmin=vmin,
        vmax=vmax,
        color_map=color_map,
        show=False,
        save=None,
    )
    plt.savefig(os.path.join(OUTDIR, fname), bbox_inches="tight", dpi=200)
    plt.close()


def save_gene_umaps(adata, genes, layer, prefix, pdf=True):
    present = upper_intersection(genes, adata.var_names)
    missing = [g for g in genes if g.upper() not in {p.upper() for p in present}]
    if missing:
        print(f"[WARN] Missing genes: {missing}")
    for g in present:
        vals = _get_gene_values(adata, g, layer=layer)
        vmax = _smart_vmax(vals)
        vmin_use, vmax_use = (ZERO_EPS, ZERO_EPS) if vmax is None else (ZERO_EPS, vmax)
        fn = f"{prefix}_{g}_{(layer or 'X')}.{'pdf' if pdf else 'png'}"
        save_umap(adata, color=g, layer=layer, fname=fn, vmin=vmin_use, vmax=vmax_use, color_map=ZERO_GRAY_VIRIDIS)


def _binarize_vector(v: np.ndarray, mode: str, layer_name: str | None = None) -> np.ndarray:
    m = (mode or "auto").lower()
    if m == "auto":
        if layer_name and (layer_name.lower().startswith("sct_")
                           or layer_name.lower().startswith("log1p")
                           or layer_name.lower().startswith("logcounts")):
            return v > 0
        return v >= 1
    if m == "gt0":
        return v > 0
    if m == "ge1":
        return v >= 1
    if m == "ge2":
        return v >= 2
    if m == "z>1":
        mu = v.mean()
        sd = v.std()
        if sd == 0:
            return np.zeros_like(v, dtype=bool)
        return (v - mu) / sd > 1.0
    raise ValueError(f"Unknown binarization mode: {mode}")


def _ga_gene_vectors(ga_adata, gene, order_index=None):
    try:
        j = ga_adata.var_names.get_loc(gene)
    except KeyError:
        n = ga_adata.n_obs
        zeros = np.zeros(n, dtype=float)
        return (zeros if order_index is None else zeros[order_index],
                zeros if order_index is None else zeros[order_index])

    Xcnt = ga_adata.X
    v_cnt = Xcnt.getcol(j).toarray().ravel() if sp.issparse(Xcnt) else np.asarray(Xcnt[:, j]).ravel()
    Xlog = ga_adata.layers["log1p"]
    v_log = Xlog.getcol(j).toarray().ravel() if sp.issparse(Xlog) else np.asarray(Xlog[:, j]).ravel()

    if order_index is None:
        return v_cnt, v_log
    return v_cnt[order_index], v_log[order_index]


def _order_cells_by_pattern(
    rna_vals, ac_cnt, me3_cnt, group_index: np.ndarray,
    rna_layer_name: str | None,
    pattern_order=PATTERN_ORDER,
    rna_rule=RNA_BIN_RULE, ac_rule=AC_BIN_RULE, me3_rule=ME3_BIN_RULE
) -> np.ndarray:
    rna_g = rna_vals[group_index]
    ac_g = ac_cnt[group_index]
    me3_g = me3_cnt[group_index]

    b_rna = _binarize_vector(rna_g, rna_rule, layer_name=rna_layer_name)
    b_ac = _binarize_vector(ac_g, ac_rule, layer_name=None)
    b_me3 = _binarize_vector(me3_g, me3_rule, layer_name=None)

    pattern = (b_rna.astype(int) << 2) | (b_ac.astype(int) << 1) | (b_me3.astype(int) << 0)
    pri_map = {p: i for i, p in enumerate(pattern_order)}
    pri = np.array([pri_map.get(int(x), len(pattern_order)) for x in pattern], dtype=int)

    mean_sig = (rna_g + ac_g + me3_g) / 3.0
    local_order = np.lexsort((-mean_sig, pri))
    return group_index[local_order]


def _group_order_and_index(adata, groupby):
    if groupby not in adata.obs:
        groupby = pick_groupby(adata)
    ser = adata.obs[groupby].astype("category") if groupby in adata.obs else pd.Series(index=adata.obs_names, data="all")
    order = list(ser.cat.categories) if str(ser.dtype) == "category" else sorted(ser.unique().tolist())
    idx = []
    for g in order:
        idx.extend(list(np.where(adata.obs[groupby].astype(str).values == str(g))[0]))
    return groupby, order, np.array(idx, dtype=int)


def _scale_rows_zero_to_one(rows):
    scaled, vmax_used = [], []
    for v in rows:
        vmax = _smart_vmax(v)
        if not (vmax and np.isfinite(vmax) and vmax > 0):
            scaled.append(np.zeros_like(v))
            vmax_used.append(0.0)
        else:
            s = np.clip(v / vmax, 0.0, 1.0)
            scaled.append(s)
            vmax_used.append(float(vmax))
    return np.vstack(scaled), vmax_used


def plot_gene_heatmap_three_rows(
    gene, adata_rna, ga_ac, ga_me3, rna_layer_for_heatmap=None,
    groupby="tonsil_id", fname_prefix="Heatmap_byGroup"
):
    groupby, group_order, _base_index = _group_order_and_index(adata_rna, groupby)

    rna_layer = resolve_rna_layer(adata_rna, RNA_LAYER_FOR_HEATMAP, RNA_LAYER_PRIORITY) \
        if rna_layer_for_heatmap is None else rna_layer_for_heatmap

    try:
        j_rna = adata_rna.var_names.get_loc(gene)
        Xr = adata_rna.layers[rna_layer] if rna_layer is not None else adata_rna.X
        rna_vec = Xr.getcol(j_rna).toarray().ravel() if sp.issparse(Xr) else np.asarray(Xr[:, j_rna]).ravel()
    except KeyError:
        rna_vec = np.zeros(adata_rna.n_obs, dtype=float)

    pos = pd.Series(index=adata_rna.obs_names, data=np.arange(adata_rna.n_obs))
    ac_cnt_rna, ac_log_rna = _ga_gene_vectors(ga_ac, gene, order_index=pos.loc[ga_ac.obs_names].values)
    me3_cnt_rna, me3_log_rna = _ga_gene_vectors(ga_me3, gene, order_index=pos.loc[ga_me3.obs_names].values)

    ser_grp = adata_rna.obs[groupby].astype(str).values
    group_names = [g for g in group_order if np.any(ser_grp == str(g))]
    refined_indices = []
    for g in group_names:
        idx_g = np.where(ser_grp == str(g))[0]
        ord_g = _order_cells_by_pattern(
            rna_vals=rna_vec,
            ac_cnt=ac_cnt_rna,
            me3_cnt=me3_cnt_rna,
            group_index=idx_g,
            rna_layer_name=rna_layer,
            pattern_order=PATTERN_ORDER,
            rna_rule=RNA_BIN_RULE,
            ac_rule=AC_BIN_RULE,
            me3_rule=ME3_BIN_RULE,
        )
        refined_indices.append(ord_g)

    order_index = np.concatenate(refined_indices, axis=0) if len(refined_indices) else _base_index

    rna_v = rna_vec[order_index]
    ac_v = ac_log_rna[order_index]
    me3_v = me3_log_rna[order_index]

    M, vmaxs = _scale_rows_zero_to_one([rna_v, ac_v, me3_v])
    n = M.shape[1]
    fig_w = max(6.0, min(18.0, n / 80.0 + 6.0))
    fig_h = 2.7

    plt.figure(figsize=(fig_w, fig_h), dpi=300)
    ax = plt.gca()
    im = ax.imshow(M, aspect="auto", interpolation="nearest",
                   cmap=ZERO_GRAY_VIRIDIS, vmin=ZERO_EPS, vmax=1.0)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["RNA", "H3K27ac", "H3K27me3"])
    ax.set_xticks([])
    ax.set_title(gene, pad=6)

    counts, labels = [], []
    series_refined = adata_rna.obs[groupby].astype(str).values[order_index]
    for g in group_names:
        c = int((series_refined == str(g)).sum())
        if c > 0:
            counts.append(c)
            labels.append(str(g))
    edges = np.cumsum([0] + counts)
    for boundary in edges[1:-1]:
        ax.vlines(boundary - 0.5, -0.5, 2.5, color="white", lw=0.6, alpha=0.85)
    centers, left = [], 0
    for c in counts:
        centers.append(left + (c - 1) / 2.0)
        left += c
    top_ax = ax.secondary_xaxis("top")
    top_ax.set_xticks(centers)
    top_ax.set_xticklabels(labels)
    plt.setp(top_ax.get_xticklabels(), rotation=35, ha="left", va="bottom", fontsize=8)
    top_ax.tick_params(axis="x", pad=2, length=0)

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("scaled intensity")

    plt.tight_layout()
    png = os.path.join(OUTDIR, f"{fname_prefix}_{gene}_RNA_AC_ME3.png")
    pdf = os.path.join(OUTDIR, f"{fname_prefix}_{gene}_RNA_AC_ME3.pdf")
    plt.savefig(png, bbox_inches="tight")
    plt.savefig(pdf, bbox_inches="tight")
    plt.close()
    print(f"[heatmap] {gene}: vmax RNA={vmaxs[0]:.3g}, ac={vmaxs[1]:.3g}, me3={vmaxs[2]:.3g} | RNA layer='{rna_layer}'")


# ------------------- Recompute RNA UMAP on overlap cells -------------------
adata_rna = recompute_rna_umap(adata_rna)

# Optional QC plots
sc.pl.pca_variance_ratio(adata_rna, n_pcs=N_PCS, log=False, show=False)
plt.savefig(os.path.join(OUTDIR, "RNA_overlap_AC_ME3_PCAVar.pdf"), bbox_inches="tight", dpi=200)
plt.close()

if "sample" in adata_rna.obs:
    save_umap(adata_rna, color="sample", fname="UMAP_RNA_overlap_AC_ME3_sample.pdf")
if "tonsil_id" in adata_rna.obs:
    save_umap(adata_rna, color="tonsil_id", fname="UMAP_RNA_overlap_AC_ME3_tonsil.pdf")
if "leiden_merged_type" in adata_rna.obs:
    save_umap(adata_rna, color="leiden_merged_type", fname="UMAP_RNA_overlap_AC_ME3_celltype.pdf")


# ------------------- Project RNA UMAP to AC/ME3 overlap cells -------------------
for ad_mod in (adata_ac, adata_me3):
    add_obs_from_rna(adata_rna, ad_mod, columns=("sample", "tonsil_id", "leiden_merged_type", "leiden"))
    ad_mod = copy_umap_from_reference(adata_rna, ad_mod)

# keep returned objects
adata_ac = copy_umap_from_reference(adata_rna, adata_ac)
adata_me3 = copy_umap_from_reference(adata_rna, adata_me3)


# ------------------- Gene activity on overlap cells -------------------
ga_ac = make_gene_activity(adata_ac, GENE_ANNO, WINDOW)
ga_me3 = make_gene_activity(adata_me3, GENE_ANNO, WINDOW)

ga_ac = copy_umap_from_reference(adata_rna, ga_ac)
ga_me3 = copy_umap_from_reference(adata_rna, ga_me3)


# ------------------- UMAPs colored by genes -------------------
rna_layer_for_umap = resolve_rna_layer(adata_rna, RNA_LAYER_FOR_HEATMAP, RNA_LAYER_PRIORITY)
if rna_layer_for_umap is None:
    print("[INFO] RNA gene UMAPs will use X (no layer).")

save_gene_umaps(adata_rna, GENES, layer=rna_layer_for_umap, prefix="UMAP_RNA_overlap_AC_ME3", pdf=True)
save_gene_umaps(ga_ac,  GENES, layer="log1p", prefix="UMAP_H3K27ac_overlap_AC_ME3_RNAembed", pdf=True)
save_gene_umaps(ga_me3, GENES, layer="log1p", prefix="UMAP_H3K27me3_overlap_AC_ME3_RNAembed", pdf=True)


# ------------------- Violins -------------------
def _save_violin(adata, keys, groupby, fname, layer=None, use_raw=None, ylabel="counts"):
    present = upper_intersection(keys, adata.var_names)
    if len(present) == 0:
        print(f"[WARN] None of {keys} were found.")
        return
    sc.pl.violin(
        adata,
        keys=present,
        groupby=groupby,
        layer=layer,
        use_raw=use_raw,
        multi_panel=True,
        rotation=30,
        show=False,
    )
    plt.suptitle(f"{fname.replace('.png','')} — {ylabel}")
    plt.savefig(os.path.join(OUTDIR, fname), bbox_inches="tight", dpi=200)
    plt.close()


rna_counts_layer, rna_use_raw = ensure_counts_layer_for_rna(adata_rna)
groupby_vln_rna = "tonsil_id" if "tonsil_id" in adata_rna.obs else pick_groupby(adata_rna)
groupby_vln_ac = "tonsil_id" if "tonsil_id" in ga_ac.obs else pick_groupby(ga_ac)
groupby_vln_me3 = "tonsil_id" if "tonsil_id" in ga_me3.obs else pick_groupby(ga_me3)

_save_violin(
    adata_rna, GENES, groupby_vln_rna,
    fname="Violin_RNA_overlap_AC_ME3_counts_byTonsil.png",
    layer=rna_counts_layer, use_raw=rna_use_raw, ylabel="RNA counts"
)
_save_violin(
    ga_ac, GENES, groupby_vln_ac,
    fname=f"Violin_H3K27ac_overlap_AC_ME3_geneActivity_pm{WINDOW//1000}kb_byTonsil.png",
    ylabel=f"H3K27ac counts in gene body ±{WINDOW//1000}kb"
)
_save_violin(
    ga_me3, GENES, groupby_vln_me3,
    fname=f"Violin_H3K27me3_overlap_AC_ME3_geneActivity_pm{WINDOW//1000}kb_byTonsil.png",
    ylabel=f"H3K27me3 counts in gene body ±{WINDOW//1000}kb"
)


# ------------------- Heatmaps -------------------
HEATMAP_DIR = os.path.join(OUTDIR, "Heatmaps")
os.makedirs(HEATMAP_DIR, exist_ok=True)

rna_layer_for_heatmap = resolve_rna_layer(adata_rna, RNA_LAYER_FOR_HEATMAP, RNA_LAYER_PRIORITY)
for g in upper_intersection(GENES, adata_rna.var_names):
    plot_gene_heatmap_three_rows(
        g,
        adata_rna=adata_rna,
        ga_ac=ga_ac,
        ga_me3=ga_me3,
        rna_layer_for_heatmap=rna_layer_for_heatmap,
        groupby="tonsil_id",
        fname_prefix="Heatmap_overlap_AC_ME3_byTonsil",
    )

print(f"Done. Figures in: {OUTDIR} (heatmaps in {HEATMAP_DIR})")
