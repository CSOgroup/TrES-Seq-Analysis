#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from anndata import AnnData
from typing import Union
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors
from matplotlib_venn import venn3

os.makedirs("figures", exist_ok=True)

# ------------------------------------------------------------
# INPUTS
# ------------------------------------------------------------
RNA_FILE = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"
AC_FILE  = "DNA_merged/MergedTonsils_H3K27ac_merged_processed_barcodeRewritten_withRNAlabels.h5ad"
ME3_FILE = "DNA_merged/MergedTonsils_H3K27me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad"

CELLTYPE_KEY = "leiden_merged_type"
TONSIL_KEY = "tonsil_id"

CELL_TYPE_COLORS = {
    "B_Cycling": "#00BFC4",
    "B_GC_DZ": "#98DF8A",
    "B_GC_LZ": "#2CA02C",
    "B_mem": "#2B8C9E",
    "B_naive": "#1F77B4",
    "B_plasma": "#76C1FF",
    "DC": "#B39B00",
    "FDC": "#7B1FA2",
    "Macrophages": "#8C564B",
    "T_CD4+": "#D62728",
    "T_CD4+_TfH": "#E377C2",
    "T_CD4+_TfH_GC": "#C05AA0",
    "T_CD4+_naive": "#FF9896",
    "T_CD8+": "#FF7F0E",
}

CELLTYPE_ORDER = [
    "B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive", "B_plasma",
    "DC", "FDC", "Macrophages",
    "T_CD4+", "T_CD4+_TfH", "T_CD4+_TfH_GC", "T_CD4+_naive", "T_CD8+",
]

TONSIL_COLORS = {
    "Tonsil1": "#1f77b4",
    "Tonsil2": "#ff7f0e",
    "Tonsil3": "#2ca02c",
    "Tonsil4": "#d62728",
}

TONSIL_ORDER = ["Tonsil1", "Tonsil2", "Tonsil3", "Tonsil4"]

# If you want to recompute RNA UMAP on the exact overlap cells:
RECOMPUTE_RNA_UMAP = True
N_HVG = 3000
N_PCS = 20
N_NEIGHBORS = 30
RNA_METRIC = "cosine"


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------
def _ensure_adata(x: Union[str, AnnData]) -> AnnData:
    if isinstance(x, AnnData):
        return x
    if isinstance(x, str):
        return sc.read_h5ad(x)
    raise TypeError("Input must be an AnnData or a path to .h5ad")


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

    hvg_kwargs = dict(
        flavor="seurat_v3",
        n_top_genes=N_HVG,
        subset=False,
    )
    if "counts" in adata.layers:
        hvg_kwargs["layer"] = "counts"
    if "sample" in adata.obs.columns:
        hvg_kwargs["batch_key"] = "sample"

    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    sc.pp.scale(adata, zero_center=False, max_value=4)
    sc.pp.pca(adata, n_comps=N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=N_NEIGHBORS, metric=RNA_METRIC)
    sc.tl.umap(adata)
    return adata


def ensure_umap_exists(adata, use_rep_priority=("X_umap", "X_spectral", "X_pca"), n_neighbors=20):
    if "X_umap" in adata.obsm:
        return adata

    if "X_spectral" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_spectral", n_neighbors=n_neighbors, metric="cosine")
        sc.tl.umap(adata)
        return adata

    if "X_pca" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=n_neighbors, metric="cosine")
        sc.tl.umap(adata)
        return adata

    raise ValueError("No existing X_umap / X_spectral / X_pca found.")


def procrustes_align_to_base(target_xy, base_xy):
    X = np.asarray(base_xy, dtype=float)
    Y = np.asarray(target_xy, dtype=float)

    Xmean = X.mean(0, keepdims=True)
    Ymean = Y.mean(0, keepdims=True)
    Xc = X - Xmean
    Yc = Y - Ymean

    M = Yc.T @ Xc
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    R = U @ Vt
    s = S.sum() / (Yc ** 2).sum()

    return s * (Yc @ R) + Xmean


def extent(arr):
    arr = np.asarray(arr)
    return float(arr[:, 0].max() - arr[:, 0].min()), float(arr[:, 1].max() - arr[:, 1].min())


def ordered_categories(series, desired_order):
    vals = pd.Index(series.astype(str).unique())
    ordered = [x for x in desired_order if x in vals]
    extras = sorted([x for x in vals if x not in ordered])
    return ordered + extras


def make_legend(ax, palette_dict, order, title):
    handles = [
        Line2D([0], [0], marker="o", linestyle="None",
               markerfacecolor=palette_dict[k], markeredgecolor="none",
               markersize=7, label=k)
        for k in order if k in palette_dict
    ]
    leg = ax.legend(
        handles=handles,
        title=title,
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        frameon=False,
        fontsize=9,
    )
    leg.get_title().set_fontsize(10)


def build_threeway_overlap_objects(rna_file, ac_file, me3_file, recompute_rna=True):
    rna = _ensure_adata(rna_file).copy()
    ac = _ensure_adata(ac_file).copy()
    me3 = _ensure_adata(me3_file).copy()

    # Venn should be based on full objects before subsetting
    rna_set = set(rna.obs_names)
    ac_set = set(ac.obs_names)
    me3_set = set(me3.obs_names)

    common = rna.obs_names.intersection(ac.obs_names).intersection(me3.obs_names)
    if len(common) == 0:
        raise ValueError("No overlapping cells across RNA / H3K27ac / H3K27me3")

    rna_ov = rna[common].copy()
    ac_ov = ac[common].copy()
    me3_ov = me3[common].copy()

    if recompute_rna:
        rna_ov = recompute_rna_umap(rna_ov)
    else:
        rna_ov = ensure_umap_exists(rna_ov)

    ac_ov = ensure_umap_exists(ac_ov)
    me3_ov = ensure_umap_exists(me3_ov)

    return (rna, ac, me3, rna_set, ac_set, me3_set, rna_ov, ac_ov, me3_ov)


def plot_venn3(rna_set, ac_set, me3_set, output_prefix):
    plt.figure(figsize=(7, 7))
    venn3(
        [rna_set, ac_set, me3_set],
        set_labels=("RNA", "H3K27ac", "H3K27me3"),
    )
    plt.title("Cell overlap: RNA vs H3K27ac vs H3K27me3")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    plt.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close()


def plot_linked_umaps_three_modalities(
    adata_rna,
    adata_ac,
    adata_me3,
    color_field,
    palette_dict,
    desired_order,
    output_prefix,
    fig_size=(18, 6),
    fig_dpi=300,
    point_size=22,
    line_width=0.25,
    line_alpha=0.20,
    gap_frac=0.35,
    lift_rna_frac=0.25,
):
    common = adata_rna.obs_names.intersection(adata_ac.obs_names).intersection(adata_me3.obs_names)
    if len(common) == 0:
        raise ValueError("No shared cells between the three overlap objects")

    a1 = adata_rna[common].copy()
    a2 = adata_ac[common].copy()
    a3 = adata_me3[common].copy()

    if color_field not in a1.obs:
        raise KeyError(f"{color_field} not found in RNA overlap object")
    if color_field not in a2.obs:
        raise KeyError(f"{color_field} not found in H3K27ac overlap object")
    if color_field not in a3.obs:
        raise KeyError(f"{color_field} not found in H3K27me3 overlap object")

    umap1 = pd.DataFrame(a1.obsm["X_umap"], index=common, columns=["x1", "y1"])
    umap2 = pd.DataFrame(a2.obsm["X_umap"], index=common, columns=["x2", "y2"])
    umap3 = pd.DataFrame(a3.obsm["X_umap"], index=common, columns=["x3", "y3"])

    # Align AC and ME3 to RNA
    umap2[["x2", "y2"]] = procrustes_align_to_base(
        umap2[["x2", "y2"]].values,
        umap1[["x1", "y1"]].values
    )
    umap3[["x3", "y3"]] = procrustes_align_to_base(
        umap3[["x3", "y3"]].values,
        umap1[["x1", "y1"]].values
    )

    w = max(
        extent(umap1[["x1", "y1"]].values)[0],
        extent(umap2[["x2", "y2"]].values)[0],
        extent(umap3[["x3", "y3"]].values)[0],
    )
    h = max(
        extent(umap1[["x1", "y1"]].values)[1],
        extent(umap2[["x2", "y2"]].values)[1],
        extent(umap3[["x3", "y3"]].values)[1],
    )
    gap = gap_frac * w

    # Arrange AC - RNA - ME3
    umap2["x2"] += (umap1["x1"].min() - umap2["x2"].max()) - gap
    umap3["x3"] += (umap1["x1"].max() - umap3["x3"].min()) + gap
    umap1["y1"] += lift_rna_frac * h

    df = pd.concat([umap1, umap2, umap3], axis=1)
    df["color_group"] = a1.obs[color_field].astype(str).values

    cats_present = ordered_categories(a1.obs[color_field], desired_order)
    fallback = "#BDBDBD"

    rgba_lines = np.vstack([
        mcolors.to_rgba(palette_dict.get(ct, fallback), alpha=line_alpha)
        for ct in df["color_group"]
    ])

    fig, ax = plt.subplots(figsize=fig_size, dpi=fig_dpi)

    seg12 = np.stack([df[["x1", "y1"]].values, df[["x2", "y2"]].values], axis=1)
    seg13 = np.stack([df[["x1", "y1"]].values, df[["x3", "y3"]].values], axis=1)

    ax.add_collection(LineCollection(seg12, colors=rgba_lines, linewidths=line_width, zorder=0))
    ax.add_collection(LineCollection(seg13, colors=rgba_lines, linewidths=line_width, zorder=0))

    panel_info = [
        ("x2", "y2", "H3K27ac", a2),
        ("x1", "y1", "RNA", a1),
        ("x3", "y3", "H3K27me3", a3),
    ]

    for x, y, title, adx in panel_info:
        for ct in cats_present:
            sub = df[df["color_group"] == ct]
            if len(sub) == 0:
                continue
            ax.scatter(
                sub[x], sub[y],
                s=point_size,
                c=palette_dict.get(ct, fallback),
                edgecolors="none",
                alpha=0.98,
                label=ct if title == "RNA" else None,
                zorder=2 if title != "RNA" else 3,
                rasterized=False,
            )

        ax.text(
            df[x].median(),
            df[y].max() + 0.08 * h,
            title,
            ha="center",
            va="bottom",
            fontsize=15,
            weight="bold",
        )

    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.03)
    ax.axis("off")

    make_legend(ax, palette_dict, cats_present, title=color_field)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    plt.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close()


# ------------------------------------------------------------
# BUILD OVERLAP OBJECTS
# ------------------------------------------------------------
(
    rna_full, ac_full, me3_full,
    rna_set, ac_set, me3_set,
    rna_ov, ac_ov, me3_ov
) = build_threeway_overlap_objects(
    RNA_FILE,
    AC_FILE,
    ME3_FILE,
    recompute_rna=RECOMPUTE_RNA_UMAP,
)

print(f"[INFO] RNA cells: {len(rna_set)}")
print(f"[INFO] H3K27ac cells: {len(ac_set)}")
print(f"[INFO] H3K27me3 cells: {len(me3_set)}")
print(f"[INFO] RNA ∩ H3K27ac ∩ H3K27me3 cells: {rna_ov.n_obs}")

# Save exact overlap objects if useful
rna_ov.write("overlap_objects/MergedTonsils_RNA_overlap_AC_ME3_reUMAP.h5ad")
ac_ov.write("overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_ME3.h5ad")
me3_ov.write("overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_AC.h5ad")

# ------------------------------------------------------------
# VENN3
# ------------------------------------------------------------
plot_venn3(
    rna_set,
    ac_set,
    me3_set,
    output_prefix="figures/MergedTonsils_RNA_H3K27ac_H3K27me3_Venn3",
)

# ------------------------------------------------------------
# LINKED UMAPS
# ------------------------------------------------------------
plot_linked_umaps_three_modalities(
    adata_rna=rna_ov,
    adata_ac=ac_ov,
    adata_me3=me3_ov,
    color_field=CELLTYPE_KEY,
    palette_dict=CELL_TYPE_COLORS,
    desired_order=CELLTYPE_ORDER,
    output_prefix="figures/LinkedUMAPs_RNA_H3K27ac_H3K27me3_byCellType",
)

plot_linked_umaps_three_modalities(
    adata_rna=rna_ov,
    adata_ac=ac_ov,
    adata_me3=me3_ov,
    color_field=TONSIL_KEY,
    palette_dict=TONSIL_COLORS,
    desired_order=TONSIL_ORDER,
    output_prefix="figures/LinkedUMAPs_RNA_H3K27ac_H3K27me3_byTonsil",
)

print("[done] figures/MergedTonsils_RNA_H3K27ac_H3K27me3_Venn3.pdf")
print("[done] figures/LinkedUMAPs_RNA_H3K27ac_H3K27me3_byCellType.pdf")
print("[done] figures/LinkedUMAPs_RNA_H3K27ac_H3K27me3_byTonsil.pdf")