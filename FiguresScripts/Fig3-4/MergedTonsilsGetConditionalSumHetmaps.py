#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ---------------------------------------------------
# INPUTS
# ---------------------------------------------------
INPUTS = {
    "RNA": {
        "path": "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad",
        "use_rep": "X_pca",
        "metric": "cosine",
        "k_values": [9, 10,11, 12,13, 14,15, 16],
        # optional manual override:
        # "manual_resolutions": {8: 0.22, 10: 0.31}
    },
    "H3K27ac": {
        "path": "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_values": [3,4, 6,7, 8,9, 10],
    },
    "H3K27me3": {
        "path": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_values": [3,4, 6,7, 8,10],
    },
    "H3K9me3": {
        "path": "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_values": [3, 4, 5, 7,8,10],
    },
}

INPUTS = {
    "H3K27ac": {
        "path": "DNA_merged/MergedTonsils_H3K27ac_merged_processed_barcodeRewritten_withRNAlabels.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_values": [3,4, 6,7, 8,9, 10],
    },
}

PREVIOUS_CSV_DIR = "Plots_overlap3_clusterSweep/csv"

CELLTYPE_KEY = "leiden_merged_type"
DROP_CELLTYPES = {"NotInRNA"}

DESIRED_ORDER = [
    "B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive", "B_plasma",
    "DC", "FDC", "Macrophages",
    "T_CD4+", "T_CD4+_TfH", "T_CD4+_TfH_GC", "T_CD4+_naive", "T_CD8+",
]

RANDOM_STATE = 0
NN=50
NP=15
OUTROOT = "Plots_overlap3_fixedClusters_multiK"
UMAP_DIR = os.path.join(OUTROOT, "umaps")
CONDPROP_DIR = os.path.join(OUTROOT, "conditional_props")
CSV_DIR = os.path.join(OUTROOT, "csv")

for d in [OUTROOT, UMAP_DIR, CONDPROP_DIR, CSV_DIR]:
    os.makedirs(d, exist_ok=True)


# ---------------------------------------------------
# HELPERS
# ---------------------------------------------------
def ensure_rep_and_umap(adata, use_rep, metric="euclidean", n_neighbors=20, n_pcs=30):
    adata = adata.copy()

    if use_rep == "X_pca" and "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=n_pcs)

    if use_rep not in adata.obsm:
        raise ValueError(f"{use_rep} not found in adata.obsm")

    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        sc.pp.neighbors(adata, use_rep=use_rep, metric=metric, n_neighbors=n_neighbors)

    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    return adata


def resolve_resolution(modality, cfg, K, previous_csv_dir):
    manual = cfg.get("manual_resolutions", {})
    if K in manual:
        return float(manual[K])

    csv_path = os.path.join(previous_csv_dir, f"{modality}_best_resolution_per_K.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"{modality}: no manual resolution for K={K}, and lookup CSV not found:\n{csv_path}"
        )

    df = pd.read_csv(csv_path)
    if "n_clusters" not in df.columns or "resolution" not in df.columns:
        raise ValueError(f"{csv_path} must contain 'n_clusters' and 'resolution' columns")

    dfk = df.loc[df["n_clusters"] == K].copy()
    if dfk.empty:
        raise ValueError(f"{modality}: K={K} not found in {csv_path}")

    sort_cols = [c for c in ["score", "resolution"] if c in dfk.columns]
    if sort_cols:
        dfk = dfk.sort_values(sort_cols, ascending=[True] * len(sort_cols))

    return float(dfk.iloc[0]["resolution"])


def conditional_props_rowsum1(adata, cluster_key, outpath, title):
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"{cluster_key} not found in adata.obs")
    if CELLTYPE_KEY not in adata.obs.columns:
        raise KeyError(f"{CELLTYPE_KEY} not found in adata.obs")

    df = pd.DataFrame({
        "cluster": adata.obs[cluster_key].astype(str).values,
        "celltype": adata.obs[CELLTYPE_KEY].astype(str).values,
    }, index=adata.obs_names)

    df = df[~df["celltype"].isin(DROP_CELLTYPES)].copy()
    if df.empty:
        print(f"[skip] {title}: no cells left after dropping {DROP_CELLTYPES}")
        return

    celltypes = df["celltype"].to_numpy()
    clusters = df["cluster"].to_numpy()

    ct_levels = np.unique(celltypes)
    cl_levels = np.unique(clusters)

    count_mat = np.zeros((len(ct_levels), len(cl_levels)), dtype=float)
    jaccard_mat = np.zeros_like(count_mat)

    idx_by_ct = {ct: np.where(celltypes == ct)[0] for ct in ct_levels}
    idx_by_cl = {cl: np.where(clusters == cl)[0] for cl in cl_levels}
    n_ct_vec = np.array([idx_by_ct[x].size for x in ct_levels], dtype=float)
    n_cl_vec = np.array([idx_by_cl[x].size for x in cl_levels], dtype=float)

    for i, ct in enumerate(ct_levels):
        idx_ct = idx_by_ct[ct]
        n_ct = n_ct_vec[i]
        for j, cl in enumerate(cl_levels):
            idx_cl = idx_by_cl[cl]
            n_cl = n_cl_vec[j]
            inter = np.intersect1d(idx_ct, idx_cl).size
            count_mat[i, j] = inter
            denom = n_ct + n_cl - inter
            jaccard_mat[i, j] = inter / denom if denom > 0 else 0.0

    with np.errstate(invalid="ignore", divide="ignore"):
        prop_mat = (count_mat.T / np.where(n_ct_vec > 0, n_ct_vec, 1)).T
    prop_mat[np.isnan(prop_mat)] = 0.0

    col_max = np.max(jaccard_mat, axis=0)
    col_idx = np.argsort(-col_max)
    cl_sorted = cl_levels[col_idx]

    row_max = np.max(jaccard_mat, axis=1)
    row_map = {x: i for i, x in enumerate(ct_levels)}
    desired_present = [x for x in DESIRED_ORDER if x in row_map]
    remaining = [x for x in ct_levels if x not in desired_present]
    remaining_sorted = sorted(remaining, key=lambda x: -row_max[row_map[x]])
    ct_sorted = np.array(desired_present + remaining_sorted, dtype=object)
    row_idx = np.array([row_map[x] for x in ct_sorted], dtype=int)

    prop_mat = prop_mat[row_idx][:, col_idx]

    ct_counts = [(df["celltype"] == x).sum() for x in ct_sorted]
    cl_counts = [(df["cluster"] == x).sum() for x in cl_sorted]

    x_labels = [f"{x} ({n})" for x, n in zip(cl_sorted, cl_counts)]
    y_labels = [f"{x} ({n})" for x, n in zip(ct_sorted, ct_counts)]

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(prop_mat, cmap="coolwarm", vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(cl_sorted)))
    ax.set_yticks(np.arange(len(ct_sorted)))
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_yticklabels(y_labels)

    for i in range(prop_mat.shape[0]):
        for j in range(prop_mat.shape[1]):
            ax.text(j, i, f"{prop_mat[i, j]:.2f}", ha="center", va="center", fontsize=8)

    ax.set_xlabel("Unbiased Leiden clusters")
    ax.set_ylabel("Transferred cell type")
    ax.set_title(title, pad=18)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)

    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------
# MAIN
# ---------------------------------------------------
summary_rows = []

for modality, cfg in INPUTS.items():
    print("\n" + "=" * 90)
    print(f"Processing {modality}")
    print("=" * 90)

    adata = sc.read_h5ad(cfg["path"])
    adata = ensure_rep_and_umap(adata,use_rep=cfg["use_rep"],metric=cfg["metric"],n_neighbors=NN,n_pcs=NP,)

    if CELLTYPE_KEY not in adata.obs.columns:
        raise KeyError(f"{CELLTYPE_KEY} not found in {cfg['path']}")

    print(f"{modality}: n_cells = {adata.n_obs}")
    print(f"{modality}: using rep = {cfg['use_rep']}")
    print(f"{modality}: K values = {cfg['k_values']}")

    for K in cfg["k_values"]:
        resolution = resolve_resolution(modality, cfg, K, PREVIOUS_CSV_DIR)
        res_str = f"{resolution:.3f}".replace(".", "_")
        cluster_key = f"leiden_K{K}_res{res_str}"

        print(f"{modality}: running K={K}, resolution={resolution:.3f}")

        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=cluster_key,
            random_state=RANDOM_STATE,
        )

        observed_k = adata.obs[cluster_key].astype(str).nunique()
        print(f"{modality}: observed clusters for requested K={K} -> {observed_k}")

        umap_pdf = os.path.join(UMAP_DIR, f"{modality}_K{K}_res{res_str}_UMAP.pdf")
        sc.pl.umap(
            adata,
            color=[cluster_key, CELLTYPE_KEY],
            ncols=2,
            show=False,
        )
        plt.savefig(umap_pdf, bbox_inches="tight")
        plt.close()

        condprop_pdf = os.path.join(
            CONDPROP_DIR,
            f"{modality}_K{K}_res{res_str}_CondProp_rowsum1.pdf"
        )
        conditional_props_rowsum1(
            adata=adata,
            cluster_key=cluster_key,
            outpath=condprop_pdf,
            title=f"{modality}: P(cluster | cell type), K={K}, res={resolution:.3f}",
        )

        cluster_sizes = (
            adata.obs[cluster_key]
            .astype(str)
            .value_counts()
            .sort_index()
            .rename_axis("cluster")
            .reset_index(name="n_cells")
        )
        cluster_sizes_csv = os.path.join(
            CSV_DIR,
            f"{modality}_K{K}_res{res_str}_cluster_sizes.csv"
        )
        cluster_sizes.to_csv(cluster_sizes_csv, index=False)

        summary_rows.append({
            "modality": modality,
            "path": cfg["path"],
            "use_rep": cfg["use_rep"],
            "metric": cfg["metric"],
            "target_k": K,
            "resolution": resolution,
            "observed_k": observed_k,
            "cluster_key": cluster_key,
            "umap_pdf": umap_pdf,
            "condprop_pdf": condprop_pdf,
            "cluster_sizes_csv": cluster_sizes_csv,
        })

        print(f"[{modality}] done: K={K}")

summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(CSV_DIR, "fixed_cluster_summary_multiK.csv")
summary_df.to_csv(summary_path, index=False)

print("\nDone.")
print(f"Summary written to: {summary_path}")