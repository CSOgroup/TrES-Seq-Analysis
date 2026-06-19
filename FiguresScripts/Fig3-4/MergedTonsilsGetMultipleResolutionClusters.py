#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from sklearn.metrics import normalized_mutual_info_score, silhouette_score
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ---------------------------------------------------
# INPUTS
# ---------------------------------------------------
INPUTS = {
    "RNA": {
        "path": "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad",
        "use_rep": "X_pca",
        "metric": "cosine",
        "k_min": 8,
        "k_max": 16,
        "highlight_k": 12,
    },
    "H3K27ac": {
        "path": "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_min": 3,
        "k_max": 11,
        "highlight_k": 8,
    },
    "H3K27me3": {
        "path": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_min": 3,
        "k_max": 10,
        "highlight_k": 6,
    },
    "H3K9me3": {
        "path": "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "use_rep": "X_spectral",
        "metric": "cosine",
        "k_min": 3,
        "k_max": 11,
        "highlight_k": 6,
    },
}

CELLTYPE_KEY = "leiden_merged_type"
DROP_CELLTYPES = {"NotInRNA"}

DESIRED_ORDER = [
    "B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive", "B_plasma",
    "DC", "FDC", "Macrophages",
    "T_CD4+", "T_CD4+_TfH", "T_CD4+_TfH_GC", "T_CD4+_naive", "T_CD8+",
]

# tighter resolution sweep
RESOLUTIONS = np.round(np.linspace(0.03, 1.20, 55), 3)

N_REPEATS = 3
RANDOM_STATE = 0
MIN_CLUSTER_SIZE = 10
SILHOUETTE_SUBSAMPLE = 100

OUTROOT = "Plots_overlap3_clusterSweep"
AUTOK_DIR = os.path.join(OUTROOT, "autoK")
UMAP_DIR = os.path.join(OUTROOT, "umaps")
JACCARD_DIR = os.path.join(OUTROOT, "jaccard")
CONDPROP_DIR = os.path.join(OUTROOT, "conditional_props")
CSV_DIR = os.path.join(OUTROOT, "csv")

for d in [OUTROOT, AUTOK_DIR, UMAP_DIR, JACCARD_DIR, CONDPROP_DIR, CSV_DIR]:
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


def leiden_key(resolution, repeat_idx):
    return f"_tmp_leiden_r{repeat_idx}_res{str(resolution).replace('.', '_')}"


def pairwise_nmi_values(label_matrix):
    vals = []
    n = len(label_matrix)
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(normalized_mutual_info_score(label_matrix[i], label_matrix[j]))
    return vals


def auto_res_candidates(
    adata,
    resolutions,
    use_rep,
    metric,
    n_repeats=10,
    random_state=0,
    min_cluster_size=10,
    sample_for_silhouette=10000,
):
    adata = adata.copy()
    sc.pp.neighbors(adata, use_rep=use_rep, metric=metric, n_neighbors=20)

    label_bank = {res: [] for res in resolutions}
    n_cells = adata.n_obs

    for r in range(n_repeats):
        seed = int(random_state + r)
        for res in resolutions:
            key = leiden_key(res, r)
            sc.tl.leiden(adata, resolution=float(res), key_added=key, random_state=seed)
            label_bank[res].append(adata.obs[key].to_numpy(dtype=str))

    X = adata.obsm[use_rep]
    if n_cells > sample_for_silhouette:
        rng = np.random.default_rng(random_state)
        idx_sil = np.sort(rng.choice(n_cells, size=sample_for_silhouette, replace=False))
        X_for_sil = X[idx_sil]
    else:
        idx_sil = slice(None)
        X_for_sil = X

    rows = []
    stab_rows = []

    for res in resolutions:
        labels0 = label_bank[res][0]
        _, counts = np.unique(labels0, return_counts=True)
        n_clusters = int(len(counts))
        min_size = int(counts.min()) if n_clusters > 0 else 0

        if n_clusters >= 2:
            labs_for_sil = labels0[idx_sil] if isinstance(idx_sil, np.ndarray) else labels0
            try:
                sil = float(silhouette_score(X_for_sil, labs_for_sil, metric="euclidean"))
            except Exception:
                sil = np.nan
        else:
            sil = np.nan

        pair_vals = pairwise_nmi_values(label_bank[res])
        mean_nmi = float(np.mean(pair_vals)) if len(pair_vals) else np.nan

        rows.append({
            "resolution": float(res),
            "n_clusters": n_clusters,
            "stability_nmi": mean_nmi,
            "silhouette": sil,
            "min_cluster_size": min_size,
        })

        stab_rows.extend([{"N_clusters": n_clusters, "Stability": v} for v in pair_vals])

    df = pd.DataFrame(rows).sort_values(["n_clusters", "resolution"]).reset_index(drop=True)
    stability_long = pd.DataFrame(stab_rows)

    df["rank_stability"] = df["stability_nmi"].rank(ascending=False, method="min")
    df["rank_silhouette"] = df["silhouette"].rank(ascending=False, method="min")
    penalty = (df["min_cluster_size"] < min_cluster_size).astype(float) * 2.0
    df["score"] = (df["rank_stability"] + df["rank_silhouette"]) / 2.0 + penalty

    return df, stability_long


def choose_best_per_k(results_df):
    if results_df.empty:
        return results_df.copy()
    best = (
        results_df
        .sort_values(["n_clusters", "score", "resolution"], ascending=[True, True, True])
        .groupby("n_clusters", as_index=False)
        .first()
        .sort_values("n_clusters")
        .reset_index(drop=True)
    )
    return best


def pick_resolution_for_K(results_df, K):
    dfK = results_df.loc[results_df["n_clusters"] == K].copy()
    if dfK.empty:
        raise ValueError(f"No rows with n_clusters == {K}")
    dfK = dfK.sort_values(["score", "resolution"], ascending=[True, True])
    return float(dfK.iloc[0]["resolution"])


def save_autok_plots(results_df, stability_long, outdir, prefix, highlight_k=None):
    os.makedirs(outdir, exist_ok=True)

    if not stability_long.empty:
        means = (
            stability_long
            .groupby("N_clusters", as_index=False)["Stability"]
            .mean()
            .rename(columns={"N_clusters": "K"})
        )

        fig, ax = plt.subplots(figsize=(5.0, 3.2))
        ax.plot(means["K"], means["Stability"], marker="o")
        ax.set_xlabel("N. clusters")
        ax.set_ylabel("Stability (NMI)")
        ax.set_title("Clustering stability")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        if highlight_k is not None and highlight_k in set(means["K"]):
            y = float(means.loc[means["K"] == highlight_k, "Stability"].iloc[0])
            ax.scatter([highlight_k], [y], s=150, color="red", zorder=5)

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{prefix}_stability_vs_K.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(outdir, f"{prefix}_stability_vs_K.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    sil = (
        results_df
        .dropna(subset=["silhouette"])
        .groupby("n_clusters", as_index=False)["silhouette"]
        .mean()
        .rename(columns={"n_clusters": "K"})
    )

    if not sil.empty:
        fig, ax = plt.subplots(figsize=(5.0, 3.2))
        ax.plot(sil["K"], sil["silhouette"], marker="o")
        ax.set_xlabel("N. clusters")
        ax.set_ylabel("Silhouette")
        ax.set_title("Silhouette vs. number of clusters")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        if highlight_k is not None and highlight_k in set(sil["K"]):
            y = float(sil.loc[sil["K"] == highlight_k, "silhouette"].iloc[0])
            ax.scatter([highlight_k], [y], s=150, color="red", zorder=5)

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{prefix}_silhouette_vs_K.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(outdir, f"{prefix}_silhouette_vs_K.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)


def ordered_celltypes(values, desired_order):
    present = list(pd.Index(values).unique())
    desired_present = [x for x in desired_order if x in present]
    remaining = [x for x in present if x not in desired_present]
    return desired_present + remaining


def compute_jaccard_matrix(celltypes, clusters):
    ct_levels = np.unique(celltypes)
    cl_levels = np.unique(clusters)

    mat = np.zeros((len(ct_levels), len(cl_levels)), dtype=float)

    for i, ct in enumerate(ct_levels):
        idx_ct = np.where(celltypes == ct)[0]
        n_ct = idx_ct.size
        for j, cl in enumerate(cl_levels):
            idx_cl = np.where(clusters == cl)[0]
            n_cl = idx_cl.size
            inter = np.intersect1d(idx_ct, idx_cl).size
            denom = n_ct + n_cl - inter
            mat[i, j] = inter / denom if denom > 0 else 0.0

    return mat, ct_levels, cl_levels


def reorder_heatmap_like_request(mat, ct_levels, cl_levels, desired_order):
    col_max = np.max(mat, axis=0)
    col_idx = np.argsort(-col_max)
    cl_sorted = cl_levels[col_idx]

    row_max = np.max(mat, axis=1)
    row_map = {x: i for i, x in enumerate(ct_levels)}

    desired_present = [x for x in desired_order if x in row_map]
    remaining = [x for x in ct_levels if x not in desired_present]
    remaining_sorted = sorted(remaining, key=lambda x: -row_max[row_map[x]])

    ct_sorted = np.array(desired_present + remaining_sorted, dtype=object)
    row_idx = np.array([row_map[x] for x in ct_sorted], dtype=int)

    mat = mat[row_idx][:, col_idx]
    return mat, ct_sorted, cl_sorted


def jaccard_types(adata, cluster_key, outpath, title):
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

    mat, ct_levels, cl_levels = compute_jaccard_matrix(
        df["celltype"].to_numpy(),
        df["cluster"].to_numpy()
    )
    mat, ct_sorted, cl_sorted = reorder_heatmap_like_request(mat, ct_levels, cl_levels, DESIRED_ORDER)

    ct_counts = [(df["celltype"] == x).sum() for x in ct_sorted]
    cl_counts = [(df["cluster"] == x).sum() for x in cl_sorted]

    x_labels = [f"{x} ({n})" for x, n in zip(cl_sorted, cl_counts)]
    y_labels = [f"{x} ({n})" for x, n in zip(ct_sorted, ct_counts)]

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(mat, cmap="coolwarm", vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(cl_sorted)))
    ax.set_yticks(np.arange(len(ct_sorted)))
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_yticklabels(y_labels)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8)

    ax.set_xlabel("Unbiased Leiden clusters")
    ax.set_ylabel("Transferred cell type")
    ax.set_title(title, pad=18)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


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
    adata = ensure_rep_and_umap(
        adata,
        use_rep=cfg["use_rep"],
        metric=cfg["metric"],
        n_neighbors=5,
        n_pcs=10,
    )

    if CELLTYPE_KEY not in adata.obs.columns:
        raise KeyError(f"{CELLTYPE_KEY} not found in {cfg['path']}")

    print(f"{modality}: n_cells = {adata.n_obs}")
    print(f"{modality}: using rep = {cfg['use_rep']}")

    res_df, stab_df = auto_res_candidates(
        adata=adata,
        resolutions=RESOLUTIONS,
        use_rep=cfg["use_rep"],
        metric=cfg["metric"],
        n_repeats=N_REPEATS,
        random_state=RANDOM_STATE,
        min_cluster_size=MIN_CLUSTER_SIZE,
        sample_for_silhouette=SILHOUETTE_SUBSAMPLE,
    )

    k_min = cfg["k_min"]
    k_max = cfg["k_max"]

    res_df = res_df[(res_df["n_clusters"] >= k_min) & (res_df["n_clusters"] <= k_max)].copy()
    if res_df.empty:
        raise ValueError(f"{modality}: no clusterings found in requested K range [{k_min}, {k_max}]")

    best_df = choose_best_per_k(res_df)
    best_csv = os.path.join(CSV_DIR, f"{modality}_best_resolution_per_K.csv")
    full_csv = os.path.join(CSV_DIR, f"{modality}_all_resolution_metrics.csv")
    best_df.to_csv(best_csv, index=False)
    res_df.to_csv(full_csv, index=False)

    save_autok_plots(
        results_df=res_df,
        stability_long=stab_df,
        outdir=AUTOK_DIR,
        prefix=modality,
        highlight_k=cfg.get("highlight_k", None),
    )

    print(f"{modality}: cluster sweep kept K = {sorted(best_df['n_clusters'].tolist())}")

    for _, row in best_df.iterrows():
        K = int(row["n_clusters"])
        res = float(row["resolution"])
        res_str = f"{res:.3f}".replace(".", "_")
        cluster_key = f"leiden_K{K}_res{res_str}"

        sc.tl.leiden(adata, resolution=res, key_added=cluster_key, random_state=RANDOM_STATE)

        # UMAP
        sc.pl.umap(
            adata,
            color=[cluster_key, CELLTYPE_KEY],
            ncols=2,
            show=False,
        )
        plt.savefig(
            os.path.join(UMAP_DIR, f"{modality}_K{K}_res{res_str}_UMAP.pdf"),
            bbox_inches="tight",
        )
        plt.close()

        # Jaccard
        jaccard_types(
            adata=adata,
            cluster_key=cluster_key,
            outpath=os.path.join(JACCARD_DIR, f"{modality}_K{K}_res{res_str}_Jaccard.pdf"),
            title=f"{modality}: Jaccard, K={K}, res={res:.3f}",
        )

        # Conditional proportions
        conditional_props_rowsum1(
            adata=adata,
            cluster_key=cluster_key,
            outpath=os.path.join(CONDPROP_DIR, f"{modality}_K{K}_res{res_str}_CondProp_rowsum1.pdf"),
            title=f"{modality}: P(cluster | cell type), K={K}, res={res:.3f}",
        )

        summary_rows.append({
            "modality": modality,
            "K": K,
            "resolution": res,
            "stability_nmi": row["stability_nmi"],
            "silhouette": row["silhouette"],
            "min_cluster_size": row["min_cluster_size"],
            "score": row["score"],
            "umap_pdf": os.path.join(UMAP_DIR, f"{modality}_K{K}_res{res_str}_UMAP.pdf"),
            "jaccard_pdf": os.path.join(JACCARD_DIR, f"{modality}_K{K}_res{res_str}_Jaccard.pdf"),
            "condprop_pdf": os.path.join(CONDPROP_DIR, f"{modality}_K{K}_res{res_str}_CondProp_rowsum1.pdf"),
        })

        print(f"[{modality}] done: K={K}, resolution={res:.3f}")

summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(CSV_DIR, "cluster_sweep_summary.csv")
summary_df.to_csv(summary_path, index=False)

print("\nDone.")
print(f"Summary written to: {summary_path}")