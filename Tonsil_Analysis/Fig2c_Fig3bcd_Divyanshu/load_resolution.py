"""Load preprocessed embeddings and recompute neighbors on the trimodal overlap."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import config


def _patch_anndata_raw():
    import anndata as ad

    if not hasattr(ad.AnnData, "Raw"):
        from anndata import Raw

        ad.AnnData.Raw = Raw  # type: ignore[attr-defined]


def ensure_celltype_column(adata) -> None:
    """Rename the stored cell-type column to ``CellType`` if needed."""
    src = config.CELLTYPE_SOURCE_KEY
    dst = config.CELLTYPE_KEY
    if dst in adata.obs.columns:
        return
    if src not in adata.obs.columns:
        raise KeyError(f"missing obs['{dst}'] (and no obs['{src}'] to rename)")
    adata.obs.rename(columns={src: dst}, inplace=True)


def three_way_cell_ids_tonsil1() -> list[str]:
    import scanpy as sc

    paths = {
        "RNA": config.tonsil1_rna_h5ad(),
        "H3K27ac": config.tonsil1_dna_h5ad("H3K27ac"),
        "H3K27me3": config.tonsil1_dna_h5ad("H3K27me3"),
    }
    ids = {}
    for mod, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        adata = sc.read_h5ad(path, backed="r")
        try:
            ids[mod] = set(map(str, adata.obs_names))
        finally:
            if getattr(adata, "file", None) is not None:
                adata.file.close()
    common = sorted(ids["RNA"] & ids["H3K27ac"] & ids["H3K27me3"])
    if len(common) != config.TONSIL1_N_CELLS:
        raise RuntimeError(
            f"Tonsil1 three-way overlap is {len(common)} cells; "
            f"expected {config.TONSIL1_N_CELLS}."
        )
    return common


def three_way_cell_ids_four_tonsils() -> list[str]:
    import scanpy as sc

    ids = {}
    for mod, path in config.MERGED_H5AD.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        adata = sc.read_h5ad(path, backed="r")
        try:
            ids[mod] = set(map(str, adata.obs_names))
        finally:
            if getattr(adata, "file", None) is not None:
                adata.file.close()
    common = sorted(ids["RNA"] & ids["H3K27ac"] & ids["H3K27me3"])
    if len(common) != config.FOUR_TONSILS_N_CELLS:
        raise RuntimeError(
            f"Four-tonsil three-way overlap is {len(common)} cells; "
            f"expected {config.FOUR_TONSILS_N_CELLS}."
        )
    return common


def _neighbors_params(adata, emb_key: str) -> dict:
    if config.NEIGHBORS_KEY not in adata.uns:
        raise KeyError(f"missing uns['{config.NEIGHBORS_KEY}']")
    params = dict(adata.uns[config.NEIGHBORS_KEY].get("params", {}))
    use_rep = params.get("use_rep")
    if use_rep is not None and str(use_rep) != emb_key:
        raise RuntimeError(
            f"neighbors use_rep={use_rep!r} but embedding is {emb_key!r}"
        )
    return {
        "n_neighbors": int(params.get("n_neighbors", 30)),
        "metric": str(params.get("metric", "euclidean")),
        "random_state": int(params.get("random_state", config.RANDOM_STATE)),
        "method": str(params.get("method", "umap")),
    }


def _recompute_neighbors(adata, emb_key: str, nn_params: dict) -> dict:
    import scanpy as sc

    for key in list(adata.obsp.keys()):
        del adata.obsp[key]
    for key in list(adata.uns.keys()):
        if str(key).startswith("neighbors"):
            del adata.uns[key]

    sc.pp.neighbors(
        adata,
        n_neighbors=nn_params["n_neighbors"],
        use_rep=emb_key,
        metric=nn_params["metric"],
        method=nn_params["method"],
        random_state=nn_params["random_state"],
        key_added=config.NEIGHBORS_KEY,
    )
    nk = adata.uns[config.NEIGHBORS_KEY]
    conn_key = nk.get("connectivities_key", config.CONNECTIVITIES_KEY)
    if conn_key != config.CONNECTIVITIES_KEY and conn_key in adata.obsp:
        adata.obsp[config.CONNECTIVITIES_KEY] = adata.obsp[conn_key]
        nk["connectivities_key"] = config.CONNECTIVITIES_KEY
    return {
        "n_neighbors": nn_params["n_neighbors"],
        "metric": nn_params["metric"],
        "use_rep": emb_key,
        "n_cells": int(adata.n_obs),
        "nnz": int(adata.obsp[config.CONNECTIVITIES_KEY].nnz),
    }


def _load_one(path: Path, cells: list[str], modality: str, note: str):
    import scanpy as sc

    adata = sc.read_h5ad(path)
    missing = [c for c in cells if c not in adata.obs_names]
    if missing:
        raise KeyError(
            f"{path.name}: {len(missing)} cells missing (e.g. {missing[0]})"
        )
    adata = adata[cells].copy()
    ensure_celltype_column(adata)

    emb_key = config.EMBEDDING_KEY[modality]
    if emb_key not in adata.obsm:
        raise KeyError(f"{path.name} missing obsm['{emb_key}']")

    nn_params = _neighbors_params(adata, emb_key)
    adata.X = None
    if adata.raw is not None:
        adata.raw = None
    for layer in list(adata.layers.keys()):
        del adata.layers[layer]
    keep_obsm = {emb_key, "X_umap"}
    for key in list(adata.obsm.keys()):
        if key not in keep_obsm:
            del adata.obsm[key]
    adata.obsm["X_embeddings"] = np.asarray(adata.obsm[emb_key], dtype=np.float64)

    nn_info = _recompute_neighbors(adata, emb_key, nn_params)
    adata.uns["resolution_neighbors"] = {**nn_info, "note": note}
    return adata, nn_info


def _bundle_from_paths(label: str, cells: list[str], paths: dict, expected: int) -> dict:
    import muon as mu

    _patch_anndata_raw()
    mods = {}
    neighbors_info = {}
    for mod in config.MODALITIES:
        adata, nn_info = _load_one(
            paths[mod],
            cells,
            mod,
            note=f"recomputed on {label} three-way overlap; embeddings unchanged",
        )
        mods[mod] = adata
        neighbors_info[mod] = nn_info

    mdata = mu.MuData(mods)
    if list(mdata.obs_names.astype(str)) != cells:
        mdata = mdata[cells].copy()

    if config.CELLTYPE_KEY not in mdata["RNA"].obs.columns:
        raise KeyError(f"RNA missing obs['{config.CELLTYPE_KEY}']")
    ct = mdata["RNA"].obs[config.CELLTYPE_KEY].astype(str)
    mdata.obs[config.CELLTYPE_KEY] = ct.values
    for mod in config.MODALITIES:
        mdata[mod].obs[config.CELLTYPE_KEY] = ct.values

    sample_counts = {}
    if config.SAMPLE_KEY in mdata["RNA"].obs.columns:
        sample = mdata["RNA"].obs[config.SAMPLE_KEY].astype(str)
        mdata.obs[config.SAMPLE_KEY] = sample.values
        for mod in config.MODALITIES:
            mdata[mod].obs[config.SAMPLE_KEY] = sample.values
        sample_counts = sample.value_counts().astype(int).to_dict()

    if int(mdata.n_obs) != expected:
        raise RuntimeError(f"{label}: n_cells={mdata.n_obs}, expected {expected}")

    return {
        "label": label,
        "mdata": mdata,
        "cells": cells,
        "n_cells": int(mdata.n_obs),
        "ct": ct.reset_index(drop=True),
        "neighbors_info": neighbors_info,
        "paths": {k: Path(v).name for k, v in paths.items()},
        "celltype_counts": ct.value_counts().astype(int).to_dict(),
        "sample_counts": sample_counts,
    }


def load_tonsil1() -> dict:
    paths = {
        "RNA": config.tonsil1_rna_h5ad(),
        "H3K27ac": config.tonsil1_dna_h5ad("H3K27ac"),
        "H3K27me3": config.tonsil1_dna_h5ad("H3K27me3"),
    }
    cells = three_way_cell_ids_tonsil1()
    return _bundle_from_paths("Tonsil1", cells, paths, config.TONSIL1_N_CELLS)


def load_four_tonsils() -> dict:
    paths = {mod: config.MERGED_H5AD[mod] for mod in config.MODALITIES}
    cells = three_way_cell_ids_four_tonsils()
    return _bundle_from_paths(
        "Tonsil1–4", cells, paths, config.FOUR_TONSILS_N_CELLS
    )


def banner(bundle: dict) -> None:
    print("=" * 72)
    print(f"Resolution sweep input — {bundle['label']} (preprocessed objects)")
    print(f"  n_cells={bundle['n_cells']}")
    if bundle.get("sample_counts"):
        print(f"  per tonsil: {bundle['sample_counts']}")
    print(f"  embeddings: {dict(config.EMBEDDING_KEY)} (unchanged)")
    print("  neighbors: recomputed on the three-way overlap")
    for mod, info in bundle.get("neighbors_info", {}).items():
        print(
            f"    {mod}: n_neighbors={info['n_neighbors']} "
            f"metric={info['metric']} nnz={info['nnz']}"
        )
    for mod, name in bundle["paths"].items():
        print(f"  {mod}: {name}")
    print(
        f"  cell types ({config.CELLTYPE_KEY}): "
        f"{len(bundle['celltype_counts'])} labels"
    )
    print("=" * 72)


def run_leiden(adata, resolution: float, key_added: str):
    import scanpy as sc

    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        directed=False,
        n_iterations=2,
        key_added=key_added,
        neighbors_key=config.NEIGHBORS_KEY,
        random_state=config.RANDOM_STATE,
    )
    return adata.obs[key_added]


def leiden_vs_celltype(bundle: dict):
    """Sweep Leiden resolutions and score each partition against cell types."""
    import pandas as pd

    from . import metrics

    mdata = bundle["mdata"]
    ct = bundle["ct"]
    rows = []
    labels_by_res = {}
    for res in config.LEIDEN_RES_SWEEP:
        res_r = round(float(res), 2)
        labs = {}
        for mod in config.MODALITIES:
            key = f"leiden_r{res_r}"
            labs[mod] = np.asarray(
                run_leiden(mdata[mod], resolution=float(res), key_added=key)
            ).astype(str)
            row = metrics.vs_celltype(labs[mod], ct)
            rows.append(
                {
                    "algorithm": "leiden",
                    "modality": mod,
                    "param": float(res_r),
                    **row,
                }
            )
        labels_by_res[res_r] = labs
        print(
            f"  r={res_r:.2f}: "
            + ", ".join(
                f"{mod} n={pd.Series(labs[mod]).nunique()} ARI="
                f"{next(r['ARI'] for r in rows if r['modality']==mod and r['param']==float(res_r)):.3f}"
                for mod in config.MODALITIES
            ),
            flush=True,
        )
    return pd.DataFrame(rows), labels_by_res


def plot_resolution_sweep(vs, title: str):
    """Display-only 2×3 resolution-sweep figure (no file write)."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    leiden = vs.query("algorithm == 'leiden'").copy()
    metric_panels = [
        ("n_clusters", "n clusters"),
        ("ARI", "ARI vs cell type"),
        ("NMI", "NMI vs cell type"),
        ("H_celltype_given_clusters", "H(cell type | clusters), bits"),
        ("H_clusters_given_celltype", "H(clusters | cell type), bits"),
        ("homogeneity_clusters_given_celltype", "homogeneity"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.2), sharex=True)
    for ax, (col, ylabel) in zip(axes.ravel(), metric_panels):
        for mod in config.MODALITIES:
            sub = leiden.loc[leiden["modality"] == mod].sort_values("param")
            ax.plot(
                sub["param"],
                sub[col],
                marker="o",
                ms=4,
                color=config.MOD_COLORS[mod],
                label=mod,
            )
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Leiden resolution")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(title, y=1.02)
    sns.despine()
    fig.tight_layout()
    # fig.savefig(...)  # disabled
    return fig, axes
