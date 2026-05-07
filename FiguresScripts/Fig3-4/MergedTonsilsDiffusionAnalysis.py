#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import gridspec

# --------------------------------------------------
# INPUT / SETTINGS
# --------------------------------------------------
RNA_PATH = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"
OUTDIR = "major_diffusion_out"
os.makedirs(OUTDIR, exist_ok=True)

GROUP_KEY = "leiden_merged_type"

# B-cell subset only, excluding plasma
B_PREFIX = "B_"
EXCLUDE_CONTAINS = ["plasma"]

# endpoints used to define the major diffusion component
END_A = "B_naive"
END_B = "B_mem"

# embedding / graph params before diffmap
N_TOP_GENES = 6000
N_PCS = 30
N_NEIGHBORS = 5
DIFFMAP_N_COMPS = 30
NEIGHBOR_METRIC = "cosine"

# output names
OUT_PREFIX = os.path.join(OUTDIR, "tonsil_allBcells")
DC_CSV = f"{OUT_PREFIX}.diffmap_DCs_per_cell.csv.gz"
MDC_CSV = f"{OUT_PREFIX}.mDC_pseudotime_per_cell.csv.gz"
PLOT_PDF = f"{OUT_PREFIX}.major_diffusion_component_celltypes.pdf"
PLOT_PNG = f"{OUT_PREFIX}.major_diffusion_component_celltypes.png"

CELL_TYPE_COLORS = {
    "B_Cycling": "#00BFC4",
    "B_GC_DZ": "#98DF8A",
    "B_GC_LZ": "#2CA02C",
    "B_mem": "#2B8C9E",
    "B_naive": "#1F77B4",
    "B_plasma": "#76C1FF",
}

# desired row order for the plot
PLOT_ORDER = [
    "B_naive",
    "B_GC_LZ",
    "B_Cycling",
    "B_GC_DZ",
    "B_mem",
]

PRETTY_LABELS = {
    "B_naive": "Naive",
    "B_GC_LZ": "GC-LZ",
    "B_Cycling": "GC-cycling",
    "B_GC_DZ": "GC-DZ",
    "B_mem": "Memory",
}


# --------------------------------------------------
# HELPERS
# --------------------------------------------------
def resolve_endpoint(label: str, labels: np.ndarray) -> str:
    """
    Resolve endpoint label against observed labels with:
      1) exact match
      2) case-insensitive exact match
      3) case-insensitive contains match
    Must resolve uniquely.
    """
    s = str(label)
    uniq = list(dict.fromkeys(labels.astype(str).tolist()))

    if s in uniq:
        return s

    low = s.lower()
    hits = [u for u in uniq if u.lower() == low]
    if len(hits) == 1:
        return hits[0]

    hits = [u for u in uniq if low in u.lower()]
    if len(hits) == 1:
        return hits[0]

    raise ValueError(
        f"Endpoint '{label}' not found uniquely. "
        f"Observed labels: {uniq}"
    )


def prepare_rna_for_diffmap(adata: sc.AnnData) -> sc.AnnData:
    """
    Use existing normalized layer if available, otherwise derive log1p from counts.
    """
    ad = adata.copy()

    if "log1p_norm" in ad.layers:
        ad.X = ad.layers["log1p_norm"].copy()
        hvg_layer = "counts" if "counts" in ad.layers else None
    elif "SCT_data" in ad.layers:
        ad.X = ad.layers["SCT_data"].copy()
        hvg_layer = "SCT_counts" if "SCT_counts" in ad.layers else None
    elif "counts" in ad.layers:
        ad.X = ad.layers["counts"].copy()
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        hvg_layer = "counts"
    else:
        hvg_layer = None

    hvg_kwargs = dict(
        flavor="seurat_v3",
        n_top_genes=N_TOP_GENES,
        subset=False,
    )
    if hvg_layer is not None and hvg_layer in ad.layers:
        hvg_kwargs["layer"] = hvg_layer

    sc.pp.highly_variable_genes(ad, **hvg_kwargs)
    sc.pp.scale(ad, zero_center=False, max_value=10)
    sc.pp.pca(ad, n_comps=max(N_PCS, 30), use_highly_variable=True)
    sc.pp.neighbors(ad, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, metric=NEIGHBOR_METRIC)
    sc.tl.diffmap(ad, n_comps=DIFFMAP_N_COMPS)

    return ad


def compute_major_diffusion_component(
    adata: sc.AnnData,
    group_key: str,
    end_a: str,
    end_b: str,
) -> tuple[sc.AnnData, int, float]:
    """
    Ignore DC0. Among DC1-DC9, choose the component with the largest absolute
    difference in mean value between endpoint groups. Flip sign if needed so that
    endpoint B has the larger mean. Store raw oriented component as 'mDC' and
    min-max scaled version as 'dc_pseudotime'.
    """
    ad = adata.copy()

    labels = ad.obs[group_key].astype(str).values
    end_a_res = resolve_endpoint(end_a, labels)
    end_b_res = resolve_endpoint(end_b, labels)

    Xdm = ad.obsm["X_diffmap"]
    if Xdm.shape[1] < 10:
        raise ValueError(f"Expected at least 10 diffusion components, got {Xdm.shape[1]}")

    mask_a = labels == end_a_res
    mask_b = labels == end_b_res

    diffs = []
    for dc in range(1, 10):  # ignore DC0
        va = Xdm[mask_a, dc]
        vb = Xdm[mask_b, dc]
        diffs.append(np.nanmean(vb) - np.nanmean(va))

    diffs = np.asarray(diffs, dtype=float)
    best_dc = int(np.argmax(np.abs(diffs)) + 1)  # back to actual DC index
    signed_diff = float(diffs[best_dc - 1])

    mdc = Xdm[:, best_dc].copy()
    if signed_diff < 0:
        mdc = -mdc
        signed_diff = -signed_diff

    dc_pseudotime = (mdc - np.nanmin(mdc)) / (np.nanmax(mdc) - np.nanmin(mdc) + 1e-12)

    ad.obs["mDC"] = mdc
    ad.obs["dc_pseudotime"] = dc_pseudotime

    print(
        f"[mDC] selected DC{best_dc} with max |Δmean| between "
        f"{end_a_res} and {end_b_res}; oriented so {end_b_res} > {end_a_res}; "
        f"Δmean={signed_diff:.6f}"
    )

    return ad, best_dc, signed_diff


def make_dc_export(adata: sc.AnnData, group_key: str, out_csv: str):
    Xdm = adata.obsm["X_diffmap"][:, :10]
    cols = [f"DC{i}" for i in range(10)]
    df = pd.DataFrame(Xdm, index=adata.obs_names, columns=cols)
    df[group_key] = adata.obs[group_key].astype(str).values
    df.to_csv(out_csv, compression="gzip")
    print(f"saved: {out_csv}")


def make_mdc_export(adata: sc.AnnData, group_key: str, out_csv: str):
    df = pd.DataFrame(index=adata.obs_names)
    df[group_key] = adata.obs[group_key].astype(str).values
    df["mDC"] = pd.to_numeric(adata.obs["mDC"], errors="coerce").values
    df["dc_pseudotime"] = pd.to_numeric(adata.obs["dc_pseudotime"], errors="coerce").values
    df.to_csv(out_csv, compression="gzip")
    print(f"saved: {out_csv}")


def plot_major_diffusion_celltype_track(
    adata: sc.AnnData,
    group_key: str,
    out_pdf: str,
    out_png: str,
):
    """
    Make the simple row-vs-major-diffusion plot:
      - cells ordered by dc_pseudotime
      - one row per cell type
      - colored vertical ticks where the cell belongs to that type
    """
    ad = adata.copy()
    ad = ad[np.argsort(ad.obs["dc_pseudotime"].to_numpy())].copy()

    present = pd.Index(ad.obs[group_key].astype(str).unique()).tolist()
    row_order = [x for x in PLOT_ORDER if x in present] + [x for x in present if x not in PLOT_ORDER]

    n_rows = len(row_order)
    n_cells = ad.n_obs
    x = np.arange(n_cells)

    fig = plt.figure(figsize=(12, max(3.5, 0.65 * n_rows + 1.6)), dpi=300)
    gs = gridspec.GridSpec(2, 1, height_ratios=[4.0, 0.9], hspace=0.05)

    ax = fig.add_subplot(gs[0])
    ax_arrow = fig.add_subplot(gs[1])

    ax.set_facecolor("white")

    labels = ad.obs[group_key].astype(str).to_numpy()

    for i, ct in enumerate(row_order):
        idx = x[labels == ct]
        if idx.size == 0:
            continue
        color = CELL_TYPE_COLORS.get(ct, "#BDBDBD")
        ax.vlines(idx, i - 0.42, i + 0.42, colors=color, linewidth=0.15)

    ax.set_xlim(-0.5, n_cells - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([PRETTY_LABELS.get(ct, ct) for ct in row_order], fontsize=16)

    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color("gray")

    ax.tick_params(axis="y", length=0)

    # bottom arrow
    ax_arrow.axis("off")
    ax_arrow.annotate(
        "",
        xy=(0.98, 0.55),
        xytext=(0.02, 0.55),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=1.6, color="black"),
    )
    ax_arrow.text(
        0.5,
        0.08,
        "Major Diffusion Component (mDC)",
        transform=ax_arrow.transAxes,
        ha="center",
        va="center",
        fontsize=18,
    )

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"saved: {out_pdf}")
    print(f"saved: {out_png}")


# --------------------------------------------------
# MAIN
# --------------------------------------------------
print(f"[load] {RNA_PATH}")
adata = sc.read_h5ad(RNA_PATH)

if GROUP_KEY not in adata.obs.columns:
    raise KeyError(f"{GROUP_KEY} not found in adata.obs")

group = adata.obs[GROUP_KEY].astype(str)

# keep B cells only, exclude plasma
keep = group.str.startswith(B_PREFIX)
for token in EXCLUDE_CONTAINS:
    keep &= ~group.str.contains(token, case=False, na=False)

adata = adata[keep].copy()
print(f"[subset] kept {adata.n_obs} B-cell non-plasma cells")

# compute diffmap
adata = prepare_rna_for_diffmap(adata)

# choose mDC exactly as requested
adata, best_dc, signed_diff = compute_major_diffusion_component(
    adata,
    group_key=GROUP_KEY,
    end_a=END_A,
    end_b=END_B,
)

# exports
make_dc_export(adata, GROUP_KEY, DC_CSV)
make_mdc_export(adata, GROUP_KEY, MDC_CSV)

# plot
plot_major_diffusion_celltype_track(
    adata,
    group_key=GROUP_KEY,
    out_pdf=PLOT_PDF,
    out_png=PLOT_PNG,
)

print("done.")