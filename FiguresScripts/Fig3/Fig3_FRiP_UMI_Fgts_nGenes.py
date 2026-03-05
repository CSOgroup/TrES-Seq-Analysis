#!/usr/bin/env python3
import os
import re
import snapatac2 as snap
import scanpy as sc
from pathlib import Path

# --- I/O root (change if needed) ---
WORKDIR = Path(".")  # set to your folder if not the CWD

# Expected inputs (auto-discovered); alternatively, list them explicitly
H5AD_GLOB = "*_*_CountFilteredShared.h5ad"

# --- helpers ---
H5AD_RE = re.compile(r"^(?P<sample>[A-Za-z0-9]+)_(?P<mark>ac|me3)_CountFilteredShared.h5ad")

MACS_DIR = WORKDIR / "macs3_out"

def peak_for(sample: str, mark: str) -> Path:
    """
    Return the matching peak path for the sample/mark in macs3_out.

    - H3K27ac → narrowPeak
    - H3K27me3 → broadPeak
    """
    if mark == "ac":
        mark_str  = "H3K27ac"
        peak_type = "narrowPeak"
    elif mark == "me3":
        mark_str  = "H3K27me3"
        peak_type = "broadPeak"
    else:
        raise ValueError(f"Unknown mark: {mark!r}")

    return MACS_DIR / f"{sample}_{mark_str}_FiltShared_peaks.{peak_type}.blacklistFiltered"


def out_name(h5ad_path: Path) -> Path:
    return h5ad_path.with_name(h5ad_path.stem + "_withFRiP.h5ad")

# --- run ---
h5ads = sorted(WORKDIR.glob(H5AD_GLOB))
if not h5ads:
    raise SystemExit(f"No inputs matched {H5AD_GLOB} in {WORKDIR.resolve()}")

for h5 in h5ads:
    m = H5AD_RE.match(h5.name)
    if not m:
        print(f"Skipping (name pattern not recognized): {h5.name}")
        continue

    sample = m.group("sample")
    mark = m.group("mark")
    peaks = peak_for(sample, mark)

    if not peaks.exists():
        print(f"!! Peak file missing for {sample} {mark}: {peaks}")
        continue

    print(f"→ Calculating FRiP for {sample} {mark}")
    print(f"   • AnnData: {h5}")
    print(f"   • Peaks:   {peaks}")

    adata = sc.read_h5ad(h5)
    snap.metrics.frip(
        adata,
        regions={"FRiP": str(peaks)},
        normalized=True,
        inplace=True,
    )

    out_h5 = out_name(h5)
    adata.write(out_h5)
    print(f"✅ Saved: {out_h5}")

print("🎉 All FRiP calculations completed.")


#!/usr/bin/env python3
# frip_celllines_violins.py
# Plot FRiP per cell line, both marks in the SAME plot (nice colors).

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# --- Colors (match your previous conventions) ---
COLORS_BY_MARK = {
    'H3K27ac': '#D55E00',
    'H3K27me3': '#5E3D99',
}

# --- File name pattern: Sc_CMD_<CELL>_(ac|me3)_CountFiltered[...].h5ad ---
H5AD_RE = re.compile(
    r"^(?P<cell>[A-Za-z0-9]+)_(?P<mark>ac|me3)_CountFilteredShared_withFRiP.h5ad$"
)

CELL_ORDER_DEFAULT = ["A20", "GM12878", "JJN2", "Karpas422", "WSU"]


def pick_frip_col(obs: pd.DataFrame) -> str:
    """Find a FRiP column in .obs, robust to naming."""
    candidates = ['FRiP', 'frip', 'FRIP', 'FRiP_score', 'frip_score']
    for c in candidates:
        if c in obs.columns:
            return c
    lower = {c.lower(): c for c in obs.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    raise KeyError(f"Could not find a FRiP column in .obs. "
                   f"Saw columns like: {list(obs.columns)[:20]}...")


def stack_one(h5: Path) -> pd.DataFrame:
    m = H5AD_RE.match(h5.name)
    if not m:
        raise ValueError(f"Unrecognized name: {h5.name}")
    cell = m.group("cell")
    mark = "H3K27ac" if m.group("mark") == "ac" else "H3K27me3"

    ad = sc.read_h5ad(h5)
    frip_col = pick_frip_col(ad.obs)
    df = ad.obs[[frip_col]].copy()
    df.rename(columns={frip_col: "FRiP"}, inplace=True)

    # Ensure numeric and clip to [0, 1]
    df["FRiP"] = pd.to_numeric(df["FRiP"], errors="coerce").clip(0, 1)
    df["FRiP_pct"] = df["FRiP"] * 100.0
    df["CellLine"] = cell
    df["Mark"] = mark
    return df


def gather_inputs(pattern: str) -> list[Path]:
    paths = sorted(Path(".").glob(pattern))
    if not paths:
        raise SystemExit(f"No files matched: {pattern}")
    return paths


def plot_combined(df: pd.DataFrame, out_prefix: str):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import os

    # Order on x-axis
    present_cells = [c for c in CELL_ORDER_DEFAULT if c in df["CellLine"].unique()]
    if not present_cells:
        present_cells = sorted(df["CellLine"].unique())

    hue_order = ["H3K27ac", "H3K27me3"]
    color_map = {m: COLORS_BY_MARK[m] for m in hue_order}

    # Precompute positions: one pair per cell line
    n = len(present_cells)
    base = np.arange(n)
    offset = 0.20                # horizontal shift for the two marks within a cell line
    pos_for = {
        "H3K27ac": base - offset,
        "H3K27me3": base + offset,
    }

    # Assemble data lists for deterministic plotting order
    data_lists = []
    positions = []
    colors = []
    for cell in present_cells:
        for mark in hue_order:
            arr = df.loc[(df["CellLine"] == cell) & (df["Mark"] == mark), "FRiP_pct"].dropna().values
            if arr.size == 0:
                continue  # gracefully handle missing mark/cell combos
            data_lists.append(arr)
            # position index for this cell/mark (same order as data_lists)
            idx = np.where(present_cells == np.array(cell))[0][0] if isinstance(present_cells, np.ndarray) else present_cells.index(cell)
            positions.append(pos_for[mark][idx])
            colors.append(color_map[mark])

    # Plot
    os.makedirs("./figures", exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))

    # Violins
    v = ax.violinplot(
        data_lists,
        positions=positions,
        widths=0.38,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body, col in zip(v["bodies"], colors):
        body.set_facecolor(col)
        body.set_edgecolor("black")
        body.set_linewidth(0.6)
        body.set_alpha(0.85)

    # Boxes (transparent so violins show through)
    bp = ax.boxplot(
        data_lists,
        positions=positions,
        widths=0.20,
        patch_artist=True,          # needed so boxes are Patch objects
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.8),
    )

    # Make boxes transparent + draw outlines on top
    for box in bp["boxes"]:
        box.set_facecolor("none")   # transparent
        box.set_edgecolor("black")
        box.set_linewidth(1.6)
        box.set_zorder(3)

    for whisk in bp["whiskers"]:
        whisk.set_color("black")
        whisk.set_linewidth(1.4)
        whisk.set_zorder(3)

    for cap in bp["caps"]:
        cap.set_color("black")
        cap.set_linewidth(1.4)
        cap.set_zorder(3)

    for med in bp["medians"]:
        med.set_zorder(4)


    # X axis ticks at the center of each pair
    ax.set_xticks(base)
    ax.set_xticklabels(present_cells)

    # Legend
    legend_handles = [Patch(facecolor=color_map[m], edgecolor="black", label=m) for m in hue_order]
    ax.legend(legend_handles, hue_order, frameon=False, loc="upper right")

    ax.set_ylabel("FRiP (%)")
    ax.set_xlabel("")
    ax.set_ylim(0, 100)
    ax.set_title("FRiP per cell line (both marks)")
    ax.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    pdf = f"./figures/{out_prefix}_FRiP_CellLines_Violin.pdf"
    png = f"./figures/{out_prefix}_FRiP_CellLines_Violin.png"
    fig.savefig(pdf)
    #fig.savefig(png, dpi=300)
    plt.close(fig)
    print(pdf)




def main():
    ap = argparse.ArgumentParser(description="Plot FRiP per cell line with both marks in one plot.")
    ap.add_argument(
        "--inputs",
        default="*_*_CountFiltered*.h5ad",
        help="Glob for input .h5ad files (default: %(default)s). Matches both ac & me3."
    )
    ap.add_argument(
        "--out",
        default="FRiP_CellLines",
        help="Output prefix (default: %(default)s)"
    )
    args = ap.parse_args()

    files = gather_inputs(args.inputs)

    dfs = []
    for h5 in files:
        try:
            dfs.append(stack_one(h5))
        except Exception as e:
            print(f"[skip] {h5.name}: {e}")

    if not dfs:
        raise SystemExit("No usable inputs contained a FRiP column.")
    df = pd.concat(dfs, ignore_index=True)

    # Keep only the two expected marks (if others slipped in)
    df = df[df["Mark"].isin(["H3K27ac", "H3K27me3"])]

    plot_combined(df, args.out)


if __name__ == "__main__":
    main()

python frip_cl.py




#!/usr/bin/env python3
# rna_celllines_violins.py
# Plot RNA UMIs & n_genes per cell line from Sc_CMD_*_rna_CountFilteredShared.h5ad
# Now annotates each cell line with n cells and the median value.

import argparse
import os
import re
from pathlib import Path


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# ---- file naming ----
H5AD_RE = re.compile(r"^(?P<cell>[A-Za-z0-9]+)_rna_CountFilteredShared.*\.h5ad$")
CELL_ORDER_DEFAULT = ["A20", "GM12878", "JJN2", "Karpas422", "WSU"]

# Accept common obs column names
ALIAS = {
    "UMIs": ["UMIs", "total_counts", "n_counts", "counts", "umi", "umis"],
    "n_genes": ["n_genes", "n_genes_by_counts", "genes_detected", "n_features", "n_features_RNA"],
}

def pick_col(obs: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in obs.columns:
            return c
    lower = {c.lower(): c for c in obs.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    raise KeyError(f"None of {candidates} found in .obs. Saw: {list(obs.columns)[:20]} ...")

def stack_one(h5: Path) -> pd.DataFrame:
    m = H5AD_RE.match(h5.name)
    if not m:
        raise ValueError(f"Unrecognized filename: {h5.name}")
    cell = m.group("cell")

    ad = sc.read_h5ad(h5)
    umi_col = pick_col(ad.obs, ALIAS["UMIs"])
    ng_col  = pick_col(ad.obs, ALIAS["n_genes"])

    df = ad.obs[[umi_col, ng_col]].copy()
    df[umi_col] = pd.to_numeric(df[umi_col], errors="coerce")
    df[ng_col]  = pd.to_numeric(df[ng_col],  errors="coerce")

    # raw + log for plotting
    df["UMIs_raw"]    = df[umi_col].values
    df["n_genes_raw"] = df[ng_col].values
    df["UMIs_log"]    = np.log10(df["UMIs_raw"] + 1.0)
    df["n_genes_log"] = np.log10(df["n_genes_raw"] + 1.0)
    df["CellLine"]    = cell

    return df[["CellLine", "UMIs_raw", "n_genes_raw", "UMIs_log", "n_genes_log"]]

def gather_inputs(globpat: str) -> list[Path]:
    files = sorted(Path(".").glob(globpat))
    if not files:
        raise SystemExit(f"No files matched: {globpat}")
    return files

def fmt_log_ticks(ax):
    """
    Use nice round ticks for log10(value + 1) violins without forcing the axis to start at 0.

    - Axis is in log10(value + 1) space.
    - We keep the current y-limits (matplotlib defaults).
    - Ticks use round numbers: 1, 2, 3, 5 × 10^k within the data range (so e.g. 2,000; 3,000; 5,000).
    """
    ymin, ymax = ax.get_ylim()

    # Convert current limits back to raw scale (approximate)
    raw_min = max(0.0, 10**ymin - 1.0)
    raw_max = max(0.0, 10**ymax - 1.0)

    # If upper bound is non-positive, nothing sensible to do -> keep defaults, just relabel
    if raw_max <= 0:
        yticks = ax.get_yticks()
        labels = [f"{int(max(0, round(10**y - 1))):,}" for y in yticks]
        ax.set_yticklabels(labels)
        return

    # Candidate "nice" values: 1/2/3/5 * 10^k within the data range
    raw_min_eff = max(raw_min, 1.0)
    exp_min = int(np.floor(np.log10(raw_min_eff)))
    exp_max = int(np.ceil(np.log10(raw_max)))

    nice_vals = []
    for exp in range(exp_min - 1, exp_max + 2):
        for m in (1, 2, 3, 5):
            v = m * (10**exp)
            if raw_min <= v <= raw_max:
                nice_vals.append(v)

    nice_vals = sorted(set(nice_vals))

    # Fallback: if range is too tight, just relabel existing ticks nicely
    if len(nice_vals) < 2:
        yticks = ax.get_yticks()
        labels = [f"{int(max(0, round(10**y - 1))):,}" for y in yticks]
        ax.set_yticklabels(labels)
        return

    # Convert to axis positions in log10(value + 1)
    ticks = [np.log10(v + 1.0) for v in nice_vals]
    labels = [f"{int(v):,}" for v in nice_vals]

    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)



def _plot_single_metric(df: pd.DataFrame, value_col: str, raw_col: str,
                        ylabel: str, title: str, out_prefix: str, suffix: str):
    # order x-axis
    present_cells = [c for c in CELL_ORDER_DEFAULT if c in df["CellLine"].unique()]
    if not present_cells:
        present_cells = sorted(df["CellLine"].unique())

    # data lists per cell (deterministic order)
    data_log_lists = [df.loc[df["CellLine"] == c, value_col].dropna().values for c in present_cells]
    data_raw_lists = [df.loc[df["CellLine"] == c, raw_col].dropna().values for c in present_cells]
    base = np.arange(len(present_cells))

    # colors per cell line (tab10 subset)
    import matplotlib as mpl
    cmap = mpl.cm.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(len(present_cells))]

    os.makedirs("./figures", exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))

    # violins
    v = ax.violinplot(
        data_log_lists,
        positions=base,
        widths=0.7,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body, col in zip(v["bodies"], colors):
        body.set_facecolor(col)
        body.set_edgecolor("black")
        body.set_linewidth(0.6)
        body.set_alpha(0.85)

    # transparent box overlays (aligned)
    bp = ax.boxplot(
        data_log_lists,
        positions=base,
        widths=0.25,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.6),
    )
    for box in bp["boxes"]:
        box.set_facecolor("none")
        box.set_edgecolor("black")
        box.set_linewidth(1.4)
        box.set_zorder(3)
    for whisk in bp["whiskers"]:
        whisk.set_color("black"); whisk.set_linewidth(1.2); whisk.set_zorder(3)
    for cap in bp["caps"]:
        cap.set_color("black"); cap.set_linewidth(1.2); cap.set_zorder(3)
    for med in bp["medians"]:
        med.set_zorder(4)

    # axes + ticks back to raw counts
    ax.set_xticks(base)
    ax.set_xticklabels(present_cells)
    fmt_log_ticks(ax)


    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)

    # ---- Annotations per cell: n and median (RAW units) ----
    ymin, ymax = ax.get_ylim()
    y_annot = ymax - 0.06 * (ymax - ymin)  # near top; same for all groups
    for i, (cell, raw_vals) in enumerate(zip(present_cells, data_raw_lists)):
        if raw_vals.size == 0:
            continue
        n = raw_vals.size
        med = float(np.median(raw_vals))
        label = f"n={n:,} • median={int(round(med)):,}"
        ax.text(base[i], y_annot, label, ha="center", va="top",
                fontsize=9.5, rotation=0, clip_on=False)

    plt.tight_layout()
    pdf = f"./figures/{out_prefix}_{suffix}.pdf"
    png = f"./figures/{out_prefix}_{suffix}.png"
    fig.savefig(pdf)
    #fig.savefig(png, dpi=300)
    plt.close(fig)
    print(pdf); print(png)

def main():
    ap = argparse.ArgumentParser(description="RNA QC violins per cell line (UMIs & n_genes) with n and median annotations.")
    ap.add_argument("--inputs", default="*_rna_CountFilteredShared*.h5ad",
                    help="Glob for input .h5ad (default: %(default)s)")
    ap.add_argument("--out", default="RNA_CellLines", help="Output prefix")
    args = ap.parse_args()

    files = gather_inputs(args.inputs)

    dfs = []
    for h5 in files:
        try:
            dfs.append(stack_one(h5))
        except Exception as e:
            print(f"[skip] {h5.name}: {e}")
    if not dfs:
        raise SystemExit("No usable inputs.")

    df = pd.concat(dfs, ignore_index=True)

    _plot_single_metric(
        df, value_col="UMIs_log", raw_col="UMIs_raw",
        ylabel="RNA UMI Count", title="RNA – UMI Counts (per cell)",
        out_prefix=args.out, suffix="RNA_UMIs_Violin"
    )
    _plot_single_metric(
        df, value_col="n_genes_log", raw_col="n_genes_raw",
        ylabel="Detected Genes per Cell", title="RNA – Detected Genes (per cell)",
        out_prefix=args.out, suffix="RNA_Genes_Violin"
    )

if __name__ == "__main__":
    main()


python lp.py



#!/usr/bin/env python3
# dna_tsse_fragments_celllines.py
# TSSE & Fragment-count violins per cell line (colors = cell lines), split by mark (ac/me3).
# Requires .obs columns: 'tsse', 'n_fragment'

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---------- Config ----------
H5AD_RE = re.compile(r"^(?P<cell>[A-Za-z0-9]+)_(?P<mark>ac|me3)_CountFilteredShared.h5ad$")
CELL_ORDER_DEFAULT = ["A20", "GM12878", "JJN2", "Karpas422", "WSU"]

# Fixed, RNA-style cell colors from tab10 in that order:
_cmap = mpl.cm.get_cmap("tab10")
CELL_COLORS = {cell: _cmap(i % 10) for i, cell in enumerate(CELL_ORDER_DEFAULT)}


def fmt_log_ticks(ax):
    """
    Use nice round ticks for log10(value + 1) violins without forcing the axis to start at 0.

    - Axis is in log10(value + 1) space.
    - We keep the current y-limits (matplotlib defaults).
    - Ticks use round numbers: 1, 2, 3, 5 × 10^k within the data range
      (so e.g. 2,000; 3,000; 5,000).
    """
    ymin, ymax = ax.get_ylim()

    # Convert current limits back to raw scale (approximate)
    raw_min = max(0.0, 10**ymin - 1.0)
    raw_max = max(0.0, 10**ymax - 1.0)

    # If upper bound is non-positive, nothing sensible to do -> keep defaults, just relabel
    if raw_max <= 0:
        yticks = ax.get_yticks()
        labels = [f"{int(max(0, round(10**y - 1))):,}" for y in yticks]
        ax.set_yticklabels(labels)
        return

    # Candidate "nice" values: 1/2/3/5 * 10^k within the data range
    raw_min_eff = max(raw_min, 1.0)
    exp_min = int(np.floor(np.log10(raw_min_eff)))
    exp_max = int(np.ceil(np.log10(raw_max)))

    nice_vals = []
    for exp in range(exp_min - 1, exp_max + 2):
        for m in (1, 2, 3, 5):
            v = m * (10**exp)
            if raw_min <= v <= raw_max:
                nice_vals.append(v)

    nice_vals = sorted(set(nice_vals))

    # Fallback: if range is too tight, just relabel existing ticks nicely
    if len(nice_vals) < 2:
        yticks = ax.get_yticks()
        labels = [f"{int(max(0, round(10**y - 1))):,}" for y in yticks]
        ax.set_yticklabels(labels)
        return

    # Convert to axis positions in log10(value + 1)
    ticks = [np.log10(v + 1.0) for v in nice_vals]
    labels = [f"{int(v):,}" for v in nice_vals]

    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)


# ---------- I/O helpers ----------
def gather_inputs(globpat: str) -> list[Path]:
    files = sorted(Path(".").glob(globpat))
    if not files:
        raise SystemExit(f"No files matched: {globpat}")
    return files


def load_stack(h5: Path) -> pd.DataFrame:
    m = H5AD_RE.match(h5.name)
    if not m:
        raise ValueError(f"Unrecognized filename: {h5.name}")
    cell = m.group("cell")
    mark = m.group("mark")  # 'ac' or 'me3'

    ad = sc.read_h5ad(h5)
    if "tsse" not in ad.obs or "n_fragment" not in ad.obs:
        raise KeyError(f"{h5.name} is missing 'tsse' or 'n_fragment' in .obs")

    tsse = pd.to_numeric(ad.obs["tsse"], errors="coerce")
    nfrag = pd.to_numeric(ad.obs["n_fragment"], errors="coerce")
    out = pd.DataFrame({
        "CellLine": cell,
        "Mark": mark,
        "tsse": tsse.values,
        "n_fragment": nfrag.values,
    }).dropna()

    # log for plotting fragments nicely
    out["log10_n_fragment"] = np.log10(out["n_fragment"] + 1.0)
    return out


# ---------- Plotting core (matplotlib-only, explicit positions) ----------
def _plot_row(ax, df_row: pd.DataFrame, feature: str, mark: str, ylabel: str,
              annotate_counts: bool, annotate_median_raw: bool, convert_ticks_from_log=False):
    """Plot one row (one mark) with violins + transparent box overlays per cell line."""
    # Order cells for this mark
    present_cells = [c for c in CELL_ORDER_DEFAULT if c in df_row["CellLine"].unique()]
    if not present_cells:
        ax.text(0.5, 0.5, f"No data for mark={mark}", ha="center", va="center", fontsize=12)
        ax.set_axis_off()
        return

    base = np.arange(len(present_cells))
    # Data lists in that order
    data = [df_row.loc[df_row["CellLine"] == c, feature].dropna().values for c in present_cells]

    # Colors per cell line (fixed mapping for consistency)
    colors = [CELL_COLORS[c] for c in present_cells]

    # Violins
    v = ax.violinplot(
        data, positions=base, widths=0.7,
        showmeans=False, showmedians=False, showextrema=False
    )
    for body, col in zip(v["bodies"], colors):
        body.set_facecolor(col)
        body.set_edgecolor("black")
        body.set_linewidth(0.6)
        body.set_alpha(0.85)

    # Transparent box overlays
    bp = ax.boxplot(
        data, positions=base, widths=0.25,
        patch_artist=True, showfliers=False,
        medianprops=dict(color="black", linewidth=1.6),
    )
    for box in bp["boxes"]:
        box.set_facecolor("none"); box.set_edgecolor("black"); box.set_linewidth(1.4); box.set_zorder(3)
    for whisk in bp["whiskers"]:
        whisk.set_color("black"); whisk.set_linewidth(1.2); whisk.set_zorder(3)
    for cap in bp["caps"]:
        cap.set_color("black"); cap.set_linewidth(1.2); cap.set_zorder(3)
    for med in bp["medians"]:
        med.set_zorder(4)

    ax.set_xticks(base)
    ax.set_xticklabels(present_cells)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.set_title(f"H3K27{mark}")

    # Pretty ticks for fragments (log -> raw, with nice round values)
    if convert_ticks_from_log:
        fmt_log_ticks(ax)

    # Add headroom for annotations
    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin
    ax.set_ylim(ymin, ymax + 0.12 * y_range)
    y_annot = ymax + 0.06 * y_range

    # Annotations per cell
    if annotate_counts or annotate_median_raw:
        for i, cell in enumerate(present_cells):
            sub = df_row.loc[df_row["CellLine"] == cell]
            n = len(sub)
            txt = []
            if annotate_counts:
                txt.append(f"n={n:,}")
            if annotate_median_raw:
                # convert median back to raw if feature is log; otherwise raw already
                if convert_ticks_from_log:
                    med_raw = int(round(float(np.median(sub["n_fragment"].values))))
                    txt.append(f"median={med_raw:,}")
                else:
                    med_val = float(np.median(sub[feature].values))
                    txt.append(f"median={med_val:,.2f}")
            if txt:
                ax.text(base[i], y_annot, " • ".join(txt), ha="center", va="bottom",
                        fontsize=9.5, clip_on=False)

    ax.grid(axis="y", alpha=0.25)


def plot_tsse_and_fragments(df: pd.DataFrame, out_prefix: str):
    import matplotlib.pyplot as plt
    import os
    os.makedirs("./figures", exist_ok=True)

    # --- TSSE (H3K27ac ONLY) ---
    df_ac = df[df["Mark"] == "ac"]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11, 5))
    _plot_row(
        ax, df_ac, feature="tsse", mark="ac", ylabel="TSSE",
        annotate_counts=True, annotate_median_raw=False,
        convert_ticks_from_log=False
    )
    fig.suptitle("TSSE per cell line (H3K27ac only)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    pdf = f"./figures/{out_prefix}_TSSE_CellLines_H3K27acOnly.pdf"
    png = f"./figures/{out_prefix}_TSSE_CellLines_H3K27acOnly.png"
    fig.savefig(pdf); fig.savefig(png, dpi=300); plt.close(fig)
    print(pdf); print(png)

    # --- Fragment counts (ac + me3) ---
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(11, 8), sharex=False)
    for r, mark in enumerate(["ac", "me3"]):
        sub = df[df["Mark"] == mark]
        _plot_row(
            axes[r], sub, feature="log10_n_fragment", mark=mark, ylabel="Fragment Count",
            annotate_counts=True, annotate_median_raw=True,
            convert_ticks_from_log=True
        )
    fig.suptitle("Fragment counts per cell line", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    pdf = f"./figures/{out_prefix}_Fragments_CellLines_byMark.pdf"
    png = f"./figures/{out_prefix}_Fragments_CellLines_byMark.png"
    fig.savefig(pdf); plt.close(fig)
    print(pdf); print(png)


def main():
    ap = argparse.ArgumentParser(description="TSSE & Fragment-count violins per cell line (ac & me3).")
    ap.add_argument("--inputs", default="*_*_CountFilteredShared.h5ad",
                    help="Glob for input .h5ad (default: %(default)s)")
    ap.add_argument("--out", default="DNA_CellLines", help="Output prefix")
    args = ap.parse_args()

    files = gather_inputs(args.inputs)

    dfs = []
    for h5 in files:
        try:
            dfs.append(load_stack(h5))
        except Exception as e:
            print(f"[skip] {h5.name}: {e}")
    if not dfs:
        raise SystemExit("No usable inputs.")

    df = pd.concat(dfs, ignore_index=True)

    # Keep only known cells to ensure stable colors
    df = df[df["CellLine"].isin(CELL_ORDER_DEFAULT)]

    plot_tsse_and_fragments(df, args.out)


if __name__ == "__main__":
    main()

python lp.py 