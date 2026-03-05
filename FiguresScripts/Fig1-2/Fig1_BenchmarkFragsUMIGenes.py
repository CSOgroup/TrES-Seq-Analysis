#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Benchmark fragments-per-cell (DNA) and UMI/genes-per-cell (RNA) across technologies.

Key fixes (from previous step):
- EXACT intersection on raw .obs_names (no custom transforms).
- Sanitize barcodes (strip) + drop duplicate obs_names with reports.
- Subset by integer positions; assert intersect counts match across modalities.

Outputs:
  - BOTH plot sets are always produced:
      *_intersect*.pdf  -> strict pass-in-all-modalities cells
      *_full*.pdf       -> per-modality filtered cells (no cross-modality restriction)
  - NEW: Combined TSSE violin for H3K27ac (both full + intersect):
      Benchmark_Combined_H3K27ac_TSSE_Violin_full.pdf
      Benchmark_Combined_H3K27ac_TSSE_Violin_intersect.pdf
  - Intersection debug files and duplicate reports as before.

Intersection rule:
  - InHouse_NanoCNT, NanoCNT, NTT:  H3K27ac ∩ H3K27me3
  - MTR (MTR folder), TrES:         H3K27ac ∩ H3K27me3 ∩ RNA
  - 10xMultiome: RNA-only; never intersected with DNA (included in _full RNA plots)
"""

import os
import glob
import argparse
import warnings
from typing import Dict, Tuple, Optional, List, Set

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

import snapatac2 as snap
from matplotlib_venn import venn2, venn3

# ------------------------- CLI -------------------------

def get_args():
    ap = argparse.ArgumentParser(
        description="Benchmark DNA fragments/cell and RNA UMIs/genes per cell across technologies (full + intersect)."
    )
    ap.add_argument("--base-dir", default=".", help="Directory containing InHouse_NanoCNT, MTR, NanoCNT, NTT, TrES, 10xMultiome")
    ap.add_argument("--threads", type=int, default=24, help="Threads for SnapATAC2")
    ap.add_argument("--figdir", default="figures", help="Output folder for figures")
    ap.add_argument("--prefix", default="Benchmark", help="Filename prefix for figures")

    # DNA thresholds
    ap.add_argument("--ac-min", type=int, default=250)
    ap.add_argument("--ac-max", type=int, default=15000)
    ap.add_argument("--ac-min-tsse", type=float, default=0.0)

    ap.add_argument("--me3-min", type=int, default=500)
    ap.add_argument("--me3-max", type=int, default=15000)
    ap.add_argument("--me3-min-tsse", type=float, default=0.0)

    # RNA thresholds
    ap.add_argument("--rna-min-umi", type=int, default=1000)
    ap.add_argument("--rna-max-umi", type=int, default=25000)
    ap.add_argument("--rna-max-mt", type=float, default=10.0)      # %
    ap.add_argument("--rna-max-hb", type=float, default=5.0)       # %
    ap.add_argument("--rna-max-top50", type=float, default=30.0)   # %

    # Partial reruns
    ap.add_argument("--skip-dna", action="store_true")
    ap.add_argument("--skip-rna", action="store_true")

    # Caching for fragment imports
    ap.add_argument("--force-reimport", action="store_true")
    ap.add_argument("--no-save-cache", action="store_true")

    return ap.parse_args()

# ------------------------- Utils -------------------------

def ensure_dir(d: str):
    if not os.path.exists(d):
        os.makedirs(d, exist_ok=True)

def pick_first(patterns: List[str]) -> Optional[str]:
    for p in patterns:
        hits = sorted(glob.glob(p))
        if hits:
            return hits[0]
    return None

def tech_palette() -> Dict[str, str]:
    return {
        "TrES": "#66c2a5",
        "MTR":  "#3288bd",     # renamed display for MTR
        "MTR":  "#3288bd",
        "InHouse_NanoCNT": "#e78ac3",
        "NanoCNT": "#a6d854",
        "NTT": "#ffd92f",
        "10xMultiome": "#e6ab02",
    }

def display_name(tech: str) -> str:
    return "MTR" if tech == "MTR" else tech

def genome_for_tech(tech: str):
    return snap.genome.mm10 if tech == "NanoCNT" else snap.genome.hg38

def add_log10_frag(adata: ad.AnnData):
    if "n_fragment" not in adata.obs:
        raise RuntimeError("AnnData missing obs['n_fragment']")
    adata.obs["log10_n_fragment"] = np.log10(adata.obs["n_fragment"] + 1)

def fmt_log_ticks(ax):
    """
    Use nice round ticks for log10(value + 1) violins without forcing the axis to start at 0.

    - Axis is in log10(value + 1) space.
    - We keep the current y-limits (matplotlib defaults).
    - Ticks use round numbers: 1, 2, 3, 5 × 10^k within the data range (so e.g. 2000, 3000, 5000).
    """
    ymin, ymax = ax.get_ylim()

    # Convert current limits back to raw scale (approximate)
    raw_min = max(0.0, 10**ymin - 1.0)
    raw_max = max(0.0, 10**ymax - 1.0)

    # If upper bound is non-positive, nothing sensible to do -> keep defaults
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


def violin_with_box(ax, data, color=None, order=None):
    sns.violinplot(data=data, x="Technology", y="Value",
                   palette=None if color is None else color,
                   inner=None, order=order, ax=ax)
    sns.boxplot(data=data, x="Technology", y="Value", order=order, ax=ax,
                width=0.2, showcaps=True, dodge=False,
                boxprops={'facecolor': 'none', 'linewidth': 1.5},
                whiskerprops={'linewidth': 1.5},
                medianprops={'color': 'black'},
                flierprops={'marker': '', 'markersize': 0})

def annotate_counts_and_median(ax, df, y_is_log=True):
    ymin, ymax = ax.get_ylim()
    headroom = 0.15 * (ymax - ymin)
    ax.set_ylim(ymin, ymax + headroom)
    ypos = ymax + 0.05 * (ymax - ymin)
    for xt, tick in enumerate(ax.get_xticklabels()):
        tech_name = tick.get_text()
        sub = df[df["Technology"] == tech_name]["Value"].values
        if sub.size == 0:
            continue
        med = float(np.median(sub))
        med_raw = int(max(0, round(10 ** med - 1))) if y_is_log else int(round(med))
        ax.text(xt, ypos, f"{sub.size} cells\n{med_raw:,} median", ha="center", va="bottom", fontsize=10)

def cached_h5ad_path(fragment_path: str) -> str:
    p = fragment_path
    if p.endswith(".gz"):
        p = p[:-3]
    root, _ext = os.path.splitext(p)
    return root + ".imported.h5ad"

def ids(A: ad.AnnData) -> Set[str]:
    return set(map(str, A.obs_names.to_list()))

# ---------- Sanitize barcodes & enforce uniqueness ----------

def _write_duplicates_report(adata: ad.AnnData, label: str, outdir: str, prefix: str):
    dup = adata.obs_names[adata.obs_names.duplicated()].to_list()
    if not dup:
        return
    vc = pd.Series(dup).value_counts().sort_values(ascending=False)
    path = os.path.join(outdir, f"{prefix}_{label}_DUPLICATE_obsnames.txt")
    with open(path, "w") as fh:
        fh.write(f"# Duplicated obs_names for {label}\n")
        fh.write("# count\tobs_name\n")
        for name, cnt in vc.items():
            fh.write(f"{cnt}\t{name}\n")
    print(f"[WARN] Duplicates detected and dropped in {label}: {len(dup)} entries. Report: {path}")

def sanitize_obsnames(adata: ad.AnnData, label: str, outdir: str, prefix: str) -> ad.AnnData:
    names = pd.Index(adata.obs_names.astype(str))
    cleaned = names.str.strip()
    if not cleaned.equals(names):
        adata = adata.copy()
        adata.obs_names = cleaned
    if (adata.obs_names == "").any():
        mask = adata.obs_names != ""
        dropped = int((~mask).sum())
        adata = adata[mask, :].copy()
        print(f"[WARN] Dropped {dropped} empty obs_names in {label}.")
    if adata.obs_names.has_duplicates:
        _write_duplicates_report(adata, label, outdir, prefix)
        keep_mask = ~adata.obs_names.duplicated(keep="first")
        adata = adata[keep_mask, :].copy()
    if adata.obs_names.has_duplicates:
        raise RuntimeError(f"{label}: obs_names still not unique after sanitization.")
    return adata

# ------------------------- DNA I/O -------------------------

def detect_dna_files(base_dir: str) -> Dict[str, Dict[str, Tuple[str, str]]]:
    mapping = {}

    def add_entry(tech, mark, path, ftype):
        mapping.setdefault(tech, {})
        mapping[tech][mark] = (path, ftype)

    # TrES DNA
    tres_dir = os.path.join(base_dir, "TrES")
    if os.path.isdir(tres_dir):
        ac_h5 = pick_first([os.path.join(tres_dir, "H2RD_Human_H3K27ac.h5ad"),
                            os.path.join(tres_dir, "*H3K27ac*.h5ad")])
        me3_h5 = pick_first([os.path.join(tres_dir, "H2RD_Human_H3K27me3.h5ad"),
                             os.path.join(tres_dir, "*H3K27me3*.h5ad")])
        if ac_h5:  add_entry("TrES", "ac", ac_h5, "h5ad")
        if me3_h5: add_entry("TrES", "me3", me3_h5, "h5ad")

    # InHouse_NanoCNT
    ih_dir = os.path.join(base_dir, "InHouse_NanoCNT")
    if os.path.isdir(ih_dir):
        ac = pick_first([os.path.join(ih_dir, "ScKDMA_H3K27ac.tsv.gz"),
                         os.path.join(ih_dir, "*H3K27ac*.tsv.gz"),
                         os.path.join(ih_dir, "*H3K27ac*.bed.gz")])
        me3 = pick_first([os.path.join(ih_dir, "ScKDMA_H3K27me3.tsv.gz"),
                          os.path.join(ih_dir, "*H3K27me3*.tsv.gz"),
                          os.path.join(ih_dir, "*H3K27me3*.bed.gz")])
        if ac and os.path.dirname(ac) == ih_dir:   add_entry("InHouse_NanoCNT", "ac", ac, "fragments")
        if me3 and os.path.dirname(me3) == ih_dir: add_entry("InHouse_NanoCNT", "me3", me3, "fragments")
        if "ac" not in mapping.get("InHouse_NanoCNT", {}):
            ac_h5 = pick_first([os.path.join(ih_dir, "filtered", "*H3K27ac*filtered*.h5ad")])
            if ac_h5: add_entry("InHouse_NanoCNT", "ac", ac_h5, "h5ad")
        if "me3" not in mapping.get("InHouse_NanoCNT", {}):
            me3_h5 = pick_first([os.path.join(ih_dir, "filtered", "*H3K27me3*filtered*.h5ad")])
            if me3_h5: add_entry("InHouse_NanoCNT", "me3", me3_h5, "h5ad")

    # MTR
    mtc_dir = os.path.join(base_dir, "MTR")
    if os.path.isdir(mtc_dir):
        ac = pick_first([os.path.join(mtc_dir, "*H3K27ac*.fragments.bed.gz"),
                         os.path.join(mtc_dir, "*H3K27ac*.fragments.tsv.gz"),
                         os.path.join(mtc_dir, "*H3K27ac*.bed.gz"),
                         os.path.join(mtc_dir, "*H3K27ac*.tsv.gz")])
        me3 = pick_first([os.path.join(mtc_dir, "*H3K27me3*.fragments.bed.gz"),
                          os.path.join(mtc_dir, "*H3K27me3*.fragments.tsv.gz"),
                          os.path.join(mtc_dir, "*H3K27me3*.bed.gz"),
                          os.path.join(mtc_dir, "*H3K27me3*.tsv.gz")])
        if ac:  add_entry("MTR", "ac", ac, "fragments")
        if me3: add_entry("MTR", "me3", me3, "fragments")

    # NanoCNT
    nc_dir = os.path.join(base_dir, "NanoCNT")
    if os.path.isdir(nc_dir):
        ac = pick_first([os.path.join(nc_dir, "*H3K27ac*frag*.tsv.gz"),
                         os.path.join(nc_dir, "*H3K27ac*.tsv.gz"),
                         os.path.join(nc_dir, "*H3K27ac*.bed.gz")])
        me3 = pick_first([os.path.join(nc_dir, "*H3K27me3*frag*.tsv.gz"),
                          os.path.join(nc_dir, "*H3K27me3*.tsv.gz"),
                          os.path.join(nc_dir, "*H3K27me3*.bed.gz")])
        if ac:  add_entry("NanoCNT", "ac", ac, "fragments")
        if me3: add_entry("NanoCNT", "me3", me3, "fragments")

    # NTT
    ntt_dir = os.path.join(base_dir, "NTT")
    if os.path.isdir(ntt_dir):
        ac = pick_first([os.path.join(ntt_dir, "*H3K27ac*.bed.gz"),
                         os.path.join(ntt_dir, "*H3K27ac*.tsv.gz")])
        me3 = pick_first([os.path.join(ntt_dir, "*H3K27me3*.bed.gz"),
                          os.path.join(ntt_dir, "*H3K27me3*.tsv.gz")])
        if ac:  add_entry("NTT", "ac", ac, "fragments")
        if me3: add_entry("NTT", "me3", me3, "fragments")

    return mapping

def load_one_dna(path: str, ftype: str, tech: str, tmpdir: str, threads: int,
                 use_cache: bool = True, force_reimport: bool = False, save_cache: bool = True) -> ad.AnnData:
    g = genome_for_tech(tech)
    if ftype == "h5ad":
        adata = sc.read(path)
        if "tsse" not in adata.obs.columns:
            try:
                snap.metrics.tsse(adata, gene_anno=g, inplace=True)
            except Exception as e:
                warnings.warn(f"Could not compute TSSE for {os.path.basename(path)}: {e}")
        if "n_fragment" not in adata.obs.columns:
            raise RuntimeError(f"{path} lacks obs['n_fragment'].")
        add_log10_frag(adata)
        return adata

    elif ftype == "fragments":
        cache_path = cached_h5ad_path(path)
        if use_cache and not force_reimport and os.path.exists(cache_path):
            print(f"  -> Using cached AnnData: {cache_path}")
            adata = sc.read(cache_path)
            if "n_fragment" not in adata.obs.columns:
                warnings.warn(f"Cached {cache_path} missing 'n_fragment'; re-importing.")
            else:
                if "tsse" not in adata.obs.columns:
                    try:
                        snap.metrics.tsse(adata, gene_anno=g, inplace=True)
                    except Exception as e:
                        warnings.warn(f"TSSE on cached failed: {e}")
                add_log10_frag(adata)
                return adata

        print(f"  -> Importing fragments (unsorted) from {path}")
        adata = snap.pp.import_fragments(
            fragment_file=path,
            chrom_sizes=g,
            min_num_fragments=0,
            sorted_by_barcode=False,
            chrM=['chrM', 'M'],
            tempdir=tmpdir,
            backend="hdf5",
            n_jobs=threads
        )
        try:
            snap.metrics.tsse(adata, gene_anno=g, inplace=True)
        except Exception as e:
            warnings.warn(f"TSSE computation failed for {os.path.basename(path)}: {e}")
        add_log10_frag(adata)

        if save_cache:
            try:
                adata.write(cache_path)
                print(f"  -> Cached imported AnnData at {cache_path}")
            except Exception as e:
                warnings.warn(f"Could not write cache {cache_path}: {e}")

        return adata

    else:
        raise ValueError("ftype must be 'h5ad' or 'fragments'.")

def filter_dna(adata: ad.AnnData, min_frags: int, max_frags: int, min_tsse: float) -> ad.AnnData:
    mask = (adata.obs["n_fragment"] >= min_frags) & (adata.obs["n_fragment"] <= max_frags)
    if "tsse" in adata.obs.columns and min_tsse > 0:
        mask = mask & (adata.obs["tsse"] >= min_tsse)
    adata_f = adata[mask].copy()
    add_log10_frag(adata_f)
    return adata_f

# ------------------------- RNA I/O -------------------------

def detect_rna_files(base_dir: str) -> Dict[str, Tuple[str, str]]:
    out = {}

    # TrES RNA
    tres_dir = os.path.join(base_dir, "TrES")
    if os.path.isdir(tres_dir):
        tres_raw = pick_first([os.path.join(tres_dir, "H2RRNA_Human.h5ad"),
                               os.path.join(tres_dir, "*RNA*.h5ad"),
                               os.path.join(tres_dir, "*scRNA*.h5ad")])
        if tres_raw and os.path.dirname(tres_raw) == tres_dir:
            out["TrES"] = (tres_raw, "h5ad")

    # MTR RNA (display as MTR)
    mtc_dir = os.path.join(base_dir, "MTR")
    if os.path.isdir(mtc_dir):
        mtc_h5 = pick_first([os.path.join(mtc_dir, "scRNA_D0_withCellNames_likeScanpy.h5ad"),
                             os.path.join(mtc_dir, "*scRNA*.h5ad"),
                             os.path.join(mtc_dir, "*.h5ad")])
        if mtc_h5 and os.path.dirname(mtc_h5) == mtc_dir:
            out["MTR"] = (mtc_h5, "h5ad")

    # 10xMultiome RNA (optional)
    mx_dir = os.path.join(base_dir, "10xMultiome")
    if os.path.isdir(mx_dir):
        tenx_h5 = pick_first([os.path.join(mx_dir, "*filtered_feature_bc_matrix.h5"),
                              os.path.join(mx_dir, "*.h5")])
        if tenx_h5:
            out["10xMultiome"] = (tenx_h5, "10x_h5")

    return out

def _subset_to_rna_if_present(adata: ad.AnnData) -> ad.AnnData:
    if "feature_types" in adata.var:
        ft = adata.var["feature_types"].astype(str)
        mask = ft.str.contains(r"Gene Expression|RNA", case=False, regex=True)
        if mask.any():
            adata = adata[:, mask.values].copy()
    return adata

def load_one_rna(path_or_dir: str, kind: str) -> ad.AnnData:
    if kind == "h5ad":
        return sc.read(path_or_dir)
    if kind == "10x_h5":
        adata = sc.read_10x_h5(path_or_dir)
        return _subset_to_rna_if_present(adata)
    raise ValueError("Unknown RNA kind (expected .h5ad or 10x_h5).")

def compute_rna_qc(adata: ad.AnnData):
    gs = adata.var_names.astype(str)
    gs_u = gs.str.upper()
    adata.var["mt"] = gs_u.str.startswith("MT-")
    adata.var["hb"] = gs_u.str.match(r"^HB[AB]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","hb"], percent_top=[50], inplace=True)
    adata.obs.rename(columns={
        "total_counts": "UMIs",
        "n_genes_by_counts": "n_genes",
        "pct_counts_mt": "pct_mt",
        "pct_counts_hb": "pct_hb",
        "pct_counts_in_top_50_genes": "pct_top50",
    }, inplace=True)

def filter_rna(adata: ad.AnnData, minUMI, maxUMI, max_mt, max_hb, max_top50) -> ad.AnnData:
    qc = adata.obs
    keep = (
        (qc["UMIs"] >= minUMI) & (qc["UMIs"] <= maxUMI) &
        (qc["pct_mt"] <= max_mt) & (qc["pct_hb"] <= max_hb) &
        (qc["pct_top50"] <= max_top50)
    )
    return adata[keep].copy()

# ------------------------- Plotting helpers -------------------------

def plot_violin_single_frag(A: ad.AnnData, tech_disp: str, mark: str, outdir: str, prefix: str, tag: str):
    title = f"{tech_disp} H3K27{mark} ({'intersect' if tag=='intersect' else 'full'})"
    outpath = os.path.join(outdir, f"{prefix}_{tech_disp}_H3K27{mark}_FragViolin_{tag}.pdf")
    fig, ax = plt.subplots(figsize=(4.8, 6))
    data = pd.DataFrame({"Technology": [tech_disp]*A.n_obs, "Value": A.obs["log10_n_fragment"].values})
    violin_with_box(ax, data, color=None)
    fmt_log_ticks(ax)
    ax.set_xlabel("")
    ax.set_ylabel("Fragments per cell")
    med_raw = int(np.median(A.obs["n_fragment"]))
    ax.set_title(f"{title}\nCells: {A.n_obs}\nMedian fragments: {med_raw:,}")
    plt.tight_layout()
    fig.savefig(outpath); plt.close(fig)

def plot_scatter_tsse(A: ad.AnnData, tech_disp: str, mark: str, outdir: str, prefix: str, tag: str):
    if "tsse" not in A.obs:
        return
    title = f"{tech_disp} H3K27{mark} ({'intersect' if tag=='intersect' else 'full'})"
    outpath = os.path.join(outdir, f"{prefix}_{tech_disp}_H3K27{mark}_TSSE_vs_LogFrags_{tag}.pdf")
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.scatter(A.obs["log10_n_fragment"], A.obs["tsse"], s=6, alpha=0.5)
    ax.set_xlabel("log10(Fragments+1)")
    ax.set_ylabel("TSSE")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outpath); plt.close(fig)

def combined_violin_frag(d: Dict[str, ad.AnnData], mark: str, outdir: str, prefix: str, tag: str):
    rows = []
    for tech_disp, A in d.items():
        for v in A.obs["log10_n_fragment"].values:
            rows.append({"Technology": tech_disp, "Value": v})
    if not rows:
        return
    df = pd.DataFrame(rows)
    pal = tech_palette()
    order = [t for t in ["TrES", "MTR", "InHouse_NanoCNT", "NanoCNT", "NTT"] if t in df["Technology"].unique()]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    sns.violinplot(data=df, x="Technology", y="Value", palette=pal, inner=None, order=order, ax=ax)
    sns.boxplot(data=df, x="Technology", y="Value", order=order, ax=ax, width=0.2,
                showcaps=True, dodge=False,
                boxprops={'facecolor':'none','linewidth':1.5},
                whiskerprops={'linewidth':1.5}, medianprops={'color':'black'},
                flierprops={'marker':'','markersize':0})
    fmt_log_ticks(ax)
    annotate_counts_and_median(ax, df, y_is_log=True)
    ax.set_xlabel("")
    ax.set_ylabel("Fragments per cell")
    title_tag = "intersected cells" if tag == "intersect" else "full (per-modality filtered)"
    ax.set_title(f"H3K27{mark} – Fragments per Cell ({title_tag})")
    plt.tight_layout()
    out = os.path.join(outdir, f"{prefix}_Combined_H3K27{mark}_Fragments_Violin_{tag}.pdf")
    fig.savefig(out); plt.close(fig)

# --- NEW: Combined TSSE violin for H3K27ac (full + intersect) ---
def combined_violin_tsse_ac(dna_sets: Dict[str, Dict[str, ad.AnnData]], outdir: str, prefix: str, tag: str):
    """
    Build a combined TSSE violin across technologies for H3K27ac only.
    Skips datasets without a 'tsse' column.
    """
    rows = []
    for tech_disp, marks in dna_sets.items():
        if "ac" not in marks:
            continue
        A = marks["ac"]
        if A is None or A.n_obs == 0 or "tsse" not in A.obs:
            continue
        tsse_vals = pd.to_numeric(A.obs["tsse"], errors="coerce").dropna().values
        rows.extend({"Technology": tech_disp, "Value": float(v)} for v in tsse_vals)

    if not rows:
        return

    df = pd.DataFrame(rows)
    pal = tech_palette()
    order = [t for t in ["TrES", "MTR", "InHouse_NanoCNT", "NanoCNT", "NTT"] if t in df["Technology"].unique()]

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    sns.violinplot(data=df, x="Technology", y="Value", palette=pal, inner=None, order=order, ax=ax)
    sns.boxplot(data=df, x="Technology", y="Value", order=order, ax=ax, width=0.2,
                showcaps=True, dodge=False,
                boxprops={'facecolor':'none','linewidth':1.5},
                whiskerprops={'linewidth':1.5}, medianprops={'color':'black'},
                flierprops={'marker':'','markersize':0})
    # custom float annotations for TSSE
    ymin, ymax = ax.get_ylim()
    headroom = 0.15 * (ymax - ymin)
    ax.set_ylim(ymin, ymax + headroom)
    ypos = ymax + 0.05 * (ymax - ymin)
    for xt, tick in enumerate(ax.get_xticklabels()):
        tech_name = tick.get_text()
        sub = df[df["Technology"] == tech_name]["Value"].values
        if sub.size == 0:
            continue
        med = float(np.median(sub))
        ax.text(xt, ypos, f"{sub.size} cells\n{med:.2f} median", ha="center", va="bottom", fontsize=10)

    ax.set_xlabel("")
    ax.set_ylabel("TSSE")
    title_tag = "intersected cells" if tag == "intersect" else "full (per-modality filtered)"
    ax.set_title(f"H3K27ac – TSSE ({title_tag})")
    plt.tight_layout()
    out = os.path.join(outdir, f"{prefix}_Combined_H3K27ac_TSSE_Violin_{tag}.pdf")
    fig.savefig(out); plt.close(fig)

def rna_violin(rna_map_for_plot: Dict[str, ad.AnnData], metric: str, outdir: str, prefix: str, tag: str):
    rows = []
    for tech_disp, A in rna_map_for_plot.items():
        if A is None or metric not in A.obs.columns:
            continue
        vals = np.log10(A.obs[metric].to_numpy() + 1)
        rows.extend({"Technology": tech_disp, "Value": v} for v in vals)
    if not rows:
        return
    df = pd.DataFrame(rows)
    pal = tech_palette()
    order = [t for t in ["TrES", "MTR", "10xMultiome"] if t in df["Technology"].unique()]
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    sns.violinplot(data=df, x="Technology", y="Value", palette=pal, inner=None, order=order, ax=ax)
    sns.boxplot(data=df, x="Technology", y="Value", order=order, ax=ax, width=0.2,
                showcaps=True, dodge=False,
                boxprops={'facecolor':'none','linewidth':1.5},
                whiskerprops={'linewidth':1.5}, medianprops={'color':'black'},
                flierprops={'marker':'','markersize':0})

    fmt_log_ticks(ax)

    # Force a 10,000 tick for n_genes
    if metric == "n_genes":
        forced_raw_tick = 10000
        forced_tick = np.log10(forced_raw_tick + 1.0)
        yticks = list(ax.get_yticks())
        if not any(np.isclose(t, forced_tick) for t in yticks):
            yticks.append(forced_tick)
            yticks = sorted(yticks)
        ax.set_yticks(yticks)
        ax.set_yticklabels(
            [f"{int(max(0, round(10**y - 1))):,}" for y in yticks]
        )

    annotate_counts_and_median(ax, df, y_is_log=True)
    ax.set_xlabel("")
    ax.set_ylabel(f"{metric} per cell")
    title_tag = "intersected cells" if tag == "intersect" else "full (per-modality filtered)"
    ax.set_title(f"RNA – {metric} per Cell ({title_tag})")
    plt.tight_layout()
    out = os.path.join(outdir, f"{prefix}_RNA_{metric}_Violin_{tag}.pdf")
    fig.savefig(out); plt.close(fig)


def venn2_ac_me3_filtered(ac: ad.AnnData, me3: ad.AnnData, tech_disp: str, outdir: str, prefix: str):
    set1 = set(ac.obs_names)
    set2 = set(me3.obs_names)
    fig = plt.figure(figsize=(6.2, 6.0))
    venn2([set1, set2], set_labels=(f"{tech_disp} H3K27ac", f"{tech_disp} H3K27me3"))
    plt.title(f"Overlap of Passing Cells – {tech_disp}")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, f"{prefix}_{tech_disp}_Venn_ac_vs_me3.pdf")); plt.close(fig)

def venn3_modalities_filtered(ac: ad.AnnData, me3: ad.AnnData, rna: ad.AnnData,
                              tech_disp: str, outdir: str, prefix: str):
    A = set(ac.obs_names); B = set(me3.obs_names); C = set(rna.obs_names)
    fig = plt.figure(figsize=(7.0, 7.0))
    venn3([A, B, C], set_labels=("H3K27ac", "H3K27me3", "RNA"))
    plt.title(f"Venn of Passing Cells – {tech_disp} (ac / me3 / RNA)")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, f"{prefix}_{tech_disp}_Venn3_ac_me3_RNA.pdf")); plt.close(fig)

# ------------------------- Main -------------------------

def main():
    args = get_args()
    ensure_dir(args.figdir)
    tmpdir = os.path.join(args.base_dir, "_tmp_snap")
    ensure_dir(tmpdir)

    skip_dna = getattr(args, "skip_dna", False)
    skip_rna = getattr(args, "skip_rna", False)

    # Detect inputs
    dna_map = detect_dna_files(args.base_dir)
    rna_files = detect_rna_files(args.base_dir)

    # ---- Load RNA (MTR/TrES; 10xMultiome optional) ----
    rna_flt_raw: Dict[str, ad.AnnData] = {}
    if not skip_rna:
        for tech in ["MTR", "TrES"]:
            if tech in rna_files:
                path_or_dir, kind = rna_files[tech]
                print(f"[RNA] Loading {tech} from {path_or_dir} ({kind})")
                R = load_one_rna(path_or_dir, kind)
                R = _subset_to_rna_if_present(R)
                compute_rna_qc(R)
                Rf = filter_rna(R, args.rna_min_umi, args.rna_max_umi,
                                args.rna_max_mt, args.rna_max_hb, args.rna_max_top50)
                Rf = sanitize_obsnames(Rf, label=f"{display_name(tech)}_RNA", outdir=args.figdir, prefix=args.prefix)
                rna_flt_raw[tech] = Rf
                print(f"  -> {display_name(tech)} RNA passing (unique IDs): {Rf.n_obs}")

    rna_10x: Optional[ad.AnnData] = None
    if (not skip_rna) and ("10xMultiome" in rna_files):
        tenx_path, tenx_kind = rna_files["10xMultiome"]
        print(f"[RNA] Loading 10xMultiome from {tenx_path} ({tenx_kind})")
        Rmx = load_one_rna(tenx_path, tenx_kind)
        Rmx = _subset_to_rna_if_present(Rmx)
        compute_rna_qc(Rmx)
        rna_10x = filter_rna(
            Rmx,
            args.rna_min_umi,
            args.rna_max_umi,
            args.rna_max_mt,
            args.rna_max_hb,
            args.rna_max_top50,
        )
        rna_10x = sanitize_obsnames(rna_10x, label="10xMultiome_RNA", outdir=args.figdir, prefix=args.prefix)
        print(f"  -> 10xMultiome RNA passing (unique IDs): {rna_10x.n_obs}")

    # ---- Load & filter DNA ----
    dna_flt_raw: Dict[str, Dict[str, ad.AnnData]] = {}
    if not skip_dna:
        for tech, marks in dna_map.items():
            dna_flt_raw.setdefault(tech, {})
            for mark, (path, ftype) in marks.items():
                print(f"[DNA] Loading {tech} {mark} from {path} ({ftype})")
                A = load_one_dna(
                    path, ftype, tech, tmpdir, args.threads,
                    use_cache=not args.force_reimport,
                    force_reimport=args.force_reimport,
                    save_cache=not args.no_save_cache
                )
                if mark == "ac":
                    Af = filter_dna(A, args.ac_min, args.ac_max, args.ac_min_tsse)
                else:
                    Af = filter_dna(A, args.me3_min, args.me3_max, args.me3_min_tsse)
                Af = sanitize_obsnames(Af, label=f"{display_name(tech)}_H3K27{mark}", outdir=args.figdir, prefix=args.prefix)
                dna_flt_raw[tech][mark] = Af
                print(f"  -> {display_name(tech)} {mark}: {Af.n_obs} cells after filters (unique IDs)")

    # ---------------- FULL maps (no intersection) ----------------
    dna_full_disp: Dict[str, Dict[str, ad.AnnData]] = {}
    rna_full_disp: Dict[str, ad.AnnData] = {}

    for tech_raw, marks in dna_flt_raw.items():
        tech_disp = display_name(tech_raw)
        dna_full_disp[tech_disp] = {mk: ad for mk, ad in marks.items()}

    for tech_raw, R in rna_flt_raw.items():
        rna_full_disp[display_name(tech_raw)] = R
    if rna_10x is not None:
        rna_full_disp["10xMultiome"] = rna_10x

    # ---------------- INTERSECT (strict) maps ----------------
    dna_intersect_disp: Dict[str, Dict[str, ad.AnnData]] = {}
    rna_intersect_disp: Dict[str, ad.AnnData] = {}

    def ids(A: ad.AnnData) -> Set[str]:
        return set(map(str, A.obs_names.to_list()))

    def compute_common_ids_for_tech(tech_raw: str) -> Set[str]:
        if tech_raw in ["InHouse_NanoCNT", "NanoCNT", "NTT"]:
            if ("ac" not in dna_flt_raw.get(tech_raw, {})) or ("me3" not in dna_flt_raw.get(tech_raw, {})):
                return set()
            return ids(dna_flt_raw[tech_raw]["ac"]) & ids(dna_flt_raw[tech_raw]["me3"])
        if tech_raw in ["MTR", "TrES"]:
            if (("ac" not in dna_flt_raw.get(tech_raw, {})) or
                ("me3" not in dna_flt_raw.get(tech_raw, {})) or
                (tech_raw not in rna_flt_raw)):
                return set()
            return ids(dna_flt_raw[tech_raw]["ac"]) & ids(dna_flt_raw[tech_raw]["me3"]) & ids(rna_flt_raw[tech_raw])
        return set()

    def debug_intersection(tech_raw: str, common: Set[str], outdir: str, prefix: str):
        tech_disp = display_name(tech_raw)
        lines = [f"TECH: {tech_disp}", f"Common n = {len(common)}"]
        if "ac" in dna_flt_raw.get(tech_raw, {}):
            aset = ids(dna_flt_raw[tech_raw]["ac"])
            lines.append(f"[ac]  have={len(aset)}  missing_from_mod={len(common - aset)}  extra_in_mod_vs_common={len(aset - common)}")
        if "me3" in dna_flt_raw.get(tech_raw, {}):
            mset = ids(dna_flt_raw[tech_raw]["me3"])
            lines.append(f"[me3] have={len(mset)}  missing_from_mod={len(common - mset)}  extra_in_mod_vs_common={len(mset - common)}")
        if tech_raw in rna_flt_raw:
            rset = ids(rna_flt_raw[tech_raw])
            lines.append(f"[rna] have={len(rset)}  missing_from_mod={len(common - rset)}  extra_in_mod_vs_common={len(rset - common)}")
        text = "\n".join(lines) + "\n"
        print(text)
        with open(os.path.join(outdir, f"{prefix}_{tech_disp}_intersection_debug.txt"), "w") as fh:
            fh.write(text)

    def subset_exact(A: ad.AnnData, common: Set[str]) -> ad.AnnData:
        name_to_pos = {name: i for i, name in enumerate(A.obs_names)}
        common_sorted = sorted(common)
        missing = [name for name in common_sorted if name not in name_to_pos]
        if missing:
            raise RuntimeError(f"{len(missing)} intersection IDs not found in AnnData (first 10): {missing[:10]}")
        pos = np.fromiter((name_to_pos[n] for n in common_sorted), dtype=int)
        return A[pos, :].copy()

    def write_post_subset_counts(tech_raw: str, tech_disp: str, outdir: str, prefix: str,
                                 ac_sub: Optional[ad.AnnData], me3_sub: Optional[ad.AnnData], rna_sub: Optional[ad.AnnData]):
        lines = [f"TECH: {tech_disp}", "Post-subset counts (should all match):"]
        if ac_sub is not None:  lines.append(f"  ac : {ac_sub.n_obs}")
        if me3_sub is not None: lines.append(f"  me3: {me3_sub.n_obs}")
        if rna_sub is not None: lines.append(f"  rna: {rna_sub.n_obs}")
        text = "\n".join(lines) + "\n"
        print(text)
        with open(os.path.join(outdir, f"{prefix}_{tech_disp}_intersection_postsubset.txt"), "w") as fh:
            fh.write(text)

    for tech_raw in sorted(dna_flt_raw.keys()):
        tech_disp = display_name(tech_raw)
        common = compute_common_ids_for_tech(tech_raw)
        print(f"[{tech_disp}] STRICT intersection size = {len(common)}")
        if len(common) == 0:
            continue

        debug_intersection(tech_raw, common, args.figdir, args.prefix)

        dna_intersect_disp.setdefault(tech_disp, {})
        ac_sub = me3_sub = rna_sub = None

        if "ac" in dna_flt_raw.get(tech_raw, {}):
            ac_sub = subset_exact(dna_flt_raw[tech_raw]["ac"], common)
            dna_intersect_disp[tech_disp]["ac"] = ac_sub
        if "me3" in dna_flt_raw.get(tech_raw, {}):
            me3_sub = subset_exact(dna_flt_raw[tech_raw]["me3"], common)
            dna_intersect_disp[tech_disp]["me3"] = me3_sub
        if tech_raw in ["MTR", "TrES"] and tech_raw in rna_flt_raw:
            rna_sub = subset_exact(rna_flt_raw[tech_raw], common)
            rna_intersect_disp[tech_disp] = rna_sub

        sizes = [x.n_obs for x in [ac_sub, me3_sub, rna_sub] if x is not None]
        if len(set(sizes)) != 1:
            raise RuntimeError(f"[{tech_disp}] Intersect mismatch after subsetting: counts={sizes} (expect identical).")

        write_post_subset_counts(tech_raw, tech_disp, args.figdir, args.prefix, ac_sub, me3_sub, rna_sub)

    # ---------------- Plot BOTH sets ----------------
    def run_plots(dna_sets: Dict[str, Dict[str, ad.AnnData]],
                  rna_sets: Dict[str, ad.AnnData],
                  tag: str):
        # Per-tech/per-mark QC plots (DNA)
        for tech_disp, marks in dna_sets.items():
            for mark, A in marks.items():
                if A is None or A.n_obs == 0:
                    continue
                plot_violin_single_frag(A, tech_disp, mark, args.figdir, args.prefix, tag)
                plot_scatter_tsse(A, tech_disp, mark, args.figdir, args.prefix, tag)
        # Combined DNA violins per mark
        for mark in ["ac", "me3"]:
            subset = {tech_disp: md[mark] for tech_disp, md in dna_sets.items()
                      if mark in md and md[mark] is not None and md[mark].n_obs > 0}
            if subset:
                combined_violin_frag(subset, mark, args.figdir, args.prefix, tag)
        # NEW: Combined H3K27ac TSSE violin
        combined_violin_tsse_ac(dna_sets, args.figdir, args.prefix, tag)
        # RNA violins
        if rna_sets:
            rna_violin(rna_sets, metric="UMIs", outdir=args.figdir, prefix=args.prefix, tag=tag)
            rna_violin(rna_sets, metric="n_genes", outdir=args.figdir, prefix=args.prefix, tag=tag)

    run_plots(dna_intersect_disp, rna_intersect_disp, tag="intersect")
    run_plots(dna_full_disp,       rna_full_disp,       tag="full")

    # Venns from FULL filtered sets
    for tech_raw in ["InHouse_NanoCNT", "NanoCNT", "NTT", "MTR", "TrES"]:
        tech_disp = display_name(tech_raw)
        dmarks = dna_flt_raw.get(tech_raw, {})
        if "ac" in dmarks and "me3" in dmarks:
            venn2_ac_me3_filtered(dmarks["ac"], dmarks["me3"], tech_disp, args.figdir, args.prefix)
    for tech_raw in ["MTR", "TrES"]:
        tech_disp = display_name(tech_raw)
        dmarks = dna_flt_raw.get(tech_raw, {})
        if ("ac" in dmarks) and ("me3" in dmarks) and (tech_raw in rna_flt_raw):
            venn3_modalities_filtered(dmarks["ac"], dmarks["me3"], rna_flt_raw[tech_raw],
                                      tech_disp, args.figdir, args.prefix)

    # ---------------- Summary ----------------
    print("\n================= SUMMARY: Thresholds & Settings =================")
    print(f"Output figures dir: {os.path.abspath(args.figdir)}")
    print(f"Threads (SnapATAC2): {args.threads}")
    print("\nDNA (H3K27ac):")
    print(f"  min_fragments = {args.ac_min}")
    print(f"  max_fragments = {args.ac_max}")
    print(f"  min_TSSE      = {args.ac_min_tsse}")
    print("\nDNA (H3K27me3):")
    print(f"  min_fragments = {args.me3_min}")
    print(f"  max_fragments = {args.me3_max}")
    print(f"  min_TSSE      = {args.me3_min_tsse}")
    print("\nRNA (for MTR/TrES; 10xMultiome RNA is optional and never intersected with DNA):")
    print(f"  min_UMI       = {args.rna_min_umi}")
    print(f"  max_UMI       = {args.rna_max_umi}")
    print(f"  max_pct_mt    = {args.rna_max_mt}%")
    print(f"  max_pct_hb    = {args.rna_max_hb}%")
    print(f"  max_pct_top50 = {args.rna_max_top50}%")
    print("\nFragment import parameter:")
    print("  sorted_by_barcode = False  (assumed UNSORTED fragments)")
    present_techs = sorted(set(list(dna_map.keys()) + list(rna_files.keys())))
    if present_techs:
        print("\nGenomes by technology (detected):")
        for tech in present_techs:
            genome = "mm10" if tech == "NanoCNT" else "hg38"
            print(f"  - {display_name(tech)}: {genome}")
    print("==================================================================")
    print("Done. Figures written to:", os.path.abspath(args.figdir))
    print("  - *_intersect*.pdf  -> strict pass-in-all-modalities cells")
    print("  - *_full*.pdf       -> per-modality filtered cells")
    print("  - Includes H3K27ac TSSE combined violins for both modes.")

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
    main()
