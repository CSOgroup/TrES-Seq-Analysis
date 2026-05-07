#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FRiP benchmarking across DNA modalities with 3-modality intersections for MTR/TrES.

Workflow per technology:
1) Load DNA (H3K27ac / H3K27me3) (h5ad or fragments) and filter.
2) For MTR & TrES, also load RNA (h5ad), compute QC, and filter.
3) Build the intersection of passing cells:
     - For MTR, TrES:  H3K27ac ∩ H3K27me3 ∩ RNA
     - For others:     H3K27ac ∩ H3K27me3 (or the single available mark)
   (Intersection uses per-tech ID normalization; MTR uses first three colon fields.)
4) Export fragments of the intersected cells to BED (gzip) using snapatac2.ex.export_fragments.
5) Call MACS3 peaks on those BEDs (-f BED; -g hs/mm; mark-specific params).
6) Compute FRiP per cell against those peaks; save .h5ad; plot violins (per mark + combined).

Assumptions:
- Fragment files are UNSORTED by barcode (sorted_by_barcode=False on import).
- Genome: NanoCNT→mm10; all others→hg38.
- MACS3 available on PATH as `macs3`.
- SnapATAC2 export signature (BED) per your version:
    snapatac2.ex.export_fragments(adata, groupby, selections=None, ids=None, ...)

Usage example:
    python frip_pipeline_3modal.py --base-dir /path/to/FullBench --figdir figures_frip

Optional knobs:
    --ac-min 250 --ac-max 15000 --ac-min-tsse 0
    --me3-min 500 --me3-max 15000 --me3-min-tsse 0
    --rna-min-umi 1000 --rna-max-umi 25000 --rna-max-mt 10 --rna-max-hb 5 --rna-max-top50 40
    --threads 24 --force-reimport --no-save-cache
"""

import os
import glob
import argparse
import warnings
import subprocess
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import snapatac2 as snap

# ------------------------- CLI -------------------------

def get_args():
    ap = argparse.ArgumentParser(description="Export intersected fragments (3-modality for MTR/TrES), call peaks, compute/plot FRiP.")
    ap.add_argument("--base-dir", default=".", help="Directory containing InHouse_NanoCNT, MTR, NanoCNT, NTT, TrES")
    ap.add_argument("--threads", type=int, default=24, help="Threads for SnapATAC2")
    ap.add_argument("--figdir", default="figures_frip", help="Output folder for figures/tables")
    ap.add_argument("--prefix", default="FRiPBench", help="Filename prefix for figures")

    # DNA thresholds
    ap.add_argument("--ac-min", type=int, default=250, help="Min fragments per cell for H3K27ac")
    ap.add_argument("--ac-max", type=int, default=15000, help="Max fragments per cell for H3K27ac")
    ap.add_argument("--ac-min-tsse", type=float, default=0.0, help="Min TSSE for H3K27ac")

    ap.add_argument("--me3-min", type=int, default=500, help="Min fragments per cell for H3K27me3")
    ap.add_argument("--me3-max", type=int, default=15000, help="Max fragments per cell for H3K27me3")
    ap.add_argument("--me3-min-tsse", type=float, default=0.0, help="Min TSSE for H3K27me3")

    # RNA thresholds (used for intersection in MTR & TrES)
    ap.add_argument("--rna-min-umi", type=int, default=1000)
    ap.add_argument("--rna-max-umi", type=int, default=25000)
    ap.add_argument("--rna-max-mt", type=float, default=10.0)      # %
    ap.add_argument("--rna-max-hb", type=float, default=5.0)       # %
    ap.add_argument("--rna-max-top50", type=float, default=30.0)   # %

    # Caching controls for fragment imports
    ap.add_argument("--force-reimport", action="store_true",
                    help="Re-import fragments even if a cached .imported.h5ad exists")
    ap.add_argument("--no-save-cache", action="store_true",
                    help="Do not write .imported.h5ad after importing fragments")

    # Reuse controls for downstream steps
    ap.add_argument("--no-reuse-bed",   action="store_true", help="Always re-export BED (default: reuse if exists)")
    ap.add_argument("--no-reuse-peaks", action="store_true", help="Always re-call MACS3 (default: reuse if exists)")
    ap.add_argument("--no-reuse-frip",  action="store_true", help="Always recompute FRiP (default: reuse if exists)")

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
    # Keep colors consistent with your count plots
    return {
        "TrES": "#66c2a5",
        "MTR":  "#3288bd",
        "InHouse_NanoCNT": "#e78ac3",
        "NanoCNT": "#a6d854",
        "NTT": "#ffd92f",
    }

def genome_for_tech(tech: str):
    """NanoCNT is mouse (mm10); all others are human (hg38)."""
    return snap.genome.mm10 if tech == "NanoCNT" else snap.genome.hg38

def macs3_gsize_for_tech(tech: str) -> str:
    """Return MACS3 -g argument."""
    return "mm" if tech == "NanoCNT" else "hs"

def cached_h5ad_path(fragment_path: str) -> str:
    """Derive a sidecar cache path for imported fragments."""
    p = fragment_path
    if p.endswith(".gz"):
        p = p[:-3]
    root, _ext = os.path.splitext(p)
    return root + ".imported.h5ad"

def add_log10_frag(adata: ad.AnnData):
    if "n_fragment" not in adata.obs:
        raise RuntimeError("AnnData missing obs['n_fragment']")
    adata.obs["log10_n_fragment"] = np.log10(adata.obs["n_fragment"] + 1)

def normalize_ids(idx: pd.Index | pd.Series, tech: str) -> pd.Series:
    """Per-tech ID normalization for cross-modality intersections."""
    s = pd.Index(idx) if isinstance(idx, pd.Index) else pd.Index(idx.astype(str))
    if tech == "MTR":
        return pd.Series(s.astype(str).str.split(":").str[:3].str.join(":"), index=s)
    return pd.Series(s.astype(str), index=s)

# ------------------------- Detect inputs -------------------------

def detect_dna_files(base_dir: str) -> Dict[str, Dict[str, Tuple[str, str]]]:
    """
    Map technologies to mark files and types:
      { 'Tech': {'ac': (path, 'h5ad'|'fragments'), 'me3': (path, type)}, ... }

    Layout expectations:
      - TrES:          DNA are top-level *.h5ad
      - InHouse_NanoCNT: prefer raw top-level *.tsv.gz (ignore ./filtered unless needed)
      - MTR:          fragments (*.fragments.bed.gz/.tsv.gz)
      - NanoCNT:      fragments (mouse)
      - NTT:          fragments
    """
    mapping = {}

    def add_entry(tech, mark, path, ftype):
        mapping.setdefault(tech, {})
        mapping[tech][mark] = (path, ftype)

    # TrES DNA: top-level h5ad
    tres_dir = os.path.join(base_dir, "TrES")
    if os.path.isdir(tres_dir):
        ac_h5 = pick_first([os.path.join(tres_dir, "H2RD_Human_H3K27ac.h5ad"),
                            os.path.join(tres_dir, "*H3K27ac*.h5ad")])
        me3_h5 = pick_first([os.path.join(tres_dir, "H2RD_Human_H3K27me3.h5ad"),
                             os.path.join(tres_dir, "*H3K27me3*.h5ad")])
        if ac_h5:  add_entry("TrES", "ac", ac_h5, "h5ad")
        if me3_h5: add_entry("TrES", "me3", me3_h5, "h5ad")

    # InHouse_NanoCNT: prefer raw fragments
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
        # fallback to filtered h5ad if needed
        if "ac" not in mapping.get("InHouse_NanoCNT", {}):
            ac_h5 = pick_first([os.path.join(ih_dir, "filtered", "*H3K27ac*filtered*.h5ad")])
            if ac_h5: add_entry("InHouse_NanoCNT", "ac", ac_h5, "h5ad")
        if "me3" not in mapping.get("InHouse_NanoCNT", {}):
            me3_h5 = pick_first([os.path.join(ih_dir, "filtered", "*H3K27me3*filtered*.h5ad")])
            if me3_h5: add_entry("InHouse_NanoCNT", "me3", me3_h5, "h5ad")

    # MTR: fragments (human)
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

    # NanoCNT: fragments (mouse)
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

    # NTT: fragments (human)
    ntt_dir = os.path.join(base_dir, "NTT")
    if os.path.isdir(ntt_dir):
        ac = pick_first([os.path.join(ntt_dir, "*H3K27ac*.bed.gz"),
                         os.path.join(ntt_dir, "*H3K27ac*.tsv.gz")])
        me3 = pick_first([os.path.join(ntt_dir, "*H3K27me3*.bed.gz"),
                          os.path.join(ntt_dir, "*H3K27me3*.tsv.gz")])
        if ac:  add_entry("NTT", "ac", ac, "fragments")
        if me3: add_entry("NTT", "me3", me3, "fragments")

    return mapping

def detect_rna_files(base_dir: str) -> Dict[str, Tuple[str, str]]:
    """
    RNA inputs used for intersection (MTR, TrES):
      - TrES: top-level *.h5ad (e.g., H2RRNA_Human.h5ad)
      - MTR:  top-level *.h5ad (e.g., scRNA_D0_withCellNames_likeScanpy.h5ad)
    """
    out = {}
    tres_dir = os.path.join(base_dir, "TrES")
    if os.path.isdir(tres_dir):
        tres_raw = pick_first([os.path.join(tres_dir, "H2RRNA_Human.h5ad"),
                               os.path.join(tres_dir, "*RNA*.h5ad"),
                               os.path.join(tres_dir, "*scRNA*.h5ad")])
        if tres_raw and os.path.dirname(tres_raw) == tres_dir:
            out["TrES"] = (tres_raw, "h5ad")

    mtc_dir = os.path.join(base_dir, "MTR")
    if os.path.isdir(mtc_dir):
        mtc_h5 = pick_first([os.path.join(mtc_dir, "scRNA_D0_withCellNames_likeScanpy.h5ad"),
                             os.path.join(mtc_dir, "*scRNA*.h5ad"),
                             os.path.join(mtc_dir, "*.h5ad")])
        if mtc_h5 and os.path.dirname(mtc_h5) == mtc_dir:
            out["MTR"] = (mtc_h5, "h5ad")

    return out

# ------------------------- Load & Filter -------------------------

def load_one_dna(path: str, ftype: str, tech: str, tmpdir: str, threads: int,
                 use_cache: bool = True, force_reimport: bool = False, save_cache: bool = True) -> ad.AnnData:
    g = genome_for_tech(tech)
    if ftype == "h5ad":
        A = sc.read(path)
        if "tsse" not in A.obs.columns:
            try:
                snap.metrics.tsse(A, gene_anno=g, inplace=True)
            except Exception as e:
                warnings.warn(f"Could not compute TSSE for {os.path.basename(path)}: {e}")
        if "n_fragment" not in A.obs.columns:
            raise RuntimeError(f"{path} lacks obs['n_fragment'].")
        add_log10_frag(A)
        return A

    elif ftype == "fragments":
        cache_path = cached_h5ad_path(path)
        if use_cache and not force_reimport and os.path.exists(cache_path):
            print(f"  -> Using cached AnnData: {cache_path}")
            A = sc.read(cache_path)
            if "n_fragment" not in A.obs.columns:
                warnings.warn(f"Cached {cache_path} missing 'n_fragment'; re-importing.")
            else:
                if "tsse" not in A.obs.columns:
                    try:
                        snap.metrics.tsse(A, gene_anno=g, inplace=True)
                    except Exception as e:
                        warnings.warn(f"TSSE on cached failed: {e}")
                add_log10_frag(A)
                return A

        print(f"  -> Importing fragments (unsorted) from {path}")
        A = snap.pp.import_fragments(
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
            snap.metrics.tsse(A, gene_anno=g, inplace=True)
        except Exception as e:
            warnings.warn(f"TSSE computation failed for {os.path.basename(path)}: {e}")
        add_log10_frag(A)

        if save_cache:
            try:
                A.write(cache_path)
                print(f"  -> Cached imported AnnData at {cache_path}")
            except Exception as e:
                warnings.warn(f"Could not write cache {cache_path}: {e}")

        return A

    else:
        raise ValueError("ftype must be 'h5ad' or 'fragments'.")

def filter_dna(A: ad.AnnData, min_frags: int, max_frags: int, min_tsse: float) -> ad.AnnData:
    mask = (A.obs["n_fragment"] >= min_frags) & (A.obs["n_fragment"] <= max_frags)
    if "tsse" in A.obs.columns and min_tsse > 0:
        mask = mask & (A.obs["tsse"] >= min_tsse)
    Af = A[mask].copy()
    add_log10_frag(Af)
    return Af

def load_one_rna(path: str, kind: str) -> ad.AnnData:
    if kind != "h5ad":
        raise ValueError("RNA loader expects h5ad for MTR/TrES.")
    return sc.read(path)

def compute_rna_qc(A: ad.AnnData):
    gs = A.var_names.astype(str)
    gs_u = gs.str.upper()
    A.var["mt"] = gs_u.str.startswith("MT-")
    A.var["hb"] = gs_u.str.match(r"^HB[AB]")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt","hb"], percent_top=[50], inplace=True)
    A.obs.rename(columns={
        "total_counts": "UMIs",
        "n_genes_by_counts": "n_genes",
        "pct_counts_mt": "pct_mt",
        "pct_counts_hb": "pct_hb",
        "pct_counts_in_top_50_genes": "pct_top50",
    }, inplace=True)

def filter_rna(A: ad.AnnData, minUMI, maxUMI, max_mt, max_hb, max_top50) -> ad.AnnData:
    qc = A.obs
    keep = (
        (qc["UMIs"] >= minUMI) & (qc["UMIs"] <= maxUMI) &
        (qc["pct_mt"] <= max_mt) & (qc["pct_hb"] <= max_hb) &
        (qc["pct_top50"] <= max_top50)
    )
    return A[keep].copy()

# ------------------------- Export, MACS3, FRiP -------------------------

def export_bed_for_cells(A: ad.AnnData, out_dir: str, base: str, reuse: bool = True) -> str:
    """
    Export fragments for all rows in A to a single BED (.bed.gz).
    Uses snapatac2.ex.export_fragments(groupby=<const>).
    Returns the path to the written BED file.
    """
    ensure_dir(out_dir)
    expected = os.path.join(out_dir, f"{base}_all.bed.gz")
    if reuse and os.path.exists(expected):
        print(f"  -> Reusing BED: {expected}")
        return expected

    group_col = "__grp__"
    A_obs = A.obs.copy()
    A.obs[group_col] = "all"
    files = snap.ex.export_fragments(
        A,
        groupby=group_col,
        out_dir=out_dir,
        prefix=f"{base}_",
        suffix=".bed.gz",
        compression="gzip",
        compression_level=1,
    )
    A.obs = A_obs
    out_path = files.get("all", expected)
    if out_path != expected and os.path.exists(out_path):
        # Standardize name for future reuse
        try:
            os.replace(out_path, expected)
            out_path = expected
        except Exception:
            pass
    if not os.path.exists(out_path):
        raise FileNotFoundError(f"Could not find exported BED for {base} in {out_dir}")
    return out_path


def call_macs3(bed_path: str, tech: str, mark: str, out_dir: str, name: str, reuse: bool = True) -> str:
    """Call peaks with MACS3 on BED fragments.

    - H3K27ac: narrow peaks, FRAG format, q=0.01
    - H3K27me3: broad peaks (broadPeak used for FRiP), FRAG format, q=0.01
    """
    ensure_dir(out_dir)
    narrow = os.path.join(out_dir, f"{name}_peaks.narrowPeak")
    broad  = os.path.join(out_dir, f"{name}_peaks.broadPeak")

    # Reuse logic: ac → narrowPeak, me3 → broadPeak
    if reuse:
        if mark == "ac" and os.path.exists(narrow):
            print(f"  -> Reusing peaks: {narrow}")
            return narrow
        if mark == "me3" and os.path.exists(broad):
            print(f"  -> Reusing peaks: {broad}")
            return broad

    gsize = macs3_gsize_for_tech(tech)

    # Base MACS3 command (shared between marks)
    cmd = [
        "macs3", "callpeak",
        "-t", bed_path,
        "-f", "FRAG",
        "-g", gsize,
        "--keep-dup", "0",
        "--nolambda",
        "--nomodel",
        "--name", name,
        "--outdir", out_dir,
        "--seed", "99",
        "--verbose", "1",
        "-q", "0.001",
    ]

    # H3K27me3: broad peaks
    if mark == "me3":
        cmd.append("--broad")

    print("  -> MACS3:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # Return the file used for FRiP
    if mark == "me3":
        if os.path.exists(broad):
            return broad
        # Fallback in case MACS3 naming changes
        if os.path.exists(narrow):
            return narrow
        raise FileNotFoundError(f"No broadPeak (or narrowPeak) found for {name} in {out_dir}")
    else:  # H3K27ac
        if os.path.exists(narrow):
            return narrow
        # Fallback if only broadPeak exists for some reason
        if os.path.exists(broad):
            return broad
        raise FileNotFoundError(f"No narrowPeak (or broadPeak) found for {name} in {out_dir}")




def compute_frip(A: ad.AnnData, peaks_file: str):
    """Annotate FRiP per cell (normalized) using the given peaks file."""
    snap.metrics.frip(
        A,
        regions={"FRiP": peaks_file},
        normalized=True,
        inplace=True,
    )
    if "FRiP" not in A.obs.columns:
        for c in ["frip", "FRIP", "FRiP_score", "frip_score"]:
            if c in A.obs.columns:
                A.obs.rename(columns={c: "FRiP"}, inplace=True)
                break

def maybe_load_frip(out_h5: str, reuse: bool) -> Optional[ad.AnnData]:
    if reuse and os.path.exists(out_h5):
        A = sc.read(out_h5)
        # Normalize FRiP column name if needed
        for c in list(A.obs.columns):
            if c.lower() == "frip" and c != "FRiP":
                A.obs.rename(columns={c: "FRiP"}, inplace=True)
        if "FRiP" in A.obs.columns:
            print(f"  -> Reusing FRiP AnnData: {out_h5}")
            return A
    return None

# ------------------------- Plotting -------------------------

def frip_df_stack(records: Dict[str, Dict[str, ad.AnnData]]) -> pd.DataFrame:
    """
    records: {tech: {'ac': AnnData (intersect), 'me3': AnnData (intersect)}}
    -> tidy DF with Technology, Mark, FRiP and FRiP_pct
    """
    rows = []
    for tech, md in records.items():
        for mark, A in md.items():
            if A is None or "FRiP" not in A.obs:
                continue
            frip = pd.to_numeric(A.obs["FRiP"], errors="coerce").clip(lower=0, upper=1)
            for v in frip.dropna().values:
                rows.append({
                    "Technology": tech,
                    "Mark": "H3K27ac" if mark == "ac" else "H3K27me3",
                    "FRiP": float(v),
                    "FRiP_pct": float(v)*100.0
                })
    return pd.DataFrame(rows)

def plot_frip_violins(df: pd.DataFrame, figdir: str, prefix: str):
    if df.empty:
        print("No FRiP data to plot.")
        return
    pal = tech_palette()
    ensure_dir(figdir)

    # Per-mark plots
    for mark in ["H3K27ac", "H3K27me3"]:
        sub = df[df["Mark"] == mark]
        if sub.empty:
            continue
        order = [t for t in ["TrES","MTR","InHouse_NanoCNT","NanoCNT","NTT"] if t in sub["Technology"].unique()]
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.violinplot(x="Technology", y="FRiP_pct", data=sub, order=order, palette=pal, inner=None, ax=ax)
        sns.boxplot(x="Technology", y="FRiP_pct", data=sub, order=order, ax=ax, width=0.2,
                    showcaps=True,
                    boxprops={'facecolor': 'none', 'linewidth': 1.5},
                    whiskerprops={'linewidth': 1.5},
                    medianprops={'color': 'black'},
                    flierprops={'marker': '', 'markersize': 0})
        ax.set_xlabel("")
        ax.set_ylabel("FRiP (%)")
        ax.set_title(f"{mark} – FRiP (intersect cells)")
        ax.set_ylim(0, 100)
        plt.tight_layout()
        out = os.path.join(figdir, f"{prefix}_{mark}_FRiP_Violin.pdf")
        fig.savefig(out); plt.close(fig); print(out)

    # Combined (marks on x, hue=technology) – keep hue_order identical for both plots
    order_mark = [m for m in ["H3K27me3", "H3K27ac"] if m in df["Mark"].unique()]
    order_tech = [t for t in ["TrES","MTR","InHouse_NanoCNT","NanoCNT","NTT"] if t in df["Technology"].unique()]
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.violinplot(
        x="Mark", y="FRiP_pct", hue="Technology",
        data=df, order=order_mark, hue_order=order_tech,
        palette=pal, inner=None, dodge=True, ax=ax,
    )
    sns.boxplot(
        x="Mark", y="FRiP_pct", hue="Technology",
        data=df, order=order_mark, hue_order=order_tech,
        dodge=True, width=0.2, ax=ax,
        showcaps=True,
        boxprops={'facecolor': 'none', 'linewidth': 1.5},
        whiskerprops={'linewidth': 1.5},
        medianprops={'color': 'black'},
        flierprops={'marker': '', 'markersize': 0},
    )
    ax.set_ylabel("FRiP (%)")
    ax.set_xlabel("")
    ax.set_title("FRiP by Mark (intersect cells)")
    ax.set_ylim(0, 100)

    # Deduplicate legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles[:len(order_tech)], labels[:len(order_tech)],
                  frameon=False, loc="upper right")

    plt.tight_layout()
    out = os.path.join(figdir, f"{prefix}_BothMarks_FRiP_Violin.pdf")
    fig.savefig(out); plt.close(fig); print(out)

# ------------------------- Main -------------------------

def main():
    args = get_args()
    ensure_dir(args.figdir)
    tmpdir = os.path.join(args.base_dir, "_tmp_snap")
    ensure_dir(tmpdir)

    reuse_bed   = not args.no_reuse_bed
    reuse_peaks = not args.no_reuse_peaks
    reuse_frip  = not args.no_reuse_frip

    dna_map = detect_dna_files(args.base_dir)
    rna_map = detect_rna_files(args.base_dir)  # only MTR/TrES present

    # Load RNA for intersection where needed (MTR/TrES)
    rna_flt: Dict[str, ad.AnnData] = {}
    for tech in ["MTR", "TrES"]:
        if tech in rna_map:
            rna_path, kind = rna_map[tech]
            print(f"[RNA] Loading {tech} from {rna_path} ({kind})")
            R = load_one_rna(rna_path, kind)
            compute_rna_qc(R)
            Rf = filter_rna(R, args.rna_min_umi, args.rna_max_umi, args.rna_max_mt, args.rna_max_hb, args.rna_max_top50)
            Rf.obs["__norm_id"] = normalize_ids(Rf.obs_names, tech).values
            rna_flt[tech] = Rf
            print(f"  -> {tech} RNA passing: {Rf.n_obs}")

    # Load & filter DNA
    dna_flt: Dict[str, Dict[str, ad.AnnData]] = {}
    dna_raw: Dict[str, Dict[str, ad.AnnData]] = {}
    for tech, marks in dna_map.items():
        dna_flt.setdefault(tech, {})
        dna_raw.setdefault(tech, {})
        for mark, (path, ftype) in marks.items():
            print(f"[DNA] Loading {tech} {mark} from {path} ({ftype})")
            A = load_one_dna(
                path, ftype, tech, tmpdir, args.threads,
                use_cache=not args.force_reimport,
                force_reimport=args.force_reimport,
                save_cache=not args.no_save_cache
            )
            A.obs["__norm_id"] = normalize_ids(A.obs_names, tech).values
            dna_raw[tech][mark] = A
            if mark == "ac":
                Af = filter_dna(A, args.ac_min, args.ac_max, args.ac_min_tsse)
            else:
                Af = filter_dna(A, args.me3_min, args.me3_max, args.me3_min_tsse)
            Af.obs["__norm_id"] = normalize_ids(Af.obs_names, tech).values
            dna_flt[tech][mark] = Af
            print(f"  -> {tech} {mark}: {Af.n_obs} cells after filters")

    # Intersections (3-modality for MTR/TrES)
    bed_dir   = os.path.join(args.figdir, "bed_for_peaks")
    peaks_dir = os.path.join(args.figdir, "macs3_peaks")
    frip_dir  = os.path.join(args.figdir, "h5ad_frip")
    ensure_dir(bed_dir); ensure_dir(peaks_dir); ensure_dir(frip_dir)

    adatas_intersect: Dict[str, Dict[str, ad.AnnData]] = {}
    summary_rows = []

    for tech, md in dna_flt.items():
        has_ac  = "ac" in md
        has_me3 = "me3" in md
        if not has_ac and not has_me3:
            continue

        # Base intersection = DNA marks
        if has_ac and has_me3:
            common = set(md["ac"].obs["__norm_id"]) & set(md["me3"].obs["__norm_id"])
        else:
            mk = "ac" if has_ac else "me3"
            common = set(md[mk].obs["__norm_id"])

        # For MTR/TrES, further require RNA pass
        if tech in rna_flt:
            common = common & set(rna_flt[tech].obs["__norm_id"])

        print(f"[{tech}] intersect cells (required modalities): {len(common)}")

        adatas_intersect.setdefault(tech, {})
        for mark in ["ac", "me3"]:
            if mark not in md or len(common) == 0:
                adatas_intersect[tech][mark] = None
                continue

            Af = md[mark]
            keep_mask = Af.obs["__norm_id"].isin(common)
            A_sub = Af[keep_mask].copy()
            adatas_intersect[tech][mark] = A_sub

            base = f"{tech}_H3K27{mark}_intersect"
            bed_path = export_bed_for_cells(A_sub, bed_dir, base, reuse=reuse_bed)
            print(f"  -> BED: {bed_path}")

            peaks_path = call_macs3(bed_path, tech, mark, peaks_dir, base, reuse=reuse_peaks)
            print(f"  -> Peaks: {peaks_path}")

            out_h5 = os.path.join(frip_dir, f"{base}_FRiP.h5ad")
            A_cached = maybe_load_frip(out_h5, reuse=reuse_frip)
            if A_cached is None:
                compute_frip(A_sub, peaks_path)
                A_sub.write(out_h5)
                print(f"  -> Saved FRiP AnnData: {out_h5}")
            else:
                A_sub = A_cached  # reuse computed FRiP
                adatas_intersect[tech][mark] = A_sub

            med = float(np.nanmedian(pd.to_numeric(A_sub.obs.get("FRiP", np.nan), errors="coerce")))
            summary_rows.append({
                "Technology": tech,
                "Mark": "H3K27ac" if mark == "ac" else "H3K27me3",
                "n_cells": A_sub.n_obs,
                "FRiP_median": med,
            })

    # Save summary CSV
    if summary_rows:
        df_sum = pd.DataFrame(summary_rows).sort_values(["Mark", "Technology"])
        out_csv = os.path.join(args.figdir, f"{args.prefix}_FRiP_summary.csv")
        df_sum.to_csv(out_csv, index=False)
        print(out_csv)

    # Plots
    frip_df = frip_df_stack(adatas_intersect)
    plot_frip_violins(frip_df, args.figdir, args.prefix)

    # Thresholds echo
    print("\n================= SUMMARY: Thresholds & Settings =================")
    print(f"Output dir: {os.path.abspath(args.figdir)}")
    print(f"Threads (SnapATAC2): {args.threads}")
    print("\nDNA (H3K27ac):")
    print(f"  min_fragments = {args.ac_min}")
    print(f"  max_fragments = {args.ac_max}")
    print(f"  min_TSSE      = {args.ac_min_tsse}")
    print("\nDNA (H3K27me3):")
    print(f"  min_fragments = {args.me3_min}")
    print(f"  max_fragments = {args.me3_max}")
    print(f"  min_TSSE      = {args.me3_min_tsse}")
    print("\nRNA (used for intersection in MTR/TrES):")
    print(f"  min_UMI       = {args.rna_min_umi}")
    print(f"  max_UMI       = {args.rna_max_umi}")
    print(f"  max_pct_mt    = {args.rna_max_mt}%")
    print(f"  max_pct_hb    = {args.rna_max_hb}%")
    print(f"  max_pct_top50 = {args.rna_max_top50}%")
    print("\nFragment import parameter:")
    print("  sorted_by_barcode = False  (assumed UNSORTED fragments)")
    # Genomes by tech
    present_techs = sorted(dna_map.keys())
    if present_techs:
        print("\nGenomes by technology (detected):")
        for tech in present_techs:
            genome = "mm10" if tech == "NanoCNT" else "hg38"
            print(f"  - {tech}: {genome}")
    print("==================================================================")
    print("Done.")

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
    main()
