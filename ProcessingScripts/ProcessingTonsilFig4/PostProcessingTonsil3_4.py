#!/usr/bin/env python3
import os
import re
import sys
import csv
import glob
import math
import argparse
import itertools
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import snapatac2 as snap
import scanpy as sc
import anndata as ad

from upsetplot import plot, from_indicators
from matplotlib_venn import venn3
from scipy.io import mmread
from matplotlib.ticker import FixedLocator, FixedFormatter


# =============================================================================
# Utilities
# =============================================================================

def ensure_dir(p: str):
    if p and (not os.path.exists(p)):
        os.makedirs(p, exist_ok=True)

def collapse_tp_run_sample(sample: str) -> str:
    """
    Collapse TP run-specific names:
      Sc_TP3D_1_S5 -> Sc_TP3D
      Sc_TP4D_2_S10 -> Sc_TP4D
      Sc_TP3r_3_S8 -> Sc_TP3r
    Only for TP3/TP4 patterns; otherwise return unchanged.
    """
    m = re.match(r"^(Sc_TP[34][Dr])_\d+_S\d+$", sample)
    if m:
        return m.group(1)
    return sample

def read_mo_map(mo_map_path: str):
    """
    Supports BOTH:
      (A) 3-col (legacy): sample  mark  mo_bc
      (B) 4-col (tonsil): sample  sb_group  mark  mo_bc

    For 4-col tonsil, we collapse run-specific samples and use a group key:
      key = {collapsed_sample}_{sb_group}
      e.g. Sc_TP3D_TP3_A

    Returns: dict key -> sorted unique marks
    """
    sample_to_marks = {}
    with open(mo_map_path, "r") as f:
        for line in f:
            line = line.strip()
            if (not line) or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 3:
                continue

            if len(parts) >= 4:
                # tonsil: sample sb_group mark mo_bc
                s_raw, sb_group, mark = parts[0], parts[1], parts[2]
                s = collapse_tp_run_sample(s_raw)
                key = f"{s}_{sb_group}"
            else:
                # legacy: sample mark mo_bc
                key, mark = parts[0], parts[1]

            sample_to_marks.setdefault(key, set()).add(mark)

    return {k: sorted(list(v)) for k, v in sample_to_marks.items()}

def detect_dna_bam(workdir: str, dna_sample: str, mark: str) -> str:
    """
    Prefer NoDup, else MarkedDup.
    Expected filenames:
      {workdir}/{dna_sample}_{mark}_NoDup.bam
      {workdir}/{dna_sample}_{mark}_MarkedDup.bam
    Example:
      Sc_TP3D_TP3_A_H3K27me3_NoDup.bam
    """
    p1 = os.path.join(workdir, f"{dna_sample}_{mark}_NoDup.bam")
    if os.path.exists(p1) and os.path.getsize(p1) > 0:
        return p1
    p2 = os.path.join(workdir, f"{dna_sample}_{mark}_MarkedDup.bam")
    if os.path.exists(p2) and os.path.getsize(p2) > 0:
        return p2
    return ""

def normalize_rna_barcode(bc: str) -> str:
    """
    STARsolo barcodes often have suffix like '-1'. DNA CBs usually do not.
    Keep part before '-' if it looks like a suffix.
    """
    if "-" in bc:
        left = bc.split("-", 1)[0]
        if 6 <= len(left) <= 128:
            return left
    return bc

def infer_rna_samples(workdir: str):
    """
    RNA samples inferred from *.filtered_cells.bam in workdir.
    """
    out = []
    for p in glob.glob(os.path.join(workdir, "*.filtered_cells.bam")):
        bn = os.path.basename(p)
        if bn.endswith(".filtered_cells.bam"):
            out.append(bn[:-len(".filtered_cells.bam")])
    out.sort()
    return out

def infer_pairing_tonsil(dna_sample: str, rna_samples: list) -> str:
    """
    Tonsil TP pairing:
      Sc_TP3D_TP3_A -> Sc_TP3r_TP3_A
      Sc_TP4D_TP4_B -> Sc_TP4r_TP4_B
    """
    # Replace the FIRST occurrence of TP{3,4}D with TP{3,4}r
    target = re.sub(r"(Sc_TP[34])D", r"\1r", dna_sample, count=1)
    if target in rna_samples:
        return target

    # If RNA samples have extra suffixes, try normalized matching
    target_n = re.sub(r"_S\d+$", "", target)
    for r in rna_samples:
        rn = re.sub(r"_S\d+$", "", r)
        if rn == target_n or rn.startswith(target_n):
            return r
    return ""

def infer_pairing_legacy(dna_sample: str, rna_samples: list) -> str:
    """
    Heuristic pairing:
      - Replace '...D' -> '...r' in FIRST token after 'Sc_' (before first underscore)
      - Ignore trailing '_S<digits>' when matching
    """
    body = dna_sample[3:] if dna_sample.startswith("Sc_") else dna_sample
    body_nos = re.sub(r"_S\d+$", "", body)

    parts = body_nos.split("_")
    if parts:
        head = parts[0]
        if head.endswith("D"):
            head = head[:-1] + "r"
        parts[0] = head

    target = "Sc_" + "_".join(parts)

    rna_norm = [(r, re.sub(r"_S\d+$", "", r)) for r in rna_samples]

    for r, rn in rna_norm:
        if rn == target:
            return r
    for r, rn in rna_norm:
        if rn.startswith(target):
            return r
    return ""


# =============================================================================
# Inline SB rewrite: DT(3nt) -> RNA(4nt)
# =============================================================================

IDX_DT = {
    "GAT": "1",
    "AGT": "2",
    "TCA": "3",
    "ACG": "4",
    "GTC": "5",
    "CTG": "6",
    "GAC": "7",
    "AGA": "8",
}

IDX_RNA = {
    "GCAT": "1",
    "CGAT": "2",
    "GACT": "3",
    "AGCT": "4",
    "CAGT": "5",
    "ACGT": "6",
    "GCTA": "7",
    "CGTA": "8",
}

NUM_TO_RNA = {v: k for k, v in IDX_RNA.items()}
DT_TO_RNA = {dt: NUM_TO_RNA[num] for dt, num in IDX_DT.items()}

INLINE_OFFSET = 1  # injected base at position 0; SB begins at position 1
LEN_DT = 3
LEN_RNA = 4

def _uniquify_names(names):
    s = pd.Series([str(x) for x in names], dtype=object)
    if s.duplicated().any():
        dup = s.duplicated(keep=False)
        s.loc[dup] = s.loc[dup] + "_" + s.loc[dup].groupby(s.loc[dup]).cumcount().astype(str)
    return pd.Index(s.astype(str).to_numpy(dtype=object))

def rewrite_inline_sb_to_rna_sb(
    adata: ad.AnnData,
    assay: str,
    inline_offset: int = INLINE_OFFSET,
    min_recognized_frac: float = 0.95,
    strict: bool = True,
) -> ad.AnnData:
    """
    Rewrite obs_names so that the inline SB segment is always the 4-nt RNA SB.

    Input barcode structure (expected):
      - DNA DT:  <BASE><DT3><rest...>   (DT3 at offset=1, length=3)
      - RNA:     <BASE><RNA4><rest...>  (RNA4 at offset=1, length=4)

    Output barcode structure:
      - DNA rewritten: <BASE><RNA4><rest...>   (DT3 replaced by RNA4, injected base kept)
      - RNA: unchanged (but validated)
    """
    if assay not in ("dna_dt", "rna"):
        raise ValueError("assay must be 'dna_dt' or 'rna'")

    orig = pd.Series([str(x) for x in adata.obs_names], index=adata.obs_names, dtype=object)

    if assay == "dna_dt":
        tag_len = LEN_DT
        tag = orig.str.slice(inline_offset, inline_offset + tag_len)
        rna_tag = tag.map(DT_TO_RNA)
        recognized = rna_tag.notna()
        frac = float(recognized.mean()) if len(tag) else 0.0

        if frac < min_recognized_frac:
            bad_n = int((~recognized).sum())
            ex = None
            if bad_n > 0:
                ex_i = (~recognized).to_numpy().nonzero()[0][0]
                ex = f"{orig.iloc[ex_i]} (DT_tag='{tag.iloc[ex_i]}')"
            msg = f"[rewrite_inline_sb_to_rna_sb] DNA: recognized_frac={frac:.3f} ({len(tag)-bad_n}/{len(tag)}). Example bad: {ex}"
            if strict:
                raise ValueError(msg)
            sys.stderr.write("WARNING: " + msg + " -> leaving DNA obs_names unchanged.\n")
            adata.obs["barcode_raw"] = orig.values
            adata.obs["inline_sb_dt"] = tag.values
            adata.obs["inline_sb_rna"] = rna_tag.astype(object).values
            adata.obs["barcode_rewritten"] = orig.values
            return adata

        base_part = orig.str.slice(0, inline_offset)
        rest_part = orig.str.slice(inline_offset + tag_len, None)
        new_names = (base_part + rna_tag + rest_part).astype(str)

        adata.obs["barcode_raw"] = orig.values
        adata.obs["inline_sb_dt"] = tag.values
        adata.obs["inline_sb_rna"] = rna_tag.values
        adata.obs["barcode_rewritten"] = new_names.values

        adata.obs_names = _uniquify_names(new_names.values)
        return adata

    else:
        tag_len = LEN_RNA
        tag = orig.str.slice(inline_offset, inline_offset + tag_len)
        recognized = tag.isin(IDX_RNA.keys())
        frac = float(recognized.mean()) if len(tag) else 0.0

        if frac < min_recognized_frac:
            bad_n = int((~recognized).sum())
            ex = None
            if bad_n > 0:
                ex_i = (~recognized).to_numpy().nonzero()[0][0]
                ex = f"{orig.iloc[ex_i]} (RNA_tag='{tag.iloc[ex_i]}')"
            msg = f"[rewrite_inline_sb_to_rna_sb] RNA: recognized_frac={frac:.3f} ({len(tag)-bad_n}/{len(tag)}). Example bad: {ex}"
            if strict:
                raise ValueError(msg)
            sys.stderr.write("WARNING: " + msg + " -> leaving RNA obs_names unchanged.\n")

        adata.obs["barcode_raw"] = orig.values
        adata.obs["inline_sb_rna"] = tag.values
        adata.obs["barcode_rewritten"] = orig.values
        return adata

def swap_injected_base_CG_inplace(adata: ad.AnnData):
    """
    Swap injected base C<->G at position 0 of obs_names.
    (Only used for PDTD RNA, kept for backward compatibility.)
    """
    orig = pd.Series([str(x) for x in adata.obs_names], index=adata.obs_names, dtype=object)

    def _swap(s: str) -> str:
        if not s:
            return s
        b = s[0]
        if b == "C":
            return "G" + s[1:]
        if b == "G":
            return "C" + s[1:]
        return s

    new = orig.map(_swap).astype(str)
    adata.obs["barcode_injected_base_swapped"] = True
    adata.obs["barcode_before_base_swap"] = orig.values
    adata.obs_names = _uniquify_names(new.values)


# =============================================================================
# STARsolo reader (raw)
# =============================================================================

def read_star_solo_filtered(workdir: str, rna_sample: str) -> ad.AnnData:
    """
    Reads:
      {workdir}/{rna_sample}.Solo.outGeneFull/raw/
        matrix.mtx
        barcodes.tsv
        features.tsv
    """
    base = os.path.join(workdir, f"{rna_sample}.Solo.outGeneFull", "filtered")
    mtx = os.path.join(base, "matrix.mtx")
    bcs = os.path.join(base, "barcodes.tsv")
    feats = os.path.join(base, "features.tsv")

    if not (os.path.exists(mtx) and os.path.exists(bcs) and os.path.exists(feats)):
        raise FileNotFoundError(
            f"Missing STARsolo output for {rna_sample}. "
            f"Expected {mtx}, {bcs}, {feats}"
        )

    X = mmread(mtx).tocsr()

    barcodes = []
    with open(bcs, "r") as f:
        for line in f:
            bc = line.strip()
            if bc:
                barcodes.append(normalize_rna_barcode(bc))

    feat_df = pd.read_csv(feats, sep="\t", header=None)
    if feat_df.shape[1] >= 2:
        gene_names = feat_df.iloc[:, 1].astype(str).tolist()
        gene_ids = feat_df.iloc[:, 0].astype(str).tolist()
        var_names = [gn if (gn and gn != "nan") else gid for gn, gid in zip(gene_names, gene_ids)]
    else:
        var_names = feat_df.iloc[:, 0].astype(str).tolist()

    if X.shape[1] == len(barcodes) and X.shape[0] == len(var_names):
        X = X.T.tocsr()
    elif not (X.shape[0] == len(barcodes) and X.shape[1] == len(var_names)):
        raise ValueError(
            f"Matrix shape {X.shape} does not match barcodes ({len(barcodes)}) "
            f"and features ({len(var_names)})."
        )

    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index(barcodes, dtype=object)
    adata.var_names = pd.Index(var_names, dtype=object)
    adata.var_names_make_unique()
    return adata


# =============================================================================
# Nice tick helpers (your version)
# =============================================================================

def _tick_candidates_12510(min_raw: float, max_raw: float):
    if not np.isfinite(min_raw) or not np.isfinite(max_raw) or max_raw <= 0:
        return []
    if max_raw < min_raw:
        min_raw, max_raw = max_raw, min_raw

    lo = max(1.0, min_raw)
    hi = max_raw
    k0 = int(math.floor(math.log10(lo)))
    k1 = int(math.ceil(math.log10(hi)))

    out = []
    for k in range(k0 - 1, k1 + 1):
        base = 10 ** k
        for m in range(1, 11):
            t = m * base
            if t >= min_raw and t <= max_raw:
                out.append(int(round(t)))
    return sorted(set(out))

def smart_count_ticks(min_raw: float, max_raw: float, max_ticks: int = 24):
    if not np.isfinite(min_raw) or not np.isfinite(max_raw):
        return []
    if max_raw < min_raw:
        min_raw, max_raw = max_raw, min_raw
    if max_raw <= 0:
        return []

    min_raw = max(0.0, float(min_raw))
    max_raw = float(max_raw)

    candidates = _tick_candidates_12510(min_raw, max_raw)
    if not candidates:
        return []
    if len(candidates) <= max_ticks:
        return candidates

    low_region_max = max(2000.0, min_raw * 50.0)
    low = [t for t in candidates if t <= low_region_max]
    high = [t for t in candidates if t > low_region_max]

    if high:
        high_budget = max(6, max_ticks // 3)
        low_budget = max_ticks - high_budget
        low_budget = max(8, low_budget)
        high_budget = max_ticks - low_budget
    else:
        low_budget = max_ticks
        high_budget = 0

    if len(low) > low_budget:
        very_low = [t for t in low if t <= 200]
        rest_low = [t for t in low if t > 200]
        out_low = []
        out_low.extend(very_low[: min(len(very_low), low_budget // 2)])

        remaining = low_budget - len(out_low)
        if remaining > 0 and rest_low:
            idx = np.unique(np.round(np.linspace(0, len(rest_low) - 1, remaining)).astype(int))
            out_low.extend([rest_low[i] for i in idx])

        low_keep = sorted(set(out_low))
    else:
        low_keep = low

    high_keep = []
    if high_budget > 0 and high:
        if len(high) <= high_budget:
            high_keep = high
        else:
            idx = np.unique(np.round(np.linspace(0, len(high) - 1, high_budget)).astype(int))
            high_keep = [high[i] for i in idx]
            if high[-1] not in high_keep:
                high_keep[-1] = high[-1]

    out = sorted(set(low_keep + high_keep))
    if len(out) > max_ticks:
        out = out[: max_ticks - 1] + [out[-1]]
        out = sorted(set(out))
    return out

def set_log_transformed_axis_ticks(
    ax,
    transform: str,
    raw_values=None,
    max_ticks: int = 24,
    drop_zero: bool = True,
    max_quantile: float = 0.99,
):
    y0, y1 = ax.get_ylim()
    if not (np.isfinite(y0) and np.isfinite(y1)):
        return
    lo, hi = (y0, y1) if y0 <= y1 else (y1, y0)

    if transform == "log10p1":
        raw_lo_axis = max(0.0, (10.0 ** lo) - 1.0)
        raw_hi_axis = max(0.0, (10.0 ** hi) - 1.0)
    elif transform == "log1p":
        raw_lo_axis = max(0.0, math.expm1(lo))
        raw_hi_axis = max(0.0, math.expm1(hi))
    else:
        raise ValueError("transform must be 'log10p1' or 'log1p'")

    raw_min_data = None
    raw_max_cap = None
    if raw_values is not None:
        rv = np.asarray(raw_values, dtype=float)
        rv = rv[np.isfinite(rv)]
        if rv.size > 0:
            pos = rv[rv > 0]
            if pos.size > 0:
                raw_min_data = float(pos.min())
            raw_max_cap = float(np.quantile(rv, max_quantile))

    tick_min = raw_lo_axis
    if raw_min_data is not None:
        tick_min = max(tick_min, raw_min_data)

    tick_max = raw_hi_axis
    if raw_max_cap is not None and raw_max_cap > 0:
        tick_max = min(tick_max, raw_max_cap)
    if tick_max <= tick_min:
        tick_max = raw_hi_axis

    raw_ticks = smart_count_ticks(tick_min, tick_max, max_ticks=max_ticks)
    if drop_zero:
        raw_ticks = [t for t in raw_ticks if t != 0]
    if not raw_ticks:
        return

    if transform == "log10p1":
        pos = [math.log10(t + 1.0) for t in raw_ticks]
    else:
        pos = [math.log1p(t) for t in raw_ticks]

    lbl = [f"{int(t):,}" for t in raw_ticks]
    ax.yaxis.set_major_locator(FixedLocator(pos))
    ax.yaxis.set_major_formatter(FixedFormatter(lbl))


# =============================================================================
# Script-B-style violin helpers
# =============================================================================

def _backtransform(vals, transform: str):
    vals = np.asarray(vals, dtype=float)
    if transform == "log10p1":
        return (10.0 ** vals) - 1.0
    elif transform == "log1p":
        return np.expm1(vals)
    else:
        raise ValueError("transform must be 'log10p1' or 'log1p'")

def bstyle_single_violin(
    adata: ad.AnnData,
    obs_key: str,
    out_pdf: str,
    title: str,
    sample_label: str,
    y_label: str,
    transform: str,
    force_plot: bool,
):
    """
    Script-B-like violin:
      - seaborn violin
      - boxplot overlay
      - median dot
      - raw median text label
      - axis tick labels back-transformed to raw counts
    """
    if (not force_plot) and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 0:
        return

    if obs_key not in adata.obs.columns:
        sys.stderr.write(f"WARNING: missing {obs_key} for violin plot: {out_pdf}\n")
        return

    vals = pd.to_numeric(adata.obs[obs_key], errors="coerce").to_numpy()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        sys.stderr.write(f"WARNING: no finite values in {obs_key} for violin plot: {out_pdf}\n")
        return

    df = pd.DataFrame(
        {
            "sample_label": [f"{sample_label}\n(n={adata.n_obs:,})"] * vals.size,
            "metric": vals,
        }
    )

    plt.figure(figsize=(4, 4))
    sns.violinplot(
        data=df,
        x="sample_label",
        y="metric",
        inner=None,
        cut=0,
        linewidth=0.5,
        color="lightgray",
    )
    sns.boxplot(
        data=df,
        x="sample_label",
        y="metric",
        whis=1.5,
        width=0.18,
        fliersize=0,
        boxprops=dict(alpha=0.6),
    )

    ax = plt.gca()

    med_plot = float(np.median(vals))
    med_raw = float(_backtransform([med_plot], transform=transform)[0])

    ax.plot(0, med_plot, marker="o", markersize=5, mec="black", mfc="white", zorder=5)
    ax.text(
        0,
        med_plot,
        f"{int(round(max(0.0, med_raw))):,}",
        ha="center",
        va="bottom",
        fontsize=8,
        zorder=6,
    )

    yticks = ax.get_yticks()
    if transform == "log10p1":
        new_labels = [f"{int(max(0, round((10.0 ** y) - 1.0))):,}" for y in yticks]
    elif transform == "log1p":
        new_labels = [f"{int(max(0, round(np.expm1(y)))):,}" for y in yticks]
    else:
        raise ValueError("transform must be 'log10p1' or 'log1p'")

    ax.set_yticks(yticks)
    ax.set_yticklabels(new_labels)

    plt.xticks(rotation=45, ha="right")
    plt.xlabel("")
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


# =============================================================================
# Robust SnapATAC2 plotting wrappers
# =============================================================================

def safe_snap_tsse_plot(adata: ad.AnnData, out_file: str, force_plot: bool):
    if (not force_plot) and os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return
    try:
        if "n_fragment" not in adata.obs.columns or "tsse" not in adata.obs.columns:
            sys.stderr.write(f"WARNING: missing n_fragment/tsse for TSSE plot: {out_file}\n")
            return

        nfrag = adata.obs["n_fragment"].to_numpy()
        tsse = adata.obs["tsse"].to_numpy()

        ok = np.isfinite(nfrag) & np.isfinite(tsse) & (nfrag > 0)
        n_ok = int(ok.sum())
        if n_ok < 50:
            sys.stderr.write(f"WARNING: too few cells with n_fragment>0 for TSSE KDE ({n_ok}); skipping {out_file}\n")
            return

        adp = adata[ok].copy()
        snap.pl.tsse(
            adp,
            min_fragment=1,
            width=750,
            height=600,
            interactive=False,
            show=False,
            out_file=out_file,
        )
    except Exception as e:
        sys.stderr.write(f"WARNING: TSSE plot failed ({out_file}): {e}\n")

def safe_snap_frag_size_plot(adata: ad.AnnData, out_file: str, force_plot: bool):
    if (not force_plot) and os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return
    try:
        snap.pl.frag_size_distr(
            adata,
            width=750,
            height=600,
            interactive=False,
            show=False,
            out_file=out_file,
        )
    except Exception as e:
        sys.stderr.write(f"WARNING: frag-size plot failed ({out_file}): {e}\n")


# =============================================================================
# DNA: fragments caching + import/metrics + plots + filtering
# =============================================================================

def dna_ensure_fragments(workdir: str, bam_path: str, exp_name: str,
                         shift_left: int, shift_right: int,
                         force_fragments: bool):
    frag_tsv = os.path.join(workdir, f"{exp_name}.fragments.tsv.gz")
    if force_fragments or (not os.path.exists(frag_tsv)) or os.path.getsize(frag_tsv) == 0:
        info = snap.pp.make_fragment_file(
            bam_path,
            frag_tsv,
            is_paired=True,
            barcode_tag="CB",
            shift_left=shift_left,
            shift_right=shift_right,
            min_mapq=20,
            chunk_size=500000000,
            chrM=["chrM", "M"],
            tempdir=workdir,
        )
        with open(os.path.join(workdir, f"{exp_name}.frag_info.csv"), "w", newline="") as f:
            w = csv.writer(f)
            for k, v in info.items():
                w.writerow([k, v])
    return frag_tsv

def dna_prepare_raw(workdir: str, exp_name: str, genome: str, threads: int):
    raw_h5ad = os.path.join(workdir, f"{exp_name}.raw.h5ad")
    frag_tsv = os.path.join(workdir, f"{exp_name}.fragments.tsv.gz")

    if os.path.exists(raw_h5ad) and os.path.getsize(raw_h5ad) > 0:
        return sc.read_h5ad(raw_h5ad)

    genome_obj = snap.genome.hg38 if genome == "hg38" else snap.genome.mm39

    adata = snap.pp.import_fragments(
        frag_tsv,
        genome_obj,
        min_num_fragments=0,
        sorted_by_barcode=True,
        chrM=["chrM", "M"],
        shift_left=0,
        shift_right=0,
        chunk_size=50000,
        tempdir=workdir,
        backend="hdf5",
        n_jobs=threads,
    )

    snap.metrics.tsse(adata, gene_anno=genome_obj, inplace=True)
    snap.metrics.frag_size_distr(adata, add_key="frag_size_distr", inplace=True)

    adata.obs["log10_n_fragment"] = np.log10(adata.obs["n_fragment"].to_numpy() + 1.0)

    adata.write(raw_h5ad)
    return adata

def dna_plot_frag_violin(adata: ad.AnnData, exp_name: str, out_pdf: str, tag: str,
                         force_plot: bool, max_ticks: int):
    if (not force_plot) and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 0:
        return

    nfrag = adata.obs["n_fragment"].to_numpy()
    med = float(np.median(nfrag)) if nfrag.size else 0.0

    sc.pl.violin(adata, "log10_n_fragment", rotation=1e-7, show=False)
    fig = plt.gcf()
    for ax in fig.axes:
        set_log_transformed_axis_ticks(
            ax,
            transform="log10p1",
            raw_values=nfrag,
            max_ticks=max_ticks,
            drop_zero=True,
            max_quantile=0.99,
        )
        ax.set_title(
            f"{exp_name} ({tag})\n"
            f"n_cells: {adata.n_obs:,}\n"
            f"Median fragments/cell: {int(round(med)):,}"
        )
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close(fig)

def dna_filter_and_save(raw_adata: ad.AnnData, exp_name: str, workdir: str,
                        minFrags: int, maxFrags: int,
                        figdir: str, force: bool):
    pre_violin = os.path.join(figdir, f"{exp_name}_FragViolin_PreFilter.pdf")
    pre_tsse = os.path.join(figdir, f"{exp_name}_tsseDensity_PreFilter.pdf")
    pre_frag = os.path.join(figdir, f"{exp_name}_frag_size_distr_PreFilter.pdf")

    post_violin = os.path.join(figdir, f"{exp_name}_FragViolin_PostFilter.pdf")
    scatter_pdf = os.path.join(figdir, f"{exp_name}_tsse_scatter.pdf")

    dna_plot_frag_violin(raw_adata, exp_name, pre_violin, tag="PRE", force_plot=force, max_ticks=34)
    safe_snap_tsse_plot(raw_adata, pre_tsse, force_plot=force)
    safe_snap_frag_size_plot(raw_adata, pre_frag, force_plot=force)

    filt_h5ad = os.path.join(workdir, f"{exp_name}.h5ad")
    if (not force) and os.path.exists(filt_h5ad) and os.path.getsize(filt_h5ad) > 0:
        filt = sc.read_h5ad(filt_h5ad)
    else:
        mask = snap.pp.filter_cells(raw_adata, min_counts=minFrags, max_counts=maxFrags, min_tsse=0, inplace=False)
        filt = raw_adata[mask].copy()
        if "log10_n_fragment" not in filt.obs.columns:
            filt.obs["log10_n_fragment"] = np.log10(filt.obs["n_fragment"].to_numpy() + 1.0)
        filt.write(filt_h5ad)

    # Script-B-style POST violin
    bstyle_single_violin(
        filt,
        obs_key="log10_n_fragment",
        out_pdf=post_violin,
        title=exp_name,
        sample_label=exp_name,
        y_label="Fragments per cell",
        transform="log10p1",
        force_plot=force,
    )

    if force or (not os.path.exists(scatter_pdf)) or os.path.getsize(scatter_pdf) == 0:
        try:
            sc.pl.scatter(filt, "log10_n_fragment", "tsse", title=exp_name, show=False)
            plt.savefig(scatter_pdf)
            plt.close()
        except Exception as e:
            sys.stderr.write(f"WARNING: scatter plot failed for {exp_name}: {e}\n")

    return filt


# =============================================================================
# RNA: raw caching + QC/filtering + plots
# =============================================================================

def rna_prepare_raw(workdir: str, rna_sample: str):
    raw_h5ad = os.path.join(workdir, f"{rna_sample}.raw.h5ad")
    if os.path.exists(raw_h5ad) and os.path.getsize(raw_h5ad) > 0:
        return sc.read_h5ad(raw_h5ad)

    adata = read_star_solo_filtered(workdir, rna_sample)
    adata.write(raw_h5ad)
    return adata

def rna_add_qc_flags(adata: ad.AnnData, species: str):
    if species == "human":
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains(r"^HB[^ (P)]", regex=True)
    elif species == "mouse":
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
        adata.var["hb"] = adata.var_names.str.contains(r"^Hb[^ (P)]", regex=True)
    else:
        raise ValueError("species must be human|mouse")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

def rna_violin_log1p_as_raw(adata: ad.AnnData, key_log1p: str,
                            out_pdf: str, title_prefix: str,
                            force_plot: bool, max_ticks: int):
    if (not force_plot) and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 0:
        return

    raw_key = "total_counts" if key_log1p == "log1p_total_counts" else "n_genes_by_counts"
    raw_vals = adata.obs[raw_key].to_numpy()
    med = float(np.median(raw_vals)) if raw_vals.size else 0.0

    sc.pl.violin(adata, [key_log1p], jitter=0.4, rotation=1e-7, show=False)
    fig = plt.gcf()
    for ax in fig.axes:
        set_log_transformed_axis_ticks(
            ax,
            transform="log1p",
            raw_values=raw_vals,
            max_ticks=max_ticks,
            drop_zero=True,
            max_quantile=0.99,
        )
        ax.set_title(
            f"{title_prefix}\n"
            f"n_cells: {adata.n_obs:,}\n"
            f"Median {raw_key}: {int(round(med)):,}"
        )
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close(fig)

def rna_percent_violin(adata: ad.AnnData, key: str, out_pdf: str, title: str, force_plot: bool):
    if (not force_plot) and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 0:
        return
    try:
        vals = adata.obs[key].to_numpy()
        med = float(np.median(vals)) if vals.size else 0.0
        sc.pl.violin(adata, [key], jitter=0.4, rotation=1e-7, show=False)
        fig = plt.gcf()
        for ax in fig.axes:
            ax.set_title(f"{title}\nn_cells: {adata.n_obs:,}\nMedian {key}: {med:.2f}%")
        plt.tight_layout()
        plt.savefig(out_pdf)
        plt.close(fig)
    except Exception as e:
        sys.stderr.write(f"WARNING: RNA percent violin failed ({key}): {e}\n")

def rna_filter_and_save(raw_adata: ad.AnnData, rna_sample: str, workdir: str, figdir: str,
                        species: str,
                        minUMI: int, maxUMI: int,
                        max_pct_mt: float, max_pct_hb: float, max_pct_50g: float,
                        force: bool):
    """
    If qc_filtered exists and not force:
      - load it
      - STILL generate plots (PRE from freshly QC-tagged raw_adata copy; POST from loaded filtered)
    If missing or force:
      - compute QC, filter, write
      - generate plots
    """
    out_h5ad = os.path.join(workdir, f"{rna_sample}.qc_filtered.h5ad")

    # Always (re)compute QC metrics on a raw copy for PRE plots + filtering logic
    adata = raw_adata.copy()
    rna_add_qc_flags(adata, species=species)

    # Plot paths
    pre_mt = os.path.join(figdir, f"{rna_sample}_RNA_pct_mt_PreFilter.pdf")
    pre_top50 = os.path.join(figdir, f"{rna_sample}_RNA_pct_top50_PreFilter.pdf")
    pre_umi = os.path.join(figdir, f"{rna_sample}_RNA_UMI_PreFilter.pdf")
    pre_genes = os.path.join(figdir, f"{rna_sample}_RNA_Genes_PreFilter.pdf")
    pre_scatter = os.path.join(figdir, f"{rna_sample}_RNA_Scatter_PreFilter.pdf")

    post_mt = os.path.join(figdir, f"{rna_sample}_RNA_pct_mt_PostFilter.pdf")
    post_top50 = os.path.join(figdir, f"{rna_sample}_RNA_pct_top50_PostFilter.pdf")
    post_umi = os.path.join(figdir, f"{rna_sample}_RNA_UMI_PostFilter.pdf")
    post_genes = os.path.join(figdir, f"{rna_sample}_RNA_Genes_PostFilter.pdf")
    post_scatter = os.path.join(figdir, f"{rna_sample}_RNA_Scatter_PostFilter.pdf")

    # Load or compute filtered
    if (not force) and os.path.exists(out_h5ad) and os.path.getsize(out_h5ad) > 0:
        ad_f = sc.read_h5ad(out_h5ad)

        # Safety: ensure required QC fields exist in cached filtered object
        need = ["pct_counts_mt", "pct_counts_hb", "pct_counts_in_top_50_genes",
                "log1p_total_counts", "log1p_n_genes_by_counts", "total_counts", "n_genes_by_counts"]
        missing = [k for k in need if k not in ad_f.obs.columns]
        if missing:
            # Recompute QC on cached filtered (cheap relative to rerunning everything)
            rna_add_qc_flags(ad_f, species=species)

    else:
        keep = (adata.obs["total_counts"] >= minUMI) & (adata.obs["total_counts"] <= maxUMI)
        keep &= (adata.obs["pct_counts_mt"] < max_pct_mt)
        keep &= (adata.obs["pct_counts_in_top_50_genes"] < max_pct_50g)
        keep &= (adata.obs["pct_counts_hb"] < max_pct_hb)

        ad_f = adata[keep].copy()
        ad_f.write(out_h5ad)

    # --- PRE plots (unchanged) ---
    rna_percent_violin(adata, "pct_counts_mt", pre_mt, f"{rna_sample} pct_counts_mt (PRE)", force_plot=force)
    rna_percent_violin(adata, "pct_counts_in_top_50_genes", pre_top50, f"{rna_sample} pct_top50 (PRE)", force_plot=force)

    rna_violin_log1p_as_raw(
        adata, "log1p_total_counts", pre_umi,
        f"{rna_sample} UMI per cell (PRE)",
        force_plot=force, max_ticks=44,
    )
    rna_violin_log1p_as_raw(
        adata, "log1p_n_genes_by_counts", pre_genes,
        f"{rna_sample} Genes per cell (PRE)",
        force_plot=force, max_ticks=36,
    )

    if force or (not os.path.exists(pre_scatter)) or os.path.getsize(pre_scatter) == 0:
        try:
            sc.pl.scatter(
                adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                title=rna_sample + " (PRE)\nColor: pct_counts_mt", show=False
            )
            plt.savefig(pre_scatter)
            plt.close()
        except Exception as e:
            sys.stderr.write(f"WARNING: RNA scatter PRE failed ({rna_sample}): {e}\n")

    # --- POST plots ---
    rna_percent_violin(ad_f, "pct_counts_mt", post_mt, f"{rna_sample} pct_counts_mt (POST)", force_plot=force)
    rna_percent_violin(ad_f, "pct_counts_in_top_50_genes", post_top50, f"{rna_sample} pct_top50 (POST)", force_plot=force)

    # Script-B-style POST count violins
    bstyle_single_violin(
        ad_f,
        obs_key="log1p_total_counts",
        out_pdf=post_umi,
        title="RNA",
        sample_label=rna_sample,
        y_label="UMI counts per cell",
        transform="log1p",
        force_plot=force,
    )
    bstyle_single_violin(
        ad_f,
        obs_key="log1p_n_genes_by_counts",
        out_pdf=post_genes,
        title="RNA",
        sample_label=rna_sample,
        y_label="Genes per cell",
        transform="log1p",
        force_plot=force,
    )

    if force or (not os.path.exists(post_scatter)) or os.path.getsize(post_scatter) == 0:
        try:
            sc.pl.scatter(
                ad_f, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                title=rna_sample + " (POST)\nColor: pct_counts_mt", show=False
            )
            plt.savefig(post_scatter)
            plt.close()
        except Exception as e:
            sys.stderr.write(f"WARNING: RNA scatter POST failed ({rna_sample}): {e}\n")

    return ad_f



# =============================================================================
# Overlap plots: UpSet + Venn
# =============================================================================

def passing_cells_upset(sets: dict, title: str, out_pdf: str, min_subset_size="0.2%"):
    all_items = set()
    for s in sets.values():
        all_items |= s
    if not all_items:
        sys.stderr.write(f"WARNING: No items for UpSet: {title}\n")
        return

    data = {k: [item in v for item in all_items] for k, v in sets.items()}
    df = pd.DataFrame(data, index=list(all_items))

    fig = plt.figure(figsize=(9, 10))
    plot(
        from_indicators(df.columns, data=df),
        sort_by="cardinality",
        show_percentages=True,
        min_subset_size=min_subset_size,
        facecolor="darkblue",
        fig=fig,
        element_size=None,
    )
    plt.title(title, fontsize=14, fontweight="bold")
    plt.savefig(out_pdf)
    plt.close(fig)

def make_venn_plots(sets: dict, title_prefix: str, out_dir: str,
                    dna_min_frags: int = None, dna_max_frags: int = None,
                    rna_min_umi: int = None, rna_max_umi: int = None,
                    rna_max_pct_mt: float = None, rna_max_pct_top50: float = None):
    keys = list(sets.keys())
    if len(keys) < 3:
        return

    # Script-B-style main 3-set Venn when RNA + 2 DNA modalities are present
    if len(keys) == 3 and "RNA" in sets:
        non_rna = [k for k in keys if k != "RNA"]
        if len(non_rna) == 2:
            k2, k3 = non_rna[0], non_rna[1]
            a, b, c = sets["RNA"], sets[k2], sets[k3]

            plt.figure(figsize=(8, 8))
            venn3(
                [a, b, c],
                set_labels=(
                    f"RNA\n{rna_min_umi}<UMIs<{rna_max_umi}\n<{rna_max_pct_mt}% MT\n<{rna_max_pct_top50}% 50TopGenes",
                    f"{k2}\n{dna_min_frags}<fragments<{dna_max_frags}",
                    f"{k3}\n{dna_min_frags}<fragments<{dna_max_frags}",
                )
            )
            plt.suptitle(
                f"Venn Diagram of Passing Cells\nfor RNA, {k2}, and {k3}",
                fontsize=16,
                fontweight="bold",
            )
            plt.tight_layout()
            out_pdf = os.path.join(out_dir, f"{title_prefix}_VennPassingCells.pdf")
            plt.savefig(out_pdf)
            plt.close()
            return

    # Fallback to original behavior for other cases
    def _save_venn3(k1, k2, k3):
        a, b, c = sets[k1], sets[k2], sets[k3]
        plt.figure(figsize=(8, 8))
        venn3(
            [a, b, c],
            set_labels=(f"{k1}\n{len(a)}", f"{k2}\n{len(b)}", f"{k3}\n{len(c)}")
        )
        plt.suptitle(f"{title_prefix}\n{k1} vs {k2} vs {k3}", fontsize=14, fontweight="bold")
        plt.tight_layout()
        out_pdf = os.path.join(out_dir, f"{title_prefix}_Venn_{k1}_{k2}_{k3}.pdf".replace("/", "_"))
        plt.savefig(out_pdf)
        plt.close()

    if len(keys) == 3:
        _save_venn3(keys[0], keys[1], keys[2])
        return

    if "RNA" in sets:
        non = [k for k in keys if k != "RNA"]
        for k2, k3 in itertools.combinations(non, 2):
            _save_venn3("RNA", k2, k3)
    else:
        combos = list(itertools.combinations(keys, 3))[:10]
        for k1, k2, k3 in combos:
            _save_venn3(k1, k2, k3)


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description="Single-cell processing across RNA + DNA marks per sample.")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--mo-map", required=True)
    ap.add_argument("--genome", required=True, choices=["hg38", "mm39"])
    ap.add_argument("--species", required=True, choices=["human", "mouse"])

    ap.add_argument("--dna-samples", default="")
    ap.add_argument("--rna-samples", default="")
    ap.add_argument("--pairs-tsv", default="")
    ap.add_argument("--threads", type=int, default=64)

    ap.add_argument("--dna-min-frags", type=int, default=400)
    ap.add_argument("--dna-max-frags", type=int, default=10000)
    ap.add_argument("--dna-shift-left", type=int, default=4)
    ap.add_argument("--dna-shift-right", type=int, default=-5)

    ap.add_argument("--rna-min-umi", type=int, default=500)
    ap.add_argument("--rna-max-umi", type=int, default=8000)
    ap.add_argument("--rna-max-pct-mt", type=float, default=5)
    ap.add_argument("--rna-max-pct-hb", type=float, default=5.0)
    ap.add_argument("--rna-max-pct-top50", type=float, default=25.0)

    ap.add_argument("--figdir", default="", help="Default: {workdir}/figures_sc")

    ap.add_argument("--force", action="store_true")
    ap.add_argument("--force-fragments", action="store_true")

    ap.add_argument("--harmonize", action="store_true",
                    help="Rewrite DT DNA inline SB (3nt) to RNA SB (4nt) in obs_names, keeping injected base.")
    ap.add_argument("--harmonize-min-frac", type=float, default=0.95)
    ap.add_argument("--harmonize-strict", action="store_true")

    args = ap.parse_args()

    workdir = args.workdir
    ensure_dir(workdir)

    figdir = args.figdir if args.figdir else os.path.join(workdir, "figures_sc")
    ensure_dir(figdir)

    sample_to_marks = read_mo_map(args.mo_map)

    def parse_list(x):
        if not x:
            return []
        if "," in x:
            return [s.strip() for s in x.split(",") if s.strip()]
        return [s.strip() for s in x.split() if s.strip()]

    dna_samples = parse_list(args.dna_samples) or sorted(sample_to_marks.keys())
    rna_samples = parse_list(args.rna_samples) or infer_rna_samples(workdir)

    pairs_map = {}
    if args.pairs_tsv:
        with open(args.pairs_tsv, "r") as f:
            for line in f:
                line = line.strip()
                if (not line) or line.startswith("#"):
                    continue
                a = line.split("\t")
                if len(a) >= 2:
                    pairs_map[a[0]] = a[1]

    for dna_sample in dna_samples:
        marks = sample_to_marks.get(dna_sample, [])
        if not marks:
            sys.stderr.write(f"WARNING: No marks in mo_map for {dna_sample}, skipping\n")
            continue

        # Prefer explicit pairs, else tonsil pairing, else legacy fallback
        rna_sample = pairs_map.get(dna_sample, "")
        if not rna_sample:
            rna_sample = infer_pairing_tonsil(dna_sample, rna_samples)
        if not rna_sample:
            rna_sample = infer_pairing_legacy(dna_sample, rna_samples)

        if not rna_sample:
            sys.stderr.write(f"WARNING: Could not pair RNA sample for {dna_sample}. Skipping RNA integration.\n")

        group_fig = os.path.join(figdir, dna_sample)
        ensure_dir(group_fig)

        print(f"\n=== GROUP {dna_sample} ===")
        print(f"DNA marks: {marks}")
        print(f"RNA paired: {rna_sample if rna_sample else 'NONE'}")

        # Keep your old ST-vs-DT heuristic (won't trigger for tonsil)
        is_st_group = ("STD" in dna_sample)

        # -------------------------
        # DNA
        # -------------------------
        dna_filtered = {}
        for mark in marks:
            bam = detect_dna_bam(workdir, dna_sample, mark)
            if not bam:
                sys.stderr.write(f"WARNING: Missing BAM for {dna_sample} {mark}, skipping mark\n")
                continue

            bn = os.path.basename(bam)
            exp_name = bn.replace(".bam", "").replace("_NoDup", "").replace("_MarkedDup", "")

            dna_ensure_fragments(
                workdir=workdir,
                bam_path=bam,
                exp_name=exp_name,
                shift_left=args.dna_shift_left,
                shift_right=args.dna_shift_right,
                force_fragments=args.force_fragments,
            )

            raw = dna_prepare_raw(
                workdir=workdir,
                exp_name=exp_name,
                genome=args.genome,
                threads=args.threads,
            )

            filt = dna_filter_and_save(
                raw_adata=raw,
                exp_name=exp_name,
                workdir=workdir,
                minFrags=args.dna_min_frags,
                maxFrags=args.dna_max_frags,
                figdir=group_fig,
                force=args.force,
            )

            dna_filtered[mark] = (filt, exp_name)

        if not dna_filtered:
            sys.stderr.write(f"WARNING: No DNA marks processed for {dna_sample}, skipping group\n")
            continue

        # -------------------------
        # RNA
        # -------------------------
        rna_filt = None
        if rna_sample:
            rna_raw = rna_prepare_raw(workdir, rna_sample)
            rna_filt = rna_filter_and_save(
                raw_adata=rna_raw,
                rna_sample=rna_sample,
                workdir=workdir,
                figdir=group_fig,
                species=args.species,
                minUMI=args.rna_min_umi,
                maxUMI=args.rna_max_umi,
                max_pct_mt=args.rna_max_pct_mt,
                max_pct_hb=args.rna_max_pct_hb,
                max_pct_50g=args.rna_max_pct_top50,
                force=args.force,
            )

        # -------------------------
        # Harmonize DT DNA -> RNA SB (tonsil DNA is DT)
        # -------------------------
        if args.harmonize and (not is_st_group):
            print(f"[HARMONIZE] Rewriting DNA DT3 -> RNA4 at offset={INLINE_OFFSET} (keeping injected base).")

            for mark, (ad_obj, exp_name) in dna_filtered.items():
                rewrite_inline_sb_to_rna_sb(
                    ad_obj,
                    assay="dna_dt",
                    inline_offset=INLINE_OFFSET,
                    min_recognized_frac=args.harmonize_min_frac,
                    strict=args.harmonize_strict,
                )
                out_h5ad = os.path.join(workdir, f"{exp_name}.h5ad")
                ad_obj.write(out_h5ad)

            if rna_filt is not None:
                rewrite_inline_sb_to_rna_sb(
                    rna_filt,
                    assay="rna",
                    inline_offset=INLINE_OFFSET,
                    min_recognized_frac=args.harmonize_min_frac,
                    strict=args.harmonize_strict,
                )

                # Backward compatibility: PDTD-only injected-base swap
                if "PDTD" in dna_sample:
                    print("[PDTD FIX] Swapping injected base C<->G in RNA ONLY for PDTD group.")
                    swap_injected_base_CG_inplace(rna_filt)

                out_h5ad = os.path.join(workdir, f"{rna_sample}.qc_filtered.h5ad")
                rna_filt.write(out_h5ad)

        elif args.harmonize and is_st_group:
            print("[HARMONIZE] ST-like group detected ('STD' in name). Leaving barcodes unchanged.")

        # -------------------------
        # Passing sets + overlap plots
        # -------------------------
        sets = {}
        for mark, (ad_obj, _) in dna_filtered.items():
            sets[mark] = set(ad_obj.obs_names.astype(str))
        if rna_filt is not None:
            sets["RNA"] = set(rna_filt.obs_names.astype(str))

        upset_pdf = os.path.join(group_fig, f"{dna_sample}_UpsetPassingCells.pdf")
        passing_cells_upset(
            sets=sets,
            title=f"Passing Cells for {dna_sample}\n" + " + ".join(list(sets.keys())),
            out_pdf=upset_pdf,
            min_subset_size="0.2%",
        )

        make_venn_plots(
            sets,
            title_prefix=dna_sample,
            out_dir=group_fig,
            dna_min_frags=args.dna_min_frags,
            dna_max_frags=args.dna_max_frags,
            rna_min_umi=args.rna_min_umi,
            rna_max_umi=args.rna_max_umi,
            rna_max_pct_mt=args.rna_max_pct_mt,
            rna_max_pct_top50=args.rna_max_pct_top50,
        )

        shared = None
        for s in sets.values():
            shared = s if shared is None else (shared & s)
        shared = shared if shared is not None else set()

        print("Passing cells per set: " + ", ".join([f"{k}={len(v)}" for k, v in sets.items()]))
        print(f"Shared across all ({len(sets)} sets): {len(shared)}")

        shared_path = os.path.join(group_fig, f"{dna_sample}_SharedCells.lst")
        with open(shared_path, "w") as f:
            for bc in sorted(shared):
                f.write(bc + "\n")

        # Save shared h5ad objects and make separate Script-B-style shared-cell violins
        for mark, (ad_obj, exp_name) in dna_filtered.items():
            ad_s = ad_obj[ad_obj.obs_names.isin(shared)].copy()
            out_h5ad = os.path.join(workdir, f"{exp_name}.shared.h5ad")
            ad_s.write(out_h5ad)

            if ad_s.n_obs > 0:
                shared_pdf = os.path.join(group_fig, f"{exp_name}_FragViolin_PostFilter_Shared.pdf")
                bstyle_single_violin(
                    ad_s,
                    obs_key="log10_n_fragment",
                    out_pdf=shared_pdf,
                    title=exp_name,
                    sample_label=exp_name + "_shared",
                    y_label="Fragments per cell",
                    transform="log10p1",
                    force_plot=args.force,
                )

        if rna_filt is not None:
            rna_s = rna_filt[rna_filt.obs_names.isin(shared)].copy()
            out_h5ad = os.path.join(workdir, f"{rna_sample}.shared.h5ad")
            rna_s.write(out_h5ad)

            if rna_s.n_obs > 0:
                shared_umi = os.path.join(group_fig, f"{rna_sample}_RNA_UMI_PostFilter_Shared.pdf")
                shared_genes = os.path.join(group_fig, f"{rna_sample}_RNA_Genes_PostFilter_Shared.pdf")

                bstyle_single_violin(
                    rna_s,
                    obs_key="log1p_total_counts",
                    out_pdf=shared_umi,
                    title="RNA",
                    sample_label=rna_sample + "_shared",
                    y_label="UMI counts per cell",
                    transform="log1p",
                    force_plot=args.force,
                )
                bstyle_single_violin(
                    rna_s,
                    obs_key="log1p_n_genes_by_counts",
                    out_pdf=shared_genes,
                    title="RNA",
                    sample_label=rna_sample + "_shared",
                    y_label="Genes per cell",
                    transform="log1p",
                    force_plot=args.force,
                )

    print("\nDONE")


if __name__ == "__main__":
    main()




python3 SingleCell_Process_TP.py \
  --workdir /mnt/dataFast/ahrmad/triseq_202601/processed/ \
  --mo-map /mnt/dataFast/ahrmad/triseq_202601/WLs/mo_map_tonsil.tsv \
  --genome hg38 \
  --species human \
  --harmonize \
  --force