import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc

sns.set(style="whitegrid")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def _first_present(A: ad.AnnData, candidates):
    for c in candidates:
        if c in A.obs.columns:
            return c
    return None

def fmt_log_ticks(ax):
    ymin, ymax = ax.get_ylim()
    raw_min = max(0.0, 10**ymin - 1.0)
    raw_max = max(0.0, 10**ymax - 1.0)

    candidates = np.array([
         500,
        1000, 2500, 5000,
        10000
    ], dtype=float)

    ticks_raw = candidates[(candidates >= raw_min) & (candidates <= raw_max)]

    if ticks_raw.size < 2:
        rmax = max(raw_max, 1.0)
        ticks_raw = np.array([0, rmax / 2, rmax])

    ticks_raw = np.unique(ticks_raw)
    ticks_log = np.log10(ticks_raw + 1.0)

    ax.set_yticks(ticks_log)
    ax.set_yticklabels([f"{int(tr):,}" for tr in ticks_raw])

def _annotate_global_count_and_median(ax, values, y_is_log=False):
    ymin, ymax = ax.get_ylim()
    headroom = 0.12 * (ymax - ymin)
    ax.set_ylim(ymin, ymax + headroom)
    ypos = ymax + 0.04 * (ymax - ymin)

    v = np.asarray(values)
    if v.size == 0:
        return

    if y_is_log:
        med_raw = int(np.median(10**v - 1))
        txt = f"{v.size} cells\n{med_raw:,} median"
    else:
        med = float(np.median(v))
        txt = f"{v.size} cells\n{med:.2f} median"

    ax.text(0, ypos, txt, ha="center", va="bottom", fontsize=10)

def violin_global(A: ad.AnnData, value_col: str, ylabel: str, outpath: str,
                  title: str, log10_plus1: bool):
    if value_col not in A.obs.columns:
        print(f"[skip] {outpath} because {value_col} not found")
        return

    vals = pd.to_numeric(A.obs[value_col], errors="coerce").fillna(0)
    y = np.log10(vals.to_numpy() + 1) if log10_plus1 else vals.to_numpy()

    df = pd.DataFrame({"group": ["all"] * len(y), "value": y})

    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    sns.violinplot(data=df, x="group", y="value", inner=None, ax=ax)
    sns.boxplot(
        data=df, x="group", y="value",
        width=0.25, showcaps=True, dodge=False,
        boxprops={'facecolor': 'none', 'linewidth': 1.5},
        whiskerprops={'linewidth': 1.5},
        medianprops={'color': 'black'},
        flierprops={'marker': '', 'markersize': 0},
        ax=ax
    )

    if log10_plus1:
        fmt_log_ticks(ax)

    _annotate_global_count_and_median(ax, y, y_is_log=log10_plus1)

    ax.set_xlabel("")
    ax.set_xticklabels(["All cells"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

def make_qc_plots_one_modality(rna_ad: ad.AnnData, dna_ad: ad.AnnData,
                               dna_label: str, outdir: str, prefix: str):
    os.makedirs(outdir, exist_ok=True)

    # RNA
    umi_col = _first_present(rna_ad, ["total_counts", "n_counts", "nCount_RNA", "umi", "counts"])
    genes_col = _first_present(rna_ad, ["n_genes_by_counts", "nFeature_RNA", "genes"])

    if umi_col is not None:
        violin_global(
            rna_ad, umi_col, "UMIs per cell",
            os.path.join(outdir, f"{prefix}_RNA_UMI_global.pdf"),
            title=f"RNA – UMIs per Cell\nCells: {rna_ad.n_obs:,}",
            log10_plus1=True
        )

    if genes_col is not None:
        violin_global(
            rna_ad, genes_col, "Genes per cell",
            os.path.join(outdir, f"{prefix}_RNA_nGenes_global.pdf"),
            title=f"RNA – Genes per Cell\nCells: {rna_ad.n_obs:,}",
            log10_plus1=True
        )

    # DNA
    frag_col = _first_present(dna_ad, ["n_fragment", "passed_filters", "nCount_ATAC"])
    tsse_col = _first_present(dna_ad, ["tsse", "tss_enrichment", "TSSE"])

    if frag_col is not None:
        violin_global(
            dna_ad, frag_col, "Fragments per cell",
            os.path.join(outdir, f"{prefix}_{dna_label}_Fragments_global.pdf"),
            title=f"{dna_label} – Fragments per Cell\nCells: {dna_ad.n_obs:,}",
            log10_plus1=True
        )

    if tsse_col is not None:
        violin_global(
            dna_ad, tsse_col, "TSSE",
            os.path.join(outdir, f"{prefix}_{dna_label}_TSSE_global.pdf"),
            title=f"{dna_label} – TSSE\nCells: {dna_ad.n_obs:,}",
            log10_plus1=False
        )

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
RNA_FULL = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"

CONFIG = {
    "H3K27ac": {
        "full_rna": RNA_FULL,
        "full_dna": "DNA_merged/MergedTonsils_H3K27ac_merged_processed_barcodeRewritten_withRNAlabels.h5ad",
        "pair_rna": "overlap_objects/MergedTonsils_RNA_overlap_H3K27ac_barcodeRewritten_Reclust.h5ad",
        "pair_dna": "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "overlap3_rna": "overlap_objects/MergedTonsils_RNA_overlap_AC_ME3_reUMAP.h5ad",
        "overlap3_dna": "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_ME3.h5ad",
    },
    "H3K27me3": {
        "full_rna": RNA_FULL,
        "full_dna": "DNA_merged/MergedTonsils_H3K27me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad",
        "pair_rna": "overlap_objects/MergedTonsils_RNA_overlap_H3K27me3_barcodeRewritten_Reclust.h5ad",
        "pair_dna": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "overlap3_rna": "overlap_objects/MergedTonsils_RNA_overlap_AC_ME3_reUMAP.h5ad",
        "overlap3_dna": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_AC.h5ad",
    },
    "H3K9me3": {
        "full_rna": RNA_FULL,
        "full_dna": "DNA_merged/MergedTonsils_H3K9me3_merged_processed_barcodeRewritten_withRNAlabels.h5ad",
        "pair_rna": "overlap_objects/MergedTonsils_RNA_overlap_H3K9me3_barcodeRewritten_Reclust.h5ad",
        "pair_dna": "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
        "overlap3_rna": None,
        "overlap3_dna": None,
    },
}

BASE_OUTDIR = "QC_by_modality"

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------
for modality, cfg in CONFIG.items():
    print("\n" + "=" * 80)
    print(f"Processing {modality}")
    print("=" * 80)

    # ---------- full ----------
    rna_full = sc.read_h5ad(cfg["full_rna"])
    dna_full = sc.read_h5ad(cfg["full_dna"])
    make_qc_plots_one_modality(
        rna_ad=rna_full,
        dna_ad=dna_full,
        dna_label=modality,
        outdir=os.path.join(BASE_OUTDIR, modality, "full"),
        prefix=f"{modality}_Full",
    )

    # ---------- pairwise RNA ∩ modality ----------
    rna_pair = sc.read_h5ad(cfg["pair_rna"])
    dna_pair = sc.read_h5ad(cfg["pair_dna"])
    make_qc_plots_one_modality(
        rna_ad=rna_pair,
        dna_ad=dna_pair,
        dna_label=modality,
        outdir=os.path.join(BASE_OUTDIR, modality, "pairwise_overlap"),
        prefix=f"{modality}_PairwiseOverlap",
    )

    # ---------- overlap3 = RNA ∩ H3K27ac ∩ H3K27me3 ----------
    if cfg["overlap3_rna"] is not None and cfg["overlap3_dna"] is not None:
        if os.path.exists(cfg["overlap3_rna"]) and os.path.exists(cfg["overlap3_dna"]):
            rna_ov3 = sc.read_h5ad(cfg["overlap3_rna"])
            dna_ov3 = sc.read_h5ad(cfg["overlap3_dna"])
            make_qc_plots_one_modality(
                rna_ad=rna_ov3,
                dna_ad=dna_ov3,
                dna_label=modality,
                outdir=os.path.join(BASE_OUTDIR, modality, "overlap3"),
                prefix=f"{modality}_Overlap3",
            )
        else:
            print(f"[skip] overlap3 files missing for {modality}")
    else:
        print(f"[skip] no overlap3 defined for {modality}")