import os
import pandas as pd
import scanpy as sc
import snapatac2 as snap

BL = "/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz"
out_dir = "bw_tracks"
os.makedirs(out_dir, exist_ok=True)

DNA_OVERLAP_FILES = {
    "H3K27ac":  "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten_Reclust.h5ad",
    "H3K27me3": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
    "H3K9me3":  "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
}

group_col = "leiden_merged_type"

for mod, path in DNA_OVERLAP_FILES.items():
    print("\n" + "=" * 80)
    print(f"Processing {mod}")
    print("=" * 80)

    adata = sc.read_h5ad(path)

    # per-cell-type tracks, excluding NotInRNA if present
    groups = pd.Series(adata.obs[group_col].astype(str)).unique().tolist()
    selections = sorted([g for g in groups if g != "NotInRNA"])

    print(f"Exporting bigWigs for {len(selections)} groups:", selections)

    snap.ex.export_coverage(
        adata,
        groupby=group_col,
        selections=selections,
        bin_size=20,
        blacklist=BL,
        normalization="RPKM",
        include_for_norm=None,
        exclude_for_norm=None,
        min_frag_length=None,
        max_frag_length=2000,
        counting_strategy="fragment",
        smooth_base=None,
        out_dir=out_dir,
        prefix=f"{mod}_",
        suffix=".bw",
        output_format="bigwig",
        compression=None,
        compression_level=None,
        tempdir=None,
        n_jobs=8,
    )

    # optional: one full track per modality across all overlap cells
    adata.obs["all_overlap_cells"] = "all_overlap_cells"
    snap.ex.export_coverage(
        adata,
        groupby="all_overlap_cells",
        selections=["all_overlap_cells"],
        bin_size=20,
        blacklist=BL,
        normalization="RPKM",
        include_for_norm=None,
        exclude_for_norm=None,
        min_frag_length=None,
        max_frag_length=2000,
        counting_strategy="fragment",
        smooth_base=None,
        out_dir=out_dir,
        prefix=f"{mod}_FULL_",
        suffix=".bw",
        output_format="bigwig",
        compression=None,
        compression_level=None,
        tempdir=None,
        n_jobs=8,
    )





import os
import pandas as pd
import scanpy as sc
import snapatac2 as snap

BL = "/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz"
out_dir = "bw_tracks_overlap3"
os.makedirs(out_dir, exist_ok=True)

DNA_OVERLAP3_FILES = {
    "H3K27ac":  "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_ME3_reUMAP_snap.h5ad",
    "H3K27me3": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_AC_reUMAP_snap.h5ad",
}

group_col = "leiden_merged_type"

for mod, path in DNA_OVERLAP3_FILES.items():
    print("\n" + "=" * 80)
    print(f"Processing {mod}")
    print("=" * 80)

    adata = sc.read_h5ad(path)

    if group_col not in adata.obs.columns:
        raise KeyError(f"{group_col} not found in {path}")

    # per-cell-type tracks, excluding NotInRNA if present
    groups = pd.Series(adata.obs[group_col].astype(str)).unique().tolist()
    selections = sorted([g for g in groups if g != "NotInRNA"])

    print(f"Exporting bigWigs for {len(selections)} groups:", selections)

    snap.ex.export_coverage(
        adata,
        groupby=group_col,
        selections=selections,
        bin_size=20,
        blacklist=BL,
        normalization="RPKM",
        include_for_norm=None,
        exclude_for_norm=None,
        min_frag_length=None,
        max_frag_length=2000,
        counting_strategy="fragment",
        smooth_base=None,
        out_dir=out_dir,
        prefix=f"Overlap3_{mod}_",
        suffix=".bw",
        output_format="bigwig",
        compression=None,
        compression_level=None,
        tempdir=None,
        n_jobs=8,
    )

    # one full track per modality across all overlap3 cells
    adata.obs["all_overlap3_cells"] = "all_overlap3_cells"
    snap.ex.export_coverage(
        adata,
        groupby="all_overlap3_cells",
        selections=["all_overlap3_cells"],
        bin_size=20,
        blacklist=BL,
        normalization="RPKM",
        include_for_norm=None,
        exclude_for_norm=None,
        min_frag_length=None,
        max_frag_length=2000,
        counting_strategy="fragment",
        smooth_base=None,
        out_dir=out_dir,
        prefix=f"Overlap3_{mod}_FULL_",
        suffix=".bw",
        output_format="bigwig",
        compression=None,
        compression_level=None,
        tempdir=None,
        n_jobs=8,
    )


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from collections import defaultdict

import pandas as pd
import scanpy as sc

# --------------------------------------------------
# INPUTS
# --------------------------------------------------
OVERLAP_RNA_H5AD = "overlap_objects/MergedTonsils_RNA_overlap_AC_ME3_reUMAP.h5ad"
GROUP_COL = "leiden_merged_type"

BAM_MAP = {
    "Tonsil1": "BAMs/Tonsil1.bam",
    "Tonsil2": "BAMs/Sc_VTR_HumanTonsil2.filtered_cells.bam",
    "Tonsil3_K27ac": "BAMs/Sc_TP3r_TP3_A.filtered_cells.bam",
    "Tonsil3_K9me3": "BAMs/Sc_TP3r_TP3_B.filtered_cells.bam",
    "Tonsil4_K27ac": "BAMs/Sc_TP4r_TP4_A.filtered_cells.bam",
    "Tonsil4_K9me3": "BAMs/Sc_TP4r_TP4_B.filtered_cells.bam",
}

OUTDIR = "RNA_overlap3_bam_bw"
TMPDIR = os.path.join(OUTDIR, "tmp")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)

BIN_SIZE = 20
N_THREADS = 72

# --------------------------------------------------
# HELPERS
# --------------------------------------------------
def run(cmd):
    print("\n[run]", " ".join(cmd))
    subprocess.run(cmd, check=True)

def sample_and_barcode(cellname: str):
    sample, bc = cellname.rsplit("_", 1)
    return sample, bc

def write_lines(path, values):
    with open(path, "w") as f:
        for v in values:
            f.write(f"{v}\n")

def subset_bam_by_cb(in_bam, barcode_txt, out_bam, threads=8):
    # Requires samtools with tag-file support
    run([
        "samtools", "view",
        "-@", str(threads),
        "-b",
        "-D", f"CB:{barcode_txt}",
        "-o", out_bam,
        in_bam,
    ])
    run(["samtools", "index", "-@", str(threads), out_bam])

def merge_bams(bam_list, out_bam, threads=8):
    if len(bam_list) == 0:
        raise ValueError(f"No BAMs to merge for {out_bam}")
    if len(bam_list) == 1:
        run(["cp", bam_list[0], out_bam])
        run(["cp", bam_list[0] + ".bai", out_bam + ".bai"])
    else:
        run(["samtools", "merge", "-@", str(threads), "-f", out_bam] + bam_list)
        run(["samtools", "index", "-@", str(threads), out_bam])

def bam_to_bw_rpkm(in_bam, out_bw, bin_size=20, threads=8):
    run([
        "bamCoverage",
        "-b", in_bam,
        "-o", out_bw,
        "--binSize", str(bin_size),
        "--normalizeUsing", "RPKM",
        "-p", str(threads),
    ])

# --------------------------------------------------
# LOAD OVERLAP RNA CELLS
# --------------------------------------------------
adata = sc.read_h5ad(OVERLAP_RNA_H5AD)

if GROUP_COL not in adata.obs.columns:
    raise KeyError(f"{GROUP_COL} not found in {OVERLAP_RNA_H5AD}")

cells = pd.Index(adata.obs_names.astype(str))
meta = pd.DataFrame(index=cells)
meta["cell"] = cells
meta["sample"], meta["barcode"] = zip(*meta["cell"].map(sample_and_barcode))
meta[GROUP_COL] = adata.obs[GROUP_COL].astype(str).values

missing_samples = sorted(set(meta["sample"]) - set(BAM_MAP))
if missing_samples:
    raise ValueError(f"These sample prefixes are not in BAM_MAP: {missing_samples}")

meta.to_csv(os.path.join(OUTDIR, "overlap3_cells_metadata.csv"), index=False)

print("\nCells per sample:")
print(meta["sample"].value_counts())

print("\nCells per cell type:")
print(meta[GROUP_COL].value_counts())

# --------------------------------------------------
# 1) SUBSET EACH SAMPLE BAM FOR ALL OVERLAP CELLS
# --------------------------------------------------
sample_subset_bams = []

for sample, subdf in meta.groupby("sample", sort=True):
    barcodes = sorted(subdf["barcode"].unique())
    barcode_txt = os.path.join(TMPDIR, f"{sample}.overlap3.barcodes.txt")
    out_bam = os.path.join(TMPDIR, f"{sample}.overlap3.bam")

    write_lines(barcode_txt, barcodes)
    subset_bam_by_cb(BAM_MAP[sample], barcode_txt, out_bam, threads=N_THREADS)
    sample_subset_bams.append(out_bam)

# --------------------------------------------------
# 2) MERGE ALL OVERLAP RNA BAMs + MAKE POOLED BW
# --------------------------------------------------
merged_all_bam = os.path.join(OUTDIR, "RNA_overlap3_all.bam")
merged_all_bw = os.path.join(OUTDIR, "RNA_overlap3_all.RPKM.bw")

merge_bams(sample_subset_bams, merged_all_bam, threads=N_THREADS)
bam_to_bw_rpkm(merged_all_bam, merged_all_bw, bin_size=BIN_SIZE, threads=N_THREADS)

# --------------------------------------------------
# 3) PER CELL TYPE: SUBSET EACH SAMPLE BAM, MERGE, MAKE BW
# --------------------------------------------------
group_summary = []

for ct, ct_df in meta.groupby(GROUP_COL, sort=True):
    safe_ct = "".join(c if c.isalnum() or c in "._-+" else "_" for c in ct)

    ct_sample_bams = []

    for sample, subdf in ct_df.groupby("sample", sort=True):
        barcodes = sorted(subdf["barcode"].unique())
        barcode_txt = os.path.join(TMPDIR, f"{safe_ct}.{sample}.barcodes.txt")
        out_bam = os.path.join(TMPDIR, f"{safe_ct}.{sample}.bam")

        write_lines(barcode_txt, barcodes)
        subset_bam_by_cb(BAM_MAP[sample], barcode_txt, out_bam, threads=N_THREADS)
        ct_sample_bams.append(out_bam)

    ct_bam = os.path.join(OUTDIR, f"RNA_overlap3_{safe_ct}.bam")
    ct_bw = os.path.join(OUTDIR, f"RNA_overlap3_{safe_ct}.RPKM.bw")

    merge_bams(ct_sample_bams, ct_bam, threads=N_THREADS)
    bam_to_bw_rpkm(ct_bam, ct_bw, bin_size=BIN_SIZE, threads=N_THREADS)

    group_summary.append({
        "cell_type": ct,
        "n_cells": int(ct_df.shape[0]),
        "n_samples_present": int(ct_df["sample"].nunique()),
        "bam": ct_bam,
        "bw": ct_bw,
    })

# --------------------------------------------------
# 4) SUMMARY
# --------------------------------------------------
summary_df = pd.DataFrame(group_summary).sort_values("n_cells", ascending=False)
summary_df.to_csv(os.path.join(OUTDIR, "RNA_overlap3_per_celltype_summary.csv"), index=False)

print("\nDone.")
print(f"Pooled BAM: {merged_all_bam}")
print(f"Pooled BW:  {merged_all_bw}")
print(summary_df)