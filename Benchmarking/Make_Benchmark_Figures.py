#!/usr/bin/env python3
"""
Requires NumPy, pandas, SciPy, Matplotlib, AnnData, Scanpy, SnapATAC2 2.10.0, pysam, MACS3, samtools, GNU sort, awk.
SnapATAC2 downloads the TSS annotations on first use.
"""
from pathlib import Path

# Editable paths and source metadata.
# Data roots
DATA_ROOT = Path('/mnt/dataFast/ahrmad/ScH2R_202409/FullBench')
NAS_ROOT = Path('/mnt/dataNas/ahrmad')

# Reference and output directories
REFERENCE_DIR = Path(__file__).resolve().parent / "references"
OUTPUT_ROOT = Path(__file__).resolve().parent / 'benchmark_figures_regenerated'
SCRATCH_ROOT = Path('/tmp')

# Input paths
INPUTS = {
    'InHouse_NanoCNT_H3K27ac': DATA_ROOT / 'InHouse_NanoCNT/ScKDMA_H3K27ac.tsv.gz',
    'InHouse_NanoCNT_H3K27me3': DATA_ROOT / 'InHouse_NanoCNT/ScKDMA_H3K27me3.tsv.gz',
    'MTRseq_H3K27ac': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27ac_NoDup.bam',
    'MTRseq_H3K27ac_index': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27ac_NoDup.bam.bai',
    'MTRseq_H3K27me3': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27me3_NoDup.bam',
    'MTRseq_H3K27me3_index': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27me3_NoDup.bam.bai',
    'MTRseq_RNA_matrix.mtx': NAS_ROOT / 'MTRseq/rna_align/MTRseq1_A.Solo.outGeneFull/filtered/matrix.mtx',
    'MTRseq_RNA_barcodes.tsv': NAS_ROOT / 'MTRseq/rna_align/MTRseq1_A.Solo.outGeneFull/filtered/barcodes.tsv',
    'MTRseq_RNA_features.tsv': NAS_ROOT / 'MTRseq/rna_align/MTRseq1_A.Solo.outGeneFull/filtered/features.tsv',
    'MTRseq_OverNight_H3K27ac': NAS_ROOT / 'MTRseq_OverNight/dna_align/MTRseq_OvN_A_H3K27ac_NoDup.bam',
    'MTRseq_OverNight_H3K27ac_index': NAS_ROOT / 'MTRseq_OverNight/dna_align/MTRseq_OvN_A_H3K27ac_NoDup.bam.bai',
    'MTRseq_OverNight_H3K27me3': NAS_ROOT / 'MTRseq_OverNight/dna_align/MTRseq_OvN_A_H3K27me3_NoDup.bam',
    'MTRseq_OverNight_H3K27me3_index': NAS_ROOT / 'MTRseq_OverNight/dna_align/MTRseq_OvN_A_H3K27me3_NoDup.bam.bai',
    'MTRseq_OverNight_RNA_matrix.mtx': NAS_ROOT / 'MTRseq_OverNight/rna_align/MTRseq_OvN_A.Solo.outGeneFull/filtered/matrix.mtx.gz',
    'MTRseq_OverNight_RNA_barcodes.tsv': NAS_ROOT / 'MTRseq_OverNight/rna_align/MTRseq_OvN_A.Solo.outGeneFull/filtered/barcodes.tsv.gz',
    'MTRseq_OverNight_RNA_features.tsv': NAS_ROOT / 'MTRseq_OverNight/rna_align/MTRseq_OvN_A.Solo.outGeneFull/filtered/features.tsv.gz',
    'TrES_H3K27ac': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27ac_NoDup.bam',
    'TrES_H3K27ac_index': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27ac_NoDup.bam.bai',
    'TrES_H3K27me3': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27me3_NoDup.bam',
    'TrES_H3K27me3_index': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27me3_NoDup.bam.bai',
    'TrES_RNA_matrix.mtx': NAS_ROOT / 'TrES_K422/rna_align/Ori_K422.Solo.outGeneFull/filtered/matrix.mtx.gz',
    'TrES_RNA_barcodes.tsv': NAS_ROOT / 'TrES_K422/rna_align/Ori_K422.Solo.outGeneFull/filtered/barcodes.tsv.gz',
    'TrES_RNA_features.tsv': NAS_ROOT / 'TrES_K422/rna_align/Ori_K422.Solo.outGeneFull/filtered/features.tsv.gz',
    'uCoTargetX_H3K27ac': NAS_ROOT / 'uCOTx/dna_align/uCOTx1_A_H3K27ac_NoDup.bam',
    'uCoTargetX_H3K27ac_index': NAS_ROOT / 'uCOTx/dna_align/uCOTx1_A_H3K27ac_NoDup.bam.bai',
    'uCoTargetX_H3K27me3': NAS_ROOT / 'uCOTx/dna_align/uCOTx1_A_H3K27me3_NoDup.bam',
    'uCoTargetX_H3K27me3_index': NAS_ROOT / 'uCOTx/dna_align/uCOTx1_A_H3K27me3_NoDup.bam.bai',
    'uCoTargetX_RNA_matrix.mtx': NAS_ROOT / 'uCOTx/rna_align/uCOTx1_A.Solo.outGeneFull/filtered/matrix.mtx.gz',
    'uCoTargetX_RNA_barcodes.tsv': NAS_ROOT / 'uCOTx/rna_align/uCOTx1_A.Solo.outGeneFull/filtered/barcodes.tsv.gz',
    'uCoTargetX_RNA_features.tsv': NAS_ROOT / 'uCOTx/rna_align/uCOTx1_A.Solo.outGeneFull/filtered/features.tsv.gz',
    'uCoTargetX_OverNight_H3K27ac': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27ac_NoDup.bam',
    'uCoTargetX_OverNight_H3K27ac_index': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27ac_NoDup.bam.bai',
    'uCoTargetX_OverNight_H3K27me3': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27me3_NoDup.bam',
    'uCoTargetX_OverNight_H3K27me3_index': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27me3_NoDup.bam.bai',
    'uCoTargetX_OverNight_RNA_matrix.mtx': NAS_ROOT / 'uCOTx_OverNight/rna_align/uCOTx_OvN_A.Solo.outGeneFull/filtered/matrix.mtx.gz',
    'uCoTargetX_OverNight_RNA_barcodes.tsv': NAS_ROOT / 'uCOTx_OverNight/rna_align/uCOTx_OvN_A.Solo.outGeneFull/filtered/barcodes.tsv.gz',
    'uCoTargetX_OverNight_RNA_features.tsv': NAS_ROOT / 'uCOTx_OverNight/rna_align/uCOTx_OvN_A.Solo.outGeneFull/filtered/features.tsv.gz',
    'public_NTT_ac': DATA_ROOT / 'NTT/GSM6538033_hek_k562_H3K27ac.bed.gz',
    'public_NTT_me3': DATA_ROOT / 'NTT/GSM6538032_hek_k562_H3K27me3.bed.gz',
    'public_NanoCNT_ac': DATA_ROOT / 'NanoCNT/GSE198467_H3K27ac_fragments.tsv.gz',
    'public_NanoCNT_me3': DATA_ROOT / 'NanoCNT/GSE198467_H3K27me3_fragments.tsv.gz',
    'public_10xMultiome_RNA': DATA_ROOT / '10xMultiome/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
    'RT_S1_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529552_Endo.S1/matrix.mtx.gz',
    'RT_S1_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529552_Endo.S1/barcodes.tsv.gz',
    'RT_S1_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529552_Endo.S1/features.tsv.gz',
    'RT_S2_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529553_Endo.S2/matrix.mtx.gz',
    'RT_S2_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529553_Endo.S2/barcodes.tsv.gz',
    'RT_S2_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529553_Endo.S2/features.tsv.gz',
    'RT_S3_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529554_Endo.S3/matrix.mtx.gz',
    'RT_S3_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529554_Endo.S3/barcodes.tsv.gz',
    'RT_S3_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529554_Endo.S3/features.tsv.gz',
    'RT_S4_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529555_Endo.S4/matrix.mtx.gz',
    'RT_S4_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529555_Endo.S4/barcodes.tsv.gz',
    'RT_S4_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529555_Endo.S4/features.tsv.gz',
    'RT_S5_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529556_Endo.S5/matrix.mtx.gz',
    'RT_S5_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529556_Endo.S5/barcodes.tsv.gz',
    'RT_S5_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529556_Endo.S5/features.tsv.gz',
    'RT_S6_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529557_Endo.S6/matrix.mtx.gz',
    'RT_S6_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529557_Endo.S6/barcodes.tsv.gz',
    'RT_S6_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529557_Endo.S6/features.tsv.gz',
    'RT_S7_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529558_Endo.S7/matrix.mtx.gz',
    'RT_S7_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529558_Endo.S7/barcodes.tsv.gz',
    'RT_S7_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529558_Endo.S7/features.tsv.gz',
    'RT_S8_matrix': DATA_ROOT / 'MTR/Technical/rna/GSM8529559_Endo.S8/matrix.mtx.gz',
    'RT_S8_barcodes': DATA_ROOT / 'MTR/Technical/rna/GSM8529559_Endo.S8/barcodes.tsv.gz',
    'RT_S8_features': DATA_ROOT / 'MTR/Technical/rna/GSM8529559_Endo.S8/features.tsv.gz',
    'K422_K9_DT_A_H3K27me3': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9A_A_H3K27me3_NoDup.bam',
    'K422_K9_DT_A_H3K27me3_index': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9A_A_H3K27me3_NoDup.bam.bai',
    'K422_K9_DT_A_H3K9me3': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9A_A_H3K9me3_NoDup.bam',
    'K422_K9_DT_A_H3K9me3_index': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9A_A_H3K9me3_NoDup.bam.bai',
    'K422_K9_DT_A_RNA_matrix.mtx': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9A_A.Solo.outGeneFull/filtered/matrix.mtx',
    'K422_K9_DT_A_RNA_barcodes.tsv': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9A_A.Solo.outGeneFull/filtered/barcodes.tsv',
    'K422_K9_DT_A_RNA_features.tsv': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9A_A.Solo.outGeneFull/filtered/features.tsv',
    'K422_K9_DT_B_H3K27me3': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9B_B_H3K27me3_NoDup.bam',
    'K422_K9_DT_B_H3K27me3_index': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9B_B_H3K27me3_NoDup.bam.bai',
    'K422_K9_DT_B_H3K9me3': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9B_B_H3K9me3_NoDup.bam',
    'K422_K9_DT_B_H3K9me3_index': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9B_B_H3K9me3_NoDup.bam.bai',
    'K422_K9_DT_B_RNA_matrix.mtx': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9B_B.Solo.outGeneFull/filtered/matrix.mtx',
    'K422_K9_DT_B_RNA_barcodes.tsv': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9B_B.Solo.outGeneFull/filtered/barcodes.tsv',
    'K422_K9_DT_B_RNA_features.tsv': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9B_B.Solo.outGeneFull/filtered/features.tsv',
    'K422_K9_DT_C_H3K27me3': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9C_C_H3K27me3_NoDup.bam',
    'K422_K9_DT_C_H3K27me3_index': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9C_C_H3K27me3_NoDup.bam.bai',
    'K422_K9_DT_C_H3K9me3': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9C_C_H3K9me3_NoDup.bam',
    'K422_K9_DT_C_H3K9me3_index': NAS_ROOT / 'TrES_K422_K9_DT/dna_align/K9C_C_H3K9me3_NoDup.bam.bai',
    'K422_K9_DT_C_RNA_matrix.mtx': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9C_C.Solo.outGeneFull/filtered/matrix.mtx',
    'K422_K9_DT_C_RNA_barcodes.tsv': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9C_C.Solo.outGeneFull/filtered/barcodes.tsv',
    'K422_K9_DT_C_RNA_features.tsv': NAS_ROOT / 'TrES_K422_K9_DT/rna_align/K9C_C.Solo.outGeneFull/filtered/features.tsv',
    'hg38_blacklist': REFERENCE_DIR / 'hg38-blacklist.v2.bed.gz',
    'MTR_S1_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529544_HIST_S1_H3K27ac.fragments.tsv.gz',
    'MTR_S2_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529545_HIST_S2_H3K27ac.fragments.tsv.gz',
    'MTR_S3_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529546_HIST_S3_H3K27ac.fragments.tsv.gz',
    'MTR_S4_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529547_HIST_S4_H3K27ac.fragments.tsv.gz',
    'MTR_S5_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529548_HIST_S5_H3K27ac.fragments.tsv.gz',
    'MTR_S6_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529549_HIST_S6_H3K27ac.fragments.tsv.gz',
    'MTR_S7_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529550_HIST_S7_H3K27ac.fragments.tsv.gz',
    'MTR_S8_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529551_HIST_S8_H3K27ac.fragments.tsv.gz',
    'MTR_S1_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529544_HIST_S1_H3K27me3.fragments.tsv.gz',
    'MTR_S2_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529545_HIST_S2_H3K27me3.fragments.tsv.gz',
    'MTR_S3_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529546_HIST_S3_H3K27me3.fragments.tsv.gz',
    'MTR_S4_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529547_HIST_S4_H3K27me3.fragments.tsv.gz',
    'MTR_S5_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529548_HIST_S5_H3K27me3.fragments.tsv.gz',
    'MTR_S6_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529549_HIST_S6_H3K27me3.fragments.tsv.gz',
    'MTR_S7_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529550_HIST_S7_H3K27me3.fragments.tsv.gz',
    'MTR_S8_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529551_HIST_S8_H3K27me3.fragments.tsv.gz',
}

# File URLs or accessions
SOURCES = {
    'InHouse_NanoCNT_H3K27ac': 'GSE324511',
    'InHouse_NanoCNT_H3K27me3': 'GSE324511',
    'MTRseq_H3K27ac': 'GSE324511',
    'MTRseq_H3K27ac_index': 'GSE324511',
    'MTRseq_H3K27me3': 'GSE324511',
    'MTRseq_H3K27me3_index': 'GSE324511',
    'MTRseq_RNA_matrix.mtx': 'GSE324511',
    'MTRseq_RNA_barcodes.tsv': 'GSE324511',
    'MTRseq_RNA_features.tsv': 'GSE324511',
    'MTRseq_OverNight_H3K27ac': 'GSE324511',
    'MTRseq_OverNight_H3K27ac_index': 'GSE324511',
    'MTRseq_OverNight_H3K27me3': 'GSE324511',
    'MTRseq_OverNight_H3K27me3_index': 'GSE324511',
    'MTRseq_OverNight_RNA_matrix.mtx': 'GSE324511',
    'MTRseq_OverNight_RNA_barcodes.tsv': 'GSE324511',
    'MTRseq_OverNight_RNA_features.tsv': 'GSE324511',
    'TrES_H3K27ac': 'GSE324511',
    'TrES_H3K27ac_index': 'GSE324511',
    'TrES_H3K27me3': 'GSE324511',
    'TrES_H3K27me3_index': 'GSE324511',
    'TrES_RNA_matrix.mtx': 'GSE324511',
    'TrES_RNA_barcodes.tsv': 'GSE324511',
    'TrES_RNA_features.tsv': 'GSE324511',
    'uCoTargetX_H3K27ac': 'GSE324511',
    'uCoTargetX_H3K27ac_index': 'GSE324511',
    'uCoTargetX_H3K27me3': 'GSE324511',
    'uCoTargetX_H3K27me3_index': 'GSE324511',
    'uCoTargetX_RNA_matrix.mtx': 'GSE324511',
    'uCoTargetX_RNA_barcodes.tsv': 'GSE324511',
    'uCoTargetX_RNA_features.tsv': 'GSE324511',
    'uCoTargetX_OverNight_H3K27ac': 'GSE324511',
    'uCoTargetX_OverNight_H3K27ac_index': 'GSE324511',
    'uCoTargetX_OverNight_H3K27me3': 'GSE324511',
    'uCoTargetX_OverNight_H3K27me3_index': 'GSE324511',
    'uCoTargetX_OverNight_RNA_matrix.mtx': 'GSE324511',
    'uCoTargetX_OverNight_RNA_barcodes.tsv': 'GSE324511',
    'uCoTargetX_OverNight_RNA_features.tsv': 'GSE324511',
    'public_NTT_ac': 'GSE212588',
    'public_NTT_me3': 'GSE212588',
    'public_NanoCNT_ac': 'GSE198467',
    'public_NanoCNT_me3': 'GSE198467',
    'public_10xMultiome_RNA': 'https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
    'RT_S1_matrix': 'GSM8529552',
    'RT_S1_barcodes': 'GSM8529552',
    'RT_S1_features': 'GSM8529552',
    'RT_S2_matrix': 'GSM8529553',
    'RT_S2_barcodes': 'GSM8529553',
    'RT_S2_features': 'GSM8529553',
    'RT_S3_matrix': 'GSM8529554',
    'RT_S3_barcodes': 'GSM8529554',
    'RT_S3_features': 'GSM8529554',
    'RT_S4_matrix': 'GSM8529555',
    'RT_S4_barcodes': 'GSM8529555',
    'RT_S4_features': 'GSM8529555',
    'RT_S5_matrix': 'GSM8529556',
    'RT_S5_barcodes': 'GSM8529556',
    'RT_S5_features': 'GSM8529556',
    'RT_S6_matrix': 'GSM8529557',
    'RT_S6_barcodes': 'GSM8529557',
    'RT_S6_features': 'GSM8529557',
    'RT_S7_matrix': 'GSM8529558',
    'RT_S7_barcodes': 'GSM8529558',
    'RT_S7_features': 'GSM8529558',
    'RT_S8_matrix': 'GSM8529559',
    'RT_S8_barcodes': 'GSM8529559',
    'RT_S8_features': 'GSM8529559',
    'K422_K9_DT_A_H3K27me3': 'GSE324511',
    'K422_K9_DT_A_H3K27me3_index': 'GSE324511',
    'K422_K9_DT_A_H3K9me3': 'GSE324511',
    'K422_K9_DT_A_H3K9me3_index': 'GSE324511',
    'K422_K9_DT_A_RNA_matrix.mtx': 'GSE324511',
    'K422_K9_DT_A_RNA_barcodes.tsv': 'GSE324511',
    'K422_K9_DT_A_RNA_features.tsv': 'GSE324511',
    'K422_K9_DT_B_H3K27me3': 'GSE324511',
    'K422_K9_DT_B_H3K27me3_index': 'GSE324511',
    'K422_K9_DT_B_H3K9me3': 'GSE324511',
    'K422_K9_DT_B_H3K9me3_index': 'GSE324511',
    'K422_K9_DT_B_RNA_matrix.mtx': 'GSE324511',
    'K422_K9_DT_B_RNA_barcodes.tsv': 'GSE324511',
    'K422_K9_DT_B_RNA_features.tsv': 'GSE324511',
    'K422_K9_DT_C_H3K27me3': 'GSE324511',
    'K422_K9_DT_C_H3K27me3_index': 'GSE324511',
    'K422_K9_DT_C_H3K9me3': 'GSE324511',
    'K422_K9_DT_C_H3K9me3_index': 'GSE324511',
    'K422_K9_DT_C_RNA_matrix.mtx': 'GSE324511',
    'K422_K9_DT_C_RNA_barcodes.tsv': 'GSE324511',
    'K422_K9_DT_C_RNA_features.tsv': 'GSE324511',
    'hg38_blacklist': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
    'MTR_S1_H3K27ac': 'GSM8529544',
    'MTR_S2_H3K27ac': 'GSM8529545',
    'MTR_S3_H3K27ac': 'GSM8529546',
    'MTR_S4_H3K27ac': 'GSM8529547',
    'MTR_S5_H3K27ac': 'GSM8529548',
    'MTR_S6_H3K27ac': 'GSM8529549',
    'MTR_S7_H3K27ac': 'GSM8529550',
    'MTR_S8_H3K27ac': 'GSM8529551',
    'MTR_S1_H3K27me3': 'GSM8529544',
    'MTR_S2_H3K27me3': 'GSM8529545',
    'MTR_S3_H3K27me3': 'GSM8529546',
    'MTR_S4_H3K27me3': 'GSM8529547',
    'MTR_S5_H3K27me3': 'GSM8529548',
    'MTR_S6_H3K27me3': 'GSM8529549',
    'MTR_S7_H3K27me3': 'GSM8529550',
    'MTR_S8_H3K27me3': 'GSM8529551',
}
# End editable block.

THREADS = 16



import argparse
import bisect
from collections import Counter
from contextlib import ExitStack
from itertools import groupby, product
import gzip
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import gaussian_kde
from scipy.optimize import brentq, least_squares
from scipy.ndimage import distance_transform_edt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from GEO_Preparation.merge_libraries import fragment_path, library_pattern, rna_paths, split_cell_id

LOG = logging.getLogger("benchmark_figures")
MARKS = ("H3K27ac", "H3K27me3")
MODS = (*MARKS, "RNA")
K422_K9_DT_MARKS = ("H3K27me3", "H3K9me3")
K422_K9_DATASET = "TrES_K422_K9_DT"
K422_K9_LIBRARIES = {"K9A_A": "A", "K9B_B": "B", "K9C_C": "C"}
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
NUCLEAR = AUTOSOMES | {"chrX", "chrY"}
MM = 1 / 25.4
INK = "#26323C"
MUTED = "#64717D"
QC = {"H3K27ac_min": 200, "H3K27me3_min": 300, "H3K9me3_min": 300,
      "dna_max": 15000, "rna_min_umi": 1000, "rna_max_umi": 25000,
      "max_pct_mt": 10, "max_pct_hb": 5, "max_pct_top50": 30}
CONFIG = {"qc": QC, "fragment_length_max": 1000, "libraries": {"A": {}, "B": {}, "C": {}}}
GROUP_IH = "benchmark_TrES_uCoTargetX_MTRseq_full"
GROUP_PUBLIC = "figures_publication_mtr_fixed_no_onecell"
GROUP_FRIP = "figures_frip_publication_mtr_fixed_no_onecell"
GROUP_K422_K9_DT = "TrES_K422_DualTag_K27_K9_showcase"
TECHS = ("TrES", "uCoTargetX", "uCoTargetX_OverNight", "MTRseq", "MTRseq_OverNight", "InHouse_NanoCNT")
INHOUSE = {tech: {"dna": {mark: f"{tech}_{mark}" for mark in MARKS},
                  **({"rna": [f"{tech}_RNA_{stem}" for stem in ("matrix.mtx", "barcodes.tsv", "features.tsv")]}
                     if tech != "InHouse_NanoCNT" else {})} for tech in TECHS}
PUBLIC = {tech: {"dna": {mark: (tuple(f"MTR_S{i}_{mark}" for i in range(1, 9)) if tech == "MTR" else
                        f"TrES_{mark}" if tech == "TrES" else f"public_{tech}_{suffix}")
                       for mark, suffix in zip(MARKS, ("ac", "me3"))}}
          for tech in ("TrES", "MTR", "NanoCNT", "NTT")}
PUBLIC["TrES"]["rna"] = INHOUSE["TrES"]["rna"]
PUBLIC["10xMultiome"] = {"rna": "public_10xMultiome_RNA"}
MTR_RT = {f"S{i}": {name: f"RT_S{i}_{name}" for name in ("matrix", "barcodes", "features")} for i in range(1, 9)}
MTR_RNA_INPUTS = {key for library in MTR_RT.values() for key in library.values()}
MTR_MIN_FEATURES, MTR_MIN_CELLS_PER_GENE = 100, 3
K422_K9_DT = {lib: {"dna": {mark: f"K422_K9_DT_{lib}_{mark}" for mark in K422_K9_DT_MARKS},
                  "rna": [f"K422_K9_DT_{lib}_RNA_{stem}" for stem in ("matrix.mtx", "barcodes.tsv", "features.tsv")]}
            for lib in ("A", "B", "C")}
EXPERIMENT_INFO = {
    "TrES": ("TrES", "not specified", "TrES"),
    "uCoTargetX": ("uCoTargetX", "1 h", "uCoTargetX (1 h binding)"),
    "uCoTargetX_OverNight": ("uCoTargetX", "overnight", "uCoTargetX (overnight binding)"),
    "MTRseq": ("MTRseq", "1 h", "MTRseq (1 h binding)"),
    "MTRseq_OverNight": ("MTRseq", "overnight", "MTRseq (overnight binding)"),
    "InHouse_NanoCNT": ("NanoCNT", "not specified", "In-house NanoCNT")}
IH_COLORS = {"TrES":"#315A79","uCoTargetX":"#238477","uCoTargetX_OverNight":"#67AD9F",
             "MTRseq":"#B5677A","MTRseq_OverNight":"#D298A7","InHouse_NanoCNT":"#B39051"}
MOD_COLORS = {"RNA":"#6679A9","H3K27ac":"#338F80","H3K27me3":"#B29660"}
PUB_COLORS = {"TrES":"#315A79","MTR":"#B5677A","NanoCNT":"#8B83A6","NTT":"#B47D59","10xMultiome":"#81A5AE"}
PUB_ORDER = list(PUB_COLORS)
LIB_COLORS = {"A":"#315A79","B":"#238477","C":"#B5677A","pooled":"#315A79"}
TAG_CANDIDATES = ("CB", "CR", "BC", "RG")
LIBRARIES = {f"GSM{8529551+i}":f"S{i}" for i in range(1,9)}
TIMEPOINTS = {"01":"D0","02":"D1","03":"D2","04":"D3"}
RNA_ID = re.compile(r"^(?:(GSM\d+|S[1-8])[_:])?(\d{2}:\d{2}:\d{2}):(\d{2})(?:-1)?$")

EXPECTED_PDFS = {
    'benchmark_TrES_uCoTargetX_MTRseq_full/QC_Venn_InHouse_NanoCNT.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/QC_Venn_MTRseq.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/QC_Venn_MTRseq_OverNight.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/QC_Venn_TrES.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/QC_Venn_uCoTargetX.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/QC_Venn_uCoTargetX_OverNight.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/RNA_detected_genes.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/RNA_UMIs.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/H3K27ac_unique_fragments.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/H3K27ac_fragment_lengths.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/H3K27me3_unique_fragments.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/H3K27me3_fragment_lengths.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/supplementary/H3K27ac_own_peak_FRiP.pdf',
    'benchmark_TrES_uCoTargetX_MTRseq_full/supplementary/H3K27me3_own_peak_FRiP.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_Combined_H3K27ac_Fragments_Violin_full.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_Combined_H3K27ac_TSSE_Violin_full.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_Combined_H3K27me3_Fragments_Violin_full.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_MTR_QC_Venn.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_NanoCNT_QC_Venn.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_NTT_QC_Venn.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_TrES_QC_Venn.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_RNA_n_genes_Violin_full.pdf',
    'figures_publication_mtr_fixed_no_onecell/Benchmark_RNA_UMIs_Violin_full.pdf',
    'figures_frip_publication_mtr_fixed_no_onecell/FRiPBench_BothMarks_full_own_FRiP_Violin.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/RNA_UMIs.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/RNA_detected_genes.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/H3K27me3_fragment_lengths.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/H3K27me3_own_peak_FRiP.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/H3K27me3_unique_fragments.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/H3K9me3_fragment_lengths.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/H3K9me3_own_peak_FRiP.pdf',
    'TrES_K422_DualTag_K27_K9_showcase/pooled/H3K9me3_unique_fragments.pdf',

}


def experiment_label(tech, multiline=False):
    label = EXPERIMENT_INFO.get(tech, (tech, "not specified", tech))[2]
    if multiline:
        return "In-house\nNanoCNT" if tech == "InHouse_NanoCNT" else label.replace(" (", "\n").removesuffix(")")
    return label


def fragments(path):
    """Yield one unique fragment row. PCR multiplicity never weights metrics."""
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{number}: fragments require chrom/start/end/barcode.")
            chrom, start, end, barcode = fields[:4]
            start, end = int(start), int(end)
            if start < 0 or end <= start or not barcode:
                raise ValueError(f"Invalid fragment at {path}:{number}")
                                                                                  
            chrom = "chr" + chrom if not chrom.startswith("chr") else chrom
            yield chrom, start, end, barcode


def merge_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(end, out[-1][1]))
        else:
            out.append((start, end))
    return out


class Regions:
    def __init__(self, bed=None, records=None):
        bychrom = {}
        if bed:
            with open_text(bed) as handle:
                for line in handle:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    f = line.split(); chrom = f[0]
                    if not chrom.startswith("chr"): chrom = "chr" + chrom
                    a, b = int(f[1]), int(f[2])
                    if a < 0 or b <= a: raise ValueError(f"Invalid BED interval: {line.strip()}")
                    bychrom.setdefault(chrom, []).append((a, b))
        if records:
            for chrom, a, b in records: bychrom.setdefault(chrom, []).append((a, b))
        self.data = {c: merge_intervals(v) for c, v in bychrom.items()}
        self.starts = {c: [a for a, _ in v] for c, v in self.data.items()}

    def overlaps(self, chrom, start, end):
        values = self.data.get(chrom)
        if not values: return False
        i = bisect.bisect_left(self.starts[chrom], end) - 1
        return i >= 0 and values[i][1] > start

    @property
    def bases(self):
        return sum(b-a for v in self.data.values() for a, b in v)

    def write(self, path):
        atomic_text(path, "".join(f"{c}\t{a}\t{b}\n" for c in sorted(self.data) for a,b in self.data[c]))


def detect_bam_tag(bam, rna_ids, explicit=None):
    import pysam
    if explicit: return explicit, {"selection": "explicit", "tag": explicit}
    values = {tag: [] for tag in ("CB", "RG", "BC")}; primary = 0
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
            primary += 1
            for tag in values:
                if read.has_tag(tag): values[tag].append(str(read.get_tag(tag)))
            if primary >= 10000: break
    rna = {re.sub(r"-1$", "", str(x)) for x in rna_ids}
    choices = []
    audit = {"sampled_primary_alignments": primary, "candidates": {}}
    for tag, v in values.items():
        unique_values = set(v); overlap = len({re.sub(r"-1$", "", x) for x in unique_values} & rna)
        audit["candidates"][tag] = {"tagged_alignments": len(v), "distinct_values": len(unique_values), "RNA_overlap": overlap}
        if len(v) >= .95*max(primary, 1) and overlap:
            choices.append((overlap/max(len(unique_values),1), len(unique_values), -("CB", "RG", "BC").index(tag), tag))
    if not choices:
        raise ValueError(f"Cannot confidently identify a cellular BAM tag in {bam}. Set barcode_tag explicitly in config. Audit: {audit}")
    tag = max(choices)[-1]; audit["tag"] = tag; return tag, audit


def align_barcodes(rna, dna_frames):
    """Only exact IDs or collision-free removal of a terminal 10x '-1'."""
    all_ids = set(rna.barcode.astype(str)) | set().union(*(set(f.source_barcode.astype(str)) for f in dna_frames.values()))
    canonical = {x: re.sub(r"-1$", "", x) for x in all_ids}
                                                                             
    for ids in [rna.barcode, *(f.source_barcode for f in dna_frames.values())]:
        changed = [canonical[str(x)] for x in ids]
        if len(set(changed)) != len(changed):
            raise ValueError("Barcode suffix normalization creates a collision. Supply consistent cell IDs before benchmarking.")
    rna = rna.copy(); rna["source_barcode"] = rna.barcode; rna["barcode"] = rna.barcode.map(canonical)
    result = {}
    for mod, frame in dna_frames.items():
        result[mod] = frame.copy(); result[mod]["barcode"] = frame.source_barcode.map(canonical)
    return rna, result


def apply_qc(rna, dna, qc):
    result={}
    if rna is not None:
        rna = rna.copy()
        rna["pass_qc"] = (rna.UMIs.between(qc["rna_min_umi"],qc["rna_max_umi"]) &
            rna.pct_mt.le(qc["max_pct_mt"]) & rna.pct_hb.le(qc["max_pct_hb"]) & rna.pct_top50.le(qc["max_pct_top50"]))
        result["RNA"]=rna
    for mark, frame in dna.items():
        frame = frame.copy(); frame["pass_qc"] = frame.unique_fragments.between(qc[mark+"_min"],qc["dna_max"])
        result[mark] = frame
    return result


def ih_style():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.titlesize":9,
        "axes.labelsize":8,"xtick.labelsize":7.3,"ytick.labelsize":7.3,"legend.fontsize":7,
        "axes.labelcolor":INK,"text.color":INK,"xtick.color":INK,"ytick.color":INK,
        "axes.edgecolor":"#8B969F","axes.linewidth":.65,"axes.spines.top":False,"axes.spines.right":False,
        "xtick.major.size":2.5,"ytick.major.size":2.5,"xtick.major.width":.65,"ytick.major.width":.65,
        "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","savefig.transparent":False,
        "figure.facecolor":"white","axes.facecolor":"white","path.simplify":False})


def ih_axis_finish(ax):
    ax.yaxis.grid(True,color="#EDF0F2",lw=.55);ax.set_axisbelow(True)
    ax.spines["bottom"].set_color("#B0B8BE")


def ih_limits(data,metric):
    values=data[metric].dropna().to_numpy()
    if not len(values):return (0,1)
    if metric.startswith("FRiP"):return (0,min(100,max(10,float(values.max()*100)*1.1)))
    if metric in ("UMIs","genes","unique_fragments"):
        return (np.log10(max(1,float(values.min())*.78)),np.log10(max(2,float(values.max())*1.28)))
    return (max(0,float(values.min())*.80),max(1,float(values.max())*1.15))


def ih_log_ticks(ax,ylim):
    ticks=[x for x in [100,200,300,500,1000,2000,3000,5000,10000,15000,20000,30000,50000] if ylim[0]<=np.log10(x)<=ylim[1]]
                                                                            
    if len(ticks)>5:ticks=ticks[::2]
    ax.set_yticks(np.log10(ticks),[f"{v/1000:g}k" if v>=1000 else str(v) for v in ticks])


def ih_violin(ax,frame,techs,metric,ylim,ylabel):
    log=metric in ("UMIs","genes","unique_fragments");groups=[]
    for tech in techs:
        v=frame.loc[frame.technology==tech,metric].dropna().to_numpy(dtype=float)
        if metric.startswith("FRiP"):v=v*100
        if log:v=np.log10(v[v>0])
        groups.append(v)
    pooled=np.concatenate(groups) if groups else np.array([])
    bandwidth=max(np.std(pooled)*max(len(pooled),2)**(-.2),.025 if log else .5)
    densities=[];grids=[]
    for v in groups:
        if len(v)>=3 and np.std(v,ddof=1)>1e-10:
            grid=np.linspace(v.min(),v.max(),150)
            kde=gaussian_kde(v,bw_method=bandwidth/np.std(v,ddof=1));d=kde(grid)
        else:grid=np.array([]);d=np.array([])
        grids.append(grid);densities.append(d)
    maximum=max([float(d.max()) for d in densities if len(d)]+[1e-12])
    labels=[]
    for x,(tech,v,grid,density) in enumerate(zip(techs,groups,grids,densities)):
        color=IH_COLORS[tech]
        if len(grid):
            width=.36*density/maximum
            ax.fill_betweenx(grid,x-width,x+width,color=color,alpha=.72,lw=.75,edgecolor=color,zorder=2)
        if len(v):
            rng=np.random.default_rng(731+x);display=rng.choice(len(v),min(len(v),180),replace=False)
            ax.scatter(x+rng.uniform(-.11,.11,len(display)),v[display],s=4,color=color,alpha=.30,
                edgecolors="none",rasterized=False,zorder=3)
            q25,median,q75=np.quantile(v,[.25,.5,.75])
            ax.plot([x,x],[q25,q75],color=INK,lw=2.0,zorder=4,solid_capstyle="round")
            ax.scatter([x],[median],s=15,facecolor="white",edgecolor=INK,lw=.8,zorder=5)
        label=experiment_label(tech,multiline=True)
        labels.append(f"{label}\nn = {len(v):,}")
    ax.set_xticks(range(len(techs)),labels);ax.tick_params(axis="x",length=0,pad=7)
    if 2<len(techs)<=4:
        ax.set_xticklabels(labels,rotation=23,ha="right",rotation_mode="anchor")
    ax.set_xlim(-.55,len(techs)-.45);ax.set_ylim(*ylim);ax.set_ylabel(ylabel,labelpad=5)
    if log:ih_log_ticks(ax,ylim)
    ih_axis_finish(ax)


def ih_length_lines(ax,hist,techs,ylim=None):
    for tech in techs:
        f=hist[(hist.technology==tech)&~hist.overflow.astype(bool)].sort_values("length_bp")
        x=f.length_bp.to_numpy();v=f.fraction.to_numpy()*100
                                                                              
        group=x//5;values=pd.DataFrame({"x":x,"v":v,"group":group}).groupby("group").agg({"x":"mean","v":"sum"})
        ax.plot(values.x,values.v,color=IH_COLORS[tech],lw=1.3,
            linestyle="--" if tech.endswith("_OverNight") else "-",label=experiment_label(tech))
    maximum=int(hist.loc[~hist.overflow.astype(bool),"length_bp"].max())
    ax.set_xlim(0,maximum);ax.set_xticks(np.linspace(0,maximum,6));ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Unique fragments (%) / 5 bp")
    if ylim:ax.set_ylim(*ylim)
    else:ax.set_ylim(bottom=0)
    ih_axis_finish(ax);ax.legend(frameon=False,loc="upper right",handlelength=2.1,labelspacing=.45)


def ih_single_violin(config,data,modality,metric,title,label,path):
    frame=data[(data.modality==modality)&(data.comparison=="full")]
    techs=[t for t in TECHS if t in set(frame.technology)]
    wide=len(techs)>4
    fig,ax=plt.subplots(figsize=((174 if wide else 112 if len(techs)>3 else 92)*MM,96*MM))
    ih_violin(ax,frame,techs,metric,ih_limits(frame,metric),label)
    ax.set_title(title,fontweight="bold",fontsize=10,pad=16)
    fig.subplots_adjust(left=.14 if wide else .24,right=.98,top=.85,bottom=.28)
    ih_save(fig,path,config)


def ih_fragment_curve(config,hist,mark,path):
    frame=hist[(hist.modality==mark)&(hist.comparison=="full")]
    techs=[t for t in TECHS if t in set(frame.technology)]
    fig,ax=plt.subplots(figsize=(174*MM,91*MM))
    ih_length_lines(ax,frame,techs);ax.set_title(mark+" fragment lengths",fontweight="bold",fontsize=10,pad=16)
    ax.get_legend().remove()
    ax.legend(frameon=False,loc="upper left",bbox_to_anchor=(1.015,1),handlelength=2.1,labelspacing=.65,fontsize=7)
    fig.subplots_adjust(left=.11,right=.64,top=.84,bottom=.23)
    fig.text(.11,.035,f"All QC-passing cells; normalization includes >{config['fragment_length_max']:,} bp fragments.",fontsize=6.4,color=MUTED)
    ih_save(fig,path,config)


def ih_qc_sets(qc,tech):
    frame=qc[qc.technology==tech]
    return {m:set(frame.loc[(frame.modality==m)&frame.pass_qc,"barcode"].astype(str)) for m in MODS if m in set(frame.modality)}


def ih_plot_dna_venn(config,sets,tech,path):
    fig,ax=plt.subplots(figsize=(105*MM,96*MM));left,right=(sets[m] for m in MARKS)
    for x,mark in [(-.5,MARKS[0]),(.5,MARKS[1])]:
        ax.add_patch(Circle((x,0),.94,facecolor=MOD_COLORS[mark],edgecolor=MOD_COLORS[mark],alpha=.17,lw=1.5))
        ax.text(x,.99,mark,ha="center",va="bottom",fontsize=9,fontweight="bold",color=MOD_COLORS[mark])
    for x,value in [(-.9,len(left-right)),(0,len(left&right)),(.9,len(right-left))]:
        ax.text(x,0,f"{value:,}",ha="center",va="center",fontsize=12,fontweight="bold" if x==0 else "normal")
    ax.set_xlim(-1.65,1.65);ax.set_ylim(-1.3,1.3);ax.set_aspect("equal");ax.axis("off")
    ax.set_title("In-house NanoCNT | QC-passing cells",color=IH_COLORS[tech],fontweight="bold",fontsize=10,pad=15)
    q=config["qc"]
    fig.text(.09,.06,f"H3K27ac: {q['H3K27ac_min']:,}-{q['dna_max']:,}; H3K27me3: {q['H3K27me3_min']:,}-{q['dna_max']:,} nuclear fragments.\nRNA not supplied. Circle areas are schematic; counts are exact.",fontsize=7,color=MUTED,linespacing=1.7)
    fig.subplots_adjust(left=.05,right=.95,top=.86,bottom=.23);ih_save(fig,path,config)


def pub_style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,'axes.labelsize':8,
        'xtick.labelsize':7.3,'ytick.labelsize':7.3,'text.color':INK,'axes.labelcolor':INK,
        'axes.edgecolor':'#8B969F','axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,
        'xtick.major.size':2.5,'ytick.major.size':2.5,'pdf.fonttype':42,'ps.fonttype':42,
        'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white','path.simplify':False})


def pub_violin(groups,title,ylabel,path,log=False,percent=False,preview=False,ax=None):
    pub_style();groups={t:np.asarray(v,dtype=float) for t,v in groups.items()}
    groups={t:v[np.isfinite(v)] for t,v in groups.items()}
    techs=[t for t in PUB_ORDER if t in groups and len(groups[t])]
    if not techs:return
    values=[np.asarray(groups[t],float) for t in techs]
    values=[v[np.isfinite(v)] for v in values]
    values=[np.log10(v[v>=0]+1) if log else v*100 if percent else v for v in values]
    pooled=np.concatenate(values);width=max(92,min(180,22*len(techs)+40))
    owns_figure=ax is None
    if owns_figure:fig,ax=plt.subplots(figsize=(width*MM,96*MM))
    else:fig=ax.figure
    bandwidth=max(np.std(pooled)*max(len(pooled),2)**(-.2),.025 if log else .5)
    densities=[]
    for v in values:
        if len(v)>2 and np.std(v,ddof=1)>1e-10:
            grid=np.linspace(v.min(),v.max(),150);d=gaussian_kde(v,bw_method=bandwidth/np.std(v,ddof=1))(grid)
        else:grid=np.array([]);d=np.array([])
        densities.append((grid,d))
    maximum=max([d.max() for grid,d in densities if len(d)]+[1e-12])
    for i,(tech,v,(grid,d)) in enumerate(zip(techs,values,densities)):
        if len(grid):ax.fill_betweenx(grid,i-.36*d/maximum,i+.36*d/maximum,color=PUB_COLORS[tech],alpha=.70,lw=.7,edgecolor=PUB_COLORS[tech])
        if len(v):
            rng=np.random.default_rng(731+i);show=rng.choice(len(v),min(len(v),180),replace=False)
            ax.scatter(i+rng.uniform(-.11,.11,len(show)),v[show],s=4,color=PUB_COLORS[tech],alpha=.28,linewidths=0,rasterized=False)
            a,median,b=np.quantile(v,[.25,.5,.75]);ax.plot([i,i],[a,b],color=INK,lw=2,zorder=4)
            ax.scatter([i],[median],facecolor='white',edgecolor=INK,s=15,lw=.8,zorder=5)
    ax.set_xticks(range(len(techs)),[f'{t}\nn = {len(v):,}' for t,v in zip(techs,values)],rotation=23 if len(techs)>2 else 0,ha='right' if len(techs)>2 else 'center')
    ax.tick_params(axis='x',length=0,pad=7);ax.set_xlim(-.55,len(techs)-.45)
    if len(pooled):
        low,high=(pooled.min()-.12,pooled.max()+.12) if log else (max(0,pooled.min()*.8),max(1,pooled.max()*1.12))
        if percent:low=0;high=min(100,max(10,high))
        ax.set_ylim(low,max(low+.1,high))
        if log:
            ticks=[v for v in (1,2,5,10,20,50,100,200,300,500,1000,2000,3000,5000,10000,15000,25000,50000) if low<=np.log10(v+1)<=high]
            if len(ticks)>6:ticks=ticks[::2]
            ax.set_yticks(np.log10(np.asarray(ticks)+1),[f'{v/1000:g}k' if v>=1000 else str(v) for v in ticks])
    ax.set_ylabel(ylabel);ax.set_title(title,fontsize=10,fontweight='bold',pad=16)
    ax.yaxis.grid(True,color='#EDF0F2',lw=.55);ax.set_axisbelow(True)
    if owns_figure:
        fig.subplots_adjust(left=.23 if width<115 else .14,right=.98,bottom=.30,top=.83)
        pub_save(fig,path,preview)
    return ax


def overlap(r1,r2,d):
    if min(r1,r2)<=0 or d>=r1+r2:return 0.
    if d<=abs(r1-r2):return math.pi*min(r1,r2)**2
    a=np.clip((d*d+r1*r1-r2*r2)/(2*d*r1),-1,1);b=np.clip((d*d+r2*r2-r1*r1)/(2*d*r2),-1,1)
    product_=max(0,(-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
    return r1*r1*math.acos(a)+r2*r2*math.acos(b)-.5*math.sqrt(product_)


def venn(sets,title,path,preview=False,unit='Cells'):
    pub_style();names=list(sets);n=len(names)
    if not 1<=n<=3:return
    scale=max([len(s) for s in sets.values()]+[1]);radii=np.sqrt(np.array([len(s) for s in sets.values()])/scale/math.pi)
    centers=np.zeros((n,2));distances={}
    for i in range(n):
        for j in range(i+1,n):
            wanted=len(sets[names[i]]&sets[names[j]])/scale
            if wanted<=0:d=radii[i]+radii[j]+.05
            elif wanted>=math.pi*min(radii[i],radii[j])**2-1e-12:d=abs(radii[i]-radii[j])
            else:d=brentq(lambda x:overlap(radii[i],radii[j],x)-wanted,abs(radii[i]-radii[j]),radii[i]+radii[j])
            distances[(i,j)]=d
    if n>=2:centers[1]=[distances[(0,1)],0]
    if n==3:
        d01=max(distances[(0,1)],1e-6);d02=distances[(0,2)];d12=distances[(1,2)]
        x=(d02*d02+d01*d01-d12*d12)/(2*d01);y=math.sqrt(max(0,d02*d02-x*x))
        target={(i,j):len(sets[names[i]]&sets[names[j]])/scale for i,j in distances}
        def unpack(v):return np.array([[0,0],[v[0],0],[v[1],v[2]]])
        def loss(v):
            c=unpack(v);return [overlap(radii[i],radii[j],np.linalg.norm(c[i]-c[j]))-area for (i,j),area in target.items()]
        fit=least_squares(loss,np.clip([d01,x,max(y,.001)],[0,-2,0],[2,2,2]),bounds=([0,-2,0],[2,2,2]))
        centers=unpack(fit.x)
    counts={};universe=set.union(*sets.values())
    for bits in product('01',repeat=n):
        mask=''.join(bits)
        if '1' not in mask:continue
        selected=set(universe)
        for bit,name in zip(mask,names):selected=selected&sets[name] if bit=='1' else selected-sets[name]
        counts[mask]=len(selected)
    colors={'H3K27ac':'#338F80','H3K27me3':'#B29660','RNA':'#6679A9'}
    # A side key keeps all seven exact intersections readable even for tiny crescents.
    fig,ax=plt.subplots(figsize=((180 if n==3 else 115)*MM,108*MM))
    active=radii>0
    if active.any():
        lo=np.min(centers[active]-radii[active,None],axis=0);hi=np.max(centers[active]+radii[active,None],axis=0)
        middle=(lo+hi)/2;span=max(max(hi-lo),.1)
        ax.set_xlim(middle[0]-.58*span,middle[0]+.58*span);ax.set_ylim(middle[1]-.58*span,middle[1]+.58*span)
        for i,name in enumerate(names):
            if radii[i]>0:
                ax.add_patch(Circle(centers[i],radii[i],facecolor=colors[name],edgecolor='none',alpha=.22))
                ax.add_patch(Circle(centers[i],radii[i],fill=False,edgecolor=colors[name],lw=1.4))
        xx,yy=np.meshgrid(np.linspace(lo[0]-.02,hi[0]+.02,250),np.linspace(lo[1]-.02,hi[1]+.02,250))
        ins=[(xx-c[0])**2+(yy-c[1])**2<r*r if r else np.zeros(xx.shape,bool) for c,r in zip(centers,radii)]
        for mask,count in counts.items():
            if not count:continue
            region=np.ones(xx.shape,bool)
            for bit,inside in zip(mask,ins):region &= inside if bit=='1' else ~inside
            if region.any():
                distance=distance_transform_edt(region);iy,ix=np.unravel_index(np.argmax(distance),distance.shape)
                if distance[iy,ix]>6:ax.text(xx[iy,ix],yy[iy,ix],f'{count:,}',ha='center',va='center',fontweight='bold' if mask=='1'*n else 'normal',fontsize=10,bbox={'facecolor':'white','alpha':.6,'edgecolor':'none','pad':.6})
    else:ax.text(.5,.5,'0',transform=ax.transAxes,ha='center',va='center',fontsize=16)
    ax.axis('off');ax.set_aspect('equal')
    handles=[Patch(facecolor=colors[name],alpha=.65,label=f'{name}   n = {len(sets[name]):,}') for name in names]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.29 if n==3 else .5,.02),frameon=False,fontsize=8)
    if n==3:
        ax.set_position([.02,.22,.50,.64]);panel=fig.add_axes([.56,.23,.41,.63]);panel.axis('off');panel.set_xlim(0,1);panel.set_ylim(0,1)
        panel.text(.02,.98,'Key intersections' if unit!='Cells' else 'Cell intersections',fontsize=9,fontweight='bold')
        xs=[.12,.40,.66]
        for x,name in zip(xs,names):panel.text(x,.85,name,fontsize=7.6,ha='center',color=colors[name],fontweight='bold')
        panel.text(.99,.85,'Keys' if unit!='Cells' else 'Cells',fontsize=8,ha='right',fontweight='bold')
        for row,(mask,count) in enumerate(counts.items()):
            y=.74-row*.09
            if row%2==0:panel.axhspan(y-.035,y+.035,color='#F4F7F8')
            for x,bit,name in zip(xs,mask,names):panel.scatter([x],[y],s=30,facecolors=colors[name] if bit=='1' else 'white',edgecolors=colors[name] if bit=='1' else '#B9C4CB',lw=.8)
            panel.text(.99,y,f'{count:,}',ha='right',va='center',fontsize=9,fontweight='bold' if mask=='111' else 'normal')
    else:
        fig.subplots_adjust(left=.05,right=.95,bottom=.24,top=.84)
        # Two-set region labels can be outside very thin or nested regions; list exact counts.
        if n==2:
            text=f'{names[0]} only: {counts["10"]:,}   |   Both: {counts["11"]:,}   |   {names[1]} only: {counts["01"]:,}'
            fig.text(.5,.23,text,ha='center',fontsize=7,color=MUTED)
    fig.suptitle(title,fontsize=10,fontweight='bold',y=.96);pub_save(fig,path,preview)


def k422_k9_dt_style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 10,
        "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": "#8B969F", "axes.linewidth": .65, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "legend.fontsize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white",
        "axes.facecolor": "white", "savefig.transparent": False, "path.simplify": False})


def k422_k9_dt_finish(ax):
    ax.yaxis.grid(True, color="#EDF0F2", lw=.55); ax.set_axisbelow(True)


def k422_k9_dt_count_label(value):
    return f"{value / 1000:g}k" if value >= 1000 else f"{value:g}"


def k422_k9_dt_count_ticks(ax, low, high):
    ticks = [x for x in (100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 15000, 25000, 50000) if low <= np.log10(x) <= high]
    if len(ticks) > 6:
        ticks = ticks[::2]
    ax.set_yticks(np.log10(ticks), [k422_k9_dt_count_label(x) for x in ticks])


def k422_k9_dt_violin(config, frame, metric, title, ylabel, stem, pooled=False, logarithmic=False):
    libraries = ["pooled"] if pooled else [lib for lib in config["libraries"] if lib in set(frame.library)]
    groups = [frame[metric].dropna().to_numpy(float)] if pooled else [frame.loc[frame.library == lib, metric].dropna().to_numpy(float) for lib in libraries]
    percent = metric.startswith("FRiP")
    groups = [g * 100 if percent else np.log10(g[g > 0]) if logarithmic else g for g in groups]
    all_values = np.concatenate(groups) if groups else np.array([])
    fig, ax = plt.subplots(figsize=((90 if pooled else 125) * MM, 101 * MM))
    bandwidth = max(np.std(all_values) * max(len(all_values), 2) ** (-.2), .025 if logarithmic else .5)
    densities = []
    for values in groups:
        if len(values) >= 3 and np.std(values, ddof=1) > 1e-10:
            grid = np.linspace(values.min(), values.max(), 170)
            density = gaussian_kde(values, bw_method=bandwidth / np.std(values, ddof=1))(grid)
        else:
            grid, density = np.array([]), np.array([])
        densities.append((grid, density))
    maximum = max([d.max() for g, d in densities if len(d)] + [1e-12])
    for i, (lib, values, (grid, density)) in enumerate(zip(libraries, groups, densities)):
        color = LIB_COLORS.get(lib, "#315A79")
        if len(grid):
            width = .34 * density / maximum
            ax.fill_betweenx(grid, i - width, i + width, facecolor=color, edgecolor=color, alpha=.67, lw=.7)
        if len(values):
            rng = np.random.default_rng(417 + i); visible = rng.choice(len(values), min(180, len(values)), replace=False)
            ax.scatter(i + rng.uniform(-.1, .1, len(visible)), values[visible], s=4, color=color, alpha=.25, linewidths=0, rasterized=False)
            a, median, b = np.quantile(values, [.25, .5, .75])
            ax.plot([i, i], [a, b], lw=2, color=INK, zorder=4, solid_capstyle="round")
            ax.scatter([i], [median], s=16, facecolor="white", edgecolor=INK, linewidths=.8, zorder=5)
    if len(all_values):
        low, high = (max(0, all_values.min() * .8), all_values.max() * 1.12 + 1) if not logarithmic else (all_values.min() - .12, all_values.max() + .12)
        if percent:
            low, high = 0, min(100, max(10, high))
        ax.set_ylim(low, max(low + .1, high))
        if logarithmic:
            k422_k9_dt_count_ticks(ax, low, high)
    labels = [("A + B + C" if pooled else "Library " + lib) + f"\nn = {len(values):,}" for lib, values in zip(libraries, groups)]
    ax.set_xticks(range(len(libraries)), labels); ax.tick_params(axis="x", length=0, pad=8)
    ax.set_xlim(-.6, max(.6, len(libraries) - .4)); ax.set_ylabel(ylabel + (" (log scale)" if logarithmic else ""))
    ax.set_title(title, fontweight="bold", pad=16); k422_k9_dt_finish(ax)
    fig.subplots_adjust(left=.24 if pooled else .18, right=.97, bottom=.25, top=.84)
    k422_k9_dt_save(fig, stem, config)


def k422_k9_dt_length_plot(config, frame, mark, stem, pooled):
    f = frame[frame.modality == mark]
    fig, ax = plt.subplots(figsize=(125 * MM, 91 * MM))
    libraries = ["pooled"] if pooled else list(config["libraries"])
    for lib in libraries:
        hist = f.groupby("length_bp").n_unique_fragments.sum() if pooled else f[f.library == lib].set_index("length_bp").n_unique_fragments
        total = hist.sum(); visible = hist[hist.index <= config["fragment_length_max"]]
        bins = visible.groupby(visible.index // 5).sum() / max(total, 1) * 100
        ax.plot(bins.index * 5 + 2, bins, lw=1.4, color=LIB_COLORS.get(lib, "#315A79"), label="A + B + C" if pooled else "Library " + lib)
    ax.set_xlim(0, config["fragment_length_max"]); ax.set_ylim(bottom=0)
    ax.set_xlabel("Fragment length (bp)"); ax.set_ylabel("Unique fragments (%) / 5 bp")
    ax.set_title(mark + " fragment lengths", fontweight="bold", pad=15); k422_k9_dt_finish(ax)
    if not pooled:
        ax.legend(frameon=False)
    fig.subplots_adjust(left=.18, right=.97, bottom=.23, top=.83)
    k422_k9_dt_save(fig, stem, config)


def fingerprint(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def canonicalize_barcodes(ids):
    """Remove a terminal -1 only when doing so cannot merge distinct barcodes."""
    original = [str(value) for value in ids]
    if any(not value or any(char.isspace() for char in value) for value in original):
        raise ValueError("Empty or whitespace-containing cellular barcode.")
    if len(set(original)) != len(original):
        raise ValueError("Duplicate cellular barcodes in one library/modality.")
    canonical = [value[:-2] if value.endswith("-1") else value for value in original]
    if len(set(canonical)) != len(canonical):
        raise ValueError("Removing terminal -1 would create cellular barcode collisions.")
    if any(not value for value in canonical):
        raise ValueError("A barcode cannot consist only of the terminal -1 suffix.")
    return canonical


def canonical_barcode(value):
    value = str(value)
    if not value or any(char.isspace() for char in value):
        raise ValueError("Empty or whitespace-containing BAM cellular barcode.")
    return value[:-2] if value.endswith("-1") else value


def _primary(read):
    return not read.is_secondary and not read.is_supplementary


def _base_pair_reason(read1, read2):
    if not (read1.is_paired and read2.is_paired):
        return "not_paired"
    if not (read1.is_read1 and read2.is_read2):
        return "invalid_mate_identity"
    if read1.is_unmapped or read2.is_unmapped or read1.mate_is_unmapped or read2.mate_is_unmapped:
        return "unmapped_pair"
    if read1.is_duplicate or read2.is_duplicate:
        return "duplicate_flag"
    return None


def fragment_from_pair(read1, read2, barcode_tag, min_mapq=20):
    """Return (fragment-or-None, reason); pairs are primary and name-grouped.

    A fragment is (chromosome, start, end, unnormalized_barcode). Coordinates
    are the outer aligned endpoints of both reads; no Tn5 shift is applied.
    """
    reason = _base_pair_reason(read1, read2)
    if reason:
        return None, reason
    if read1.reference_id != read2.reference_id:
        return None, "interchromosomal"
    if min(read1.mapping_quality, read2.mapping_quality) < min_mapq:
        return None, "low_mapq_either_mate"
    if not read1.has_tag(barcode_tag) or not read2.has_tag(barcode_tag):
        return None, "missing_barcode_tag"
    barcode1, barcode2 = str(read1.get_tag(barcode_tag)), str(read2.get_tag(barcode_tag))
    if barcode1 != barcode2:
        raise ValueError(f"Mates carry different {barcode_tag} cellular barcodes for {read1.query_name}.")
    canonical_barcode(barcode1)                                                        
    chrom = str(read1.reference_name)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    if chrom not in NUCLEAR:
        return None, "non_nuclear_contig"
    if read1.reference_end is None or read2.reference_end is None:
        return None, "missing_alignment_endpoint"
    start = min(read1.reference_start, read2.reference_start)
    end = max(read1.reference_end, read2.reference_end)
    if start < 0 or end <= start:
        return None, "invalid_fragment_span"
    return (chrom, int(start), int(end), barcode1), "retained"


def _detect_barcode_tag(bam_path, rna_barcodes, explicit=None, max_records=100000):
    import pysam
    if explicit:
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", str(explicit)):
            raise ValueError("barcode_tag must be an explicit two-character SAM tag.")
        return str(explicit), {"mode": "explicit", "tag": str(explicit)}
    rna_ids = set(canonicalize_barcodes(list(rna_barcodes)))
    if not rna_ids:
        raise ValueError("Cellular BAM tag autodetection requires RNA barcodes; set barcode_tag explicitly.")
    present = Counter()
    matched = Counter()
    diversity = {tag: set() for tag in TAG_CANDIDATES}
    same_as_cb = Counter()
    sampled = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if not _primary(read) or read.is_unmapped:
                continue
            sampled += 1
            for tag in TAG_CANDIDATES:
                if read.has_tag(tag):
                    value = canonical_barcode(read.get_tag(tag))
                    present[tag] += 1
                    diversity[tag].add(value)
                    matched[tag] += value in rna_ids
                    if read.has_tag("CB") and canonical_barcode(read.get_tag("CB")) == value:
                        same_as_cb[tag] += 1
            if sampled >= max_records:
                break
    audit = {tag: {"present_records": int(present[tag]), "matching_records": int(matched[tag]),
                   "matching_fraction": float(matched[tag] / max(present[tag], 1)),
                   "unique_values": len(diversity[tag])} for tag in TAG_CANDIDATES}
                                                                               
                                                                               
    min_ids = min(3, len(rna_ids))
    candidates = [tag for tag in TAG_CANDIDATES
                  if len(diversity[tag] & rna_ids) >= min_ids
                  and matched[tag] >= min(20, max(1, sampled // 10))
                  and matched[tag] / max(present[tag], 1) >= .01]
    if not candidates:
        raise ValueError("No cellular CB/CR/BC/RG tag can be verified against RNA barcodes. "
                         "Set barcode_tag explicitly after inspecting BAM tags. Audit: " + json.dumps(audit))
    candidates.sort(key=lambda tag: matched[tag], reverse=True)
    if len(candidates) > 1:
        best, second = candidates[:2]
        if matched[best] < matched[second] * 1.5:
                                                                              
                                                                       
            if best == "CB" and present[best] == present[second] == same_as_cb[second]:
                pass
            else:
                raise ValueError("Ambiguous cellular BAM tags; set barcode_tag explicitly. Audit: " + json.dumps(audit))
    return candidates[0], {"mode": "RNA_overlap", "sampled_records": sampled, "candidates": audit}


def extract_fragments(bam_path, output_gz, rna_barcodes, barcode_tag=None,
                      min_mapq=20, threads=4, tmpdir=None):
    """Name-collate one NoDup BAM and stream both-mate-verified fragments.

    Mapped unique fragments count primary nonduplicate mapped pairs before
    nuclear/MAPQ/barcode filtering; this is distinct from QC fragments. The
    retained fragments are nuclear and both mates pass MAPQ. All barcode
    identities are preserved even if no RNA cell passed upstream cell calling.
    """
    import pysam
    bam_path, output_gz = Path(bam_path), Path(output_gz)
    if not bam_path.is_file():
        raise FileNotFoundError(bam_path)
    if not shutil.which("samtools"):
        raise RuntimeError("samtools is required for bounded-memory name collation.")
    if not 0 <= min_mapq <= 255 or threads < 1:
        raise ValueError("Invalid MAPQ or thread count.")
    tag, tag_audit = _detect_barcode_tag(bam_path, rna_barcodes, barcode_tag)
    output_gz.parent.mkdir(parents=True, exist_ok=True)
    partial = output_gz.with_name(output_gz.name + ".writing")
    counts = Counter()
    exclusions = Counter()
    mapped_unique = groups = 0
    originals = {}
    try:
        with tempfile.TemporaryDirectory(prefix="k27k9-collate-", dir=tmpdir) as temp:
            log_path = Path(temp) / "samtools.log"
            prefix = Path(temp) / "collate"
            command = ["samtools", "collate", "-u", "-O", "-@", str(threads), str(bam_path), str(prefix)]
            with log_path.open("wb") as logfile:
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=logfile)
                try:
                    with pysam.AlignmentFile(process.stdout, "rb") as handle, gzip.open(partial, "wt", encoding="utf-8", compresslevel=3) as target:
                        for query_name, members in groupby(handle.fetch(until_eof=True), lambda read: read.query_name):
                            groups += 1
                            primary = [read for read in members if _primary(read)]
                            first = [read for read in primary if read.is_read1]
                            second = [read for read in primary if read.is_read2]
                            if len(first) > 1 or len(second) > 1:
                                raise ValueError(f"More than one primary alignment per mate for QNAME {query_name}; cannot identify a unique fragment.")
                            if len(first) != 1 or len(second) != 1:
                                exclusions["orphan_primary_mate"] += 1
                                continue
                            read1, read2 = first[0], second[0]
                            if _base_pair_reason(read1, read2) is None:
                                mapped_unique += 1
                            fragment, reason = fragment_from_pair(read1, read2, tag, min_mapq)
                            exclusions[reason] += 1
                            if fragment is None:
                                continue
                            chrom, start, end, original = fragment
                            canonical = canonical_barcode(original)
                            if canonical in originals and originals[canonical] != original:
                                raise ValueError(f"BAM barcode normalization collision: {originals[canonical]!r} and {original!r}.")
                            originals[canonical] = original
                            target.write(f"{chrom}\t{start}\t{end}\t{canonical}\n")
                            counts[canonical] += 1
                    if process.stdout:
                        process.stdout.close()
                    return_code = process.wait()
                    if return_code:
                        raise RuntimeError(f"samtools collate failed ({return_code}): {log_path.read_text(errors='replace')[-4000:]}")
                except BaseException:
                    if process.poll() is None:
                        process.terminate()
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    raise
        partial.replace(output_gz)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    frame = pd.DataFrame({"unique_fragments": [counts[key] for key in sorted(counts)]},
                         index=pd.Index(sorted(counts), name="barcode"))
    frame["unique_fragments"] = frame.unique_fragments.astype(np.int64)
    metadata = {"source": fingerprint(bam_path), "barcode_tag": tag, "tag_audit": tag_audit,
                "query_name_groups": groups, "unique_mapped_fragments": mapped_unique,
                "nuclear_mapq_fragments": int(sum(counts.values())), "min_mapq_both_mates": min_mapq,
                "exclusions_and_retained": dict(exclusions), "fragment_coordinates": "outer_aligned_endpoints_no_Tn5_shift",
                "mapped_count_definition": "primary paired both mapped neither duplicate, before benchmark MAPQ/contig/barcode filters"}
    return frame, metadata


def parse_ids(values,library=None):
    rows=[]
    for value in map(str,values):
        match=RNA_ID.fullmatch(value.strip())
        if match is None:raise ValueError(f'Unrecognized original MTR RNA barcode: {value!r}; retain RT and library metadata.')
        prefix,barcode,rt=match.groups()
        supplied=LIBRARIES.get(prefix,prefix) if prefix else library
        if supplied not in LIBRARIES.values():raise ValueError(f'Missing/unknown MTR sublibrary for {value!r}')
        if library is not None and supplied!=library:raise ValueError(f'Conflicting MTR sublibrary for {value!r}: {supplied} versus {library}')
        rows.append((barcode+':'+supplied,rt,supplied))
    return pd.DataFrame(rows,columns=['physical_id','rt_code','sublibrary'])


def assign_timepoints(counts,min_fraction=.5):
    grouped=counts.groupby(['physical_id','rt_code'],sort=True).rna_matrix_counts.sum().unstack(fill_value=0)
    total=grouped.sum(axis=1)
    dominant=grouped.idxmax(axis=1)
    top=grouped.max(axis=1)
    fraction=top.div(total.where(total>0)).fillna(0)
    assigned=(total>0)&(fraction>min_fraction)
    timepoint=dominant.map(TIMEPOINTS).where(assigned,'unassigned').fillna('unassigned')
    frame=pd.DataFrame({'physical_id':grouped.index,'dominant_rt':dominant.to_numpy(),
                        'rna_matrix_counts_all_rt':total.to_numpy(),'dominant_rt_counts':top.to_numpy(),
                        'dominant_rt_fraction':fraction.to_numpy(),'timepoint':timepoint.to_numpy()})
    frame['assignment_status']=np.where(total.to_numpy()==0,'no_RNA_counts',
        np.where(~assigned.to_numpy(),'no_RT_majority',np.where(timepoint.to_numpy()=='unassigned','unknown_RT','assigned')))
    for rt in grouped.columns:frame['rt_'+rt+'_counts']=grouped[rt].to_numpy()
    return frame


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def save_pdf(fig, path, pad, creator=None):
    path = Path(path).with_suffix(".pdf")
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"Creator": creator, "Title": path.stem.replace("_", " ")} if creator else None
    fig.savefig(path, format="pdf", bbox_inches="tight", pad_inches=pad, metadata=metadata)
    plt.close(fig)


def ih_save(fig, path, config):
    save_pdf(fig, path, .055, "K422 focused benchmark")


def pub_save(fig, path, preview=False):
    save_pdf(fig, path, .07)


def k422_k9_dt_save(fig, path, config):
    save_pdf(fig, path, .08)


def assemble_mtr(mark, work, day0):
    """Pool published fragments for the Day 0 cells assigned from GEO RNA counts."""
    work = Path(work)
    target = work / f'MTR_{mark}.fragments.tsv'
    if target.is_file(): return target
    if not day0: raise ValueError('No Day 0 cells in the GEO RNA matrices')
    cells = work / 'MTR_day0_ids.txt'
    if not cells.is_file(): cells.write_text('\n'.join(sorted(day0)) + '\n')
    pending = target.with_suffix('.pending')
    # Keep all five fields, then restore cell/chromosome/coordinate order within each sublibrary.
    program = 'NR==FNR{keep[$1]=1;next} ($4 in keep){print $0}'
    with pending.open('wb') as out:
        for i in range(1, 9):
            raw = work / f'MTR_{mark}_S{i}.selected.tsv'
            with gzip.open(input_path(f'MTR_S{i}_{mark}'), 'rb') as source, raw.open('wb') as filtered:
                process = subprocess.Popen(['awk', '-F', '\t', program, str(cells), '-'],
                                           stdin=subprocess.PIPE, stdout=filtered)
                try:
                    shutil.copyfileobj(source, process.stdin, 1024*1024)
                    process.stdin.close()
                    if process.wait(): raise RuntimeError('MTR assembly cell filtering failed')
                finally:
                    if process.poll() is None: process.kill(); process.wait()
            subprocess.run(['sort', '-s', '-t', '\t', '-k4,4', '-k1,1', '-k2,2n', '-k3,3n',
                            '-S', '1G', '--parallel=2', '-T', str(work), str(raw)],
                           stdout=out, check=True, env=dict(os.environ, LC_ALL='C'))
            raw.unlink()
    pending.replace(target)
    return target


def input_path(key, replacements=None):
    path = INPUTS[key] if replacements is None else replacements.get(key, INPUTS[key])
    return path.expanduser().resolve()


def dataset_inputs(root):
    replacements = {}
    if root is not None:
        for dataset, source in INHOUSE.items():
            for mark, key in source["dna"].items():
                replacements[key] = fragment_path(root, dataset, mark)
            if "rna" in source:
                matrix, features, barcodes = rna_paths(root, dataset)
                replacements.update(zip(source["rna"], (matrix, barcodes, features)))
    return replacements


def mex(keys, unique_names=True, replacements=None):
    return read_mex((input_path(key, replacements) for key in keys), unique_names)


def read_mex(paths, unique_names=True):
    import anndata as ad
    matrix, barcodes, features = paths
    with open_text(matrix) as handle:
        x = sparse.csr_matrix(mmread(handle)).T.tocsr()
    with open_text(barcodes) as handle:
        ids = [line.rstrip("\r\n").split("\t")[0] for line in handle]
    with open_text(features) as handle:
        rows = [line.rstrip("\r\n").split("\t") for line in handle]
    names = [row[1] if len(row) > 1 and row[1] else row[0] for row in rows]
    if x.shape != (len(ids), len(names)):
        raise ValueError(f"MEX dimensions disagree: {matrix}")
    x.sum_duplicates(); x.eliminate_zeros()
    if not np.isfinite(x.data).all() or np.any(x.data < 0) or np.any(x.data != np.rint(x.data)):
        raise ValueError(f"Expected integer RNA counts: {matrix}")
    A = ad.AnnData(x, obs=pd.DataFrame(index=ids), var=pd.DataFrame(index=names))
    if unique_names:
        A.var_names_make_unique()
    return A


def rna_frame(A, scope):
    import scanpy as sc
    names = pd.Index(A.var_names.astype(str)).str.upper()
    if scope == "show":
        x = sparse.csr_matrix(A.X); x.sum_duplicates(); x.eliminate_zeros()
        total = np.asarray(x.sum(axis=1)).ravel()
        hb = names.str.match(r"^HB(?:A[12]?|B|D|E1?|G[12]?|M|Q1?|Z)$")
        mt = names.str.startswith("MT-")
        top = np.array([np.partition(v, len(v)-50)[-50:].sum() if len(v)>50 else v.sum()
                        for v in (x.data[x.indptr[i]:x.indptr[i+1]] for i in range(x.shape[0]))])
        f = pd.DataFrame({"UMIs": total, "genes": np.diff(x.indptr),
                          "pct_mt": np.asarray(x[:,mt].sum(axis=1)).ravel()*100/np.maximum(total,1),
                          "pct_hb": np.asarray(x[:,hb].sum(axis=1)).ravel()*100/np.maximum(total,1),
                          "pct_top50": top*100/np.maximum(total,1)})
    else:
        A.var["mt"] = names.str.startswith("MT-")
        A.var["hb"] = names.str.match(r"^HB[AB]" if scope == "public" else r"^HB[ABDEGQZ]")
        sc.pp.calculate_qc_metrics(A, qc_vars=["mt","hb"], percent_top=[50], inplace=True)
        f = A.obs[["total_counts","n_genes_by_counts","pct_counts_mt","pct_counts_hb","pct_counts_in_top_50_genes"]].reset_index(drop=True)
        f.columns = ["UMIs", "genes", "pct_mt", "pct_hb", "pct_top50"]
    f["source_barcode"] = pd.Index(A.obs_names).astype(str).to_numpy()
    f["barcode"] = f.source_barcode
    f["pass_qc"] = (f.UMIs.between(1000,25000) & f.pct_mt.le(10) & f.pct_hb.le(5) & f.pct_top50.le(30))
    return f


def normalized_ids(values):
    values = list(map(str, values))
    result = [re.sub(r"-1$", "", x) for x in values]
    if len(set(result)) != len(result):
        raise ValueError("Terminal -1 normalization would merge cells")
    return result


def histone_fragments(key, rna_ids, work, memo, replacements=None):
    source = input_path(key, replacements)
    if source.suffix != ".bam":
        return source
    if key in memo:
        return memo[key]
    import snapatac2 as snap
    folder = work / key; folder.mkdir()
    path = folder / "fragments.bed.gz"
    tag, _ = detect_bam_tag(source, rna_ids)
    snap.pp.make_fragment_file(str(source), str(path), barcode_tag=tag,
                              shift_left=0, shift_right=0, min_mapq=20,
                              chunk_size=5000000, tempdir=str(folder))
    memo[key] = path
    return path


# Public: native FRAG, q=.001. In-house: blacklist-filtered BEDPE, q=.01.
def own_peaks(bed, mark, folder, public=False, mouse=False):
    folder.mkdir(parents=True, exist_ok=True)
    broad = mark != "H3K27ac"
    cmd = ["macs3", "callpeak", "-t", str(bed), "-f", "FRAG" if public else "BEDPE",
           "-g", "mm" if mouse else "hs", "--keep-dup", "0" if public else "all",
           "-q", ".001" if public else ".01", "--nolambda", "-n", "cohort", "--outdir", str(folder)]
    if public:
        cmd += ["--nomodel", "--seed", "99", "--verbose", "1"]
    if broad:
        cmd += ["--broad"]
        if not public:
            cmd += ["--broad-cutoff", ".1"]
    with (folder / "macs3.log").open("w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    regions = Regions(folder / ("cohort_peaks." + ("broadPeak" if broad else "narrowPeak")))
    if not public:
        blacklist = Regions(input_path("hg38_blacklist"))
        regions = Regions(records=[(c,a,b) for c,intervals in regions.data.items() for a,b in intervals
                                  if c in AUTOSOMES and not blacklist.overlaps(c,a,b)])
    merged = folder / "regions.bed"; regions.write(merged)
    return merged


def selected_fragments(source, selected, path, peak_input, prefix=""):
    blacklist = Regions(input_path("hg38_blacklist"))
    mapping = dict(zip(selected.source_barcode, selected.barcode))
    with gzip.open(path,"wt") as out, gzip.open(peak_input,"wt") as call:
        for c,a,b,barcode in fragments(source):
            if c in NUCLEAR and barcode in mapping:
                out.write(f"{c}\t{a}\t{b}\t{prefix}{mapping[barcode]}\t1\n")
                if c in AUTOSOMES and not blacklist.overlaps(c,a,b):
                    call.write(f"{c}\t{a}\t{b}\n")


def score_fragments(path, frame, peak, library=None):
    regions = Regions(peak)
    counts = {barcode:[0,0] for barcode in frame.barcode}
    hist = np.zeros(1002, dtype=np.int64)
    for c,a,b,barcode in fragments(path):
        counts[barcode][0] += 1
        counts[barcode][1] += regions.overlaps(c,a,b)
        hist[min(b-a,1001)] += 1
    observed = np.array([counts[b][0] for b in frame.barcode])
    if not np.array_equal(observed, frame.unique_fragments.to_numpy()):
        raise ValueError("FRiP denominators disagree with QC counts")
    frame = frame.copy()
    frame["FRiP_own"] = np.array([counts[b][1] for b in frame.barcode])/observed
    lengths = pd.DataFrame({"length_bp":np.arange(1002), "n_unique_fragments":hist,
                            "fraction":hist/max(hist.sum(),1), "overflow":np.arange(1002)==1001})
    lengths["modality"] = frame.modality.iloc[0]
    if library:
        lengths["library"] = library
    else:
        lengths["technology"] = frame.technology.iloc[0]
        lengths["comparison"] = "full"
    return frame, lengths


def compute_inhouse(work, memo, replacements=None):
    qc_frames=[]; selected=[]; histograms=[]
    for tech in TECHS:
        source = INHOUSE[tech]
        LOG.info("In-house | %s", tech)
        rna = rna_frame(mex(source["rna"], replacements=replacements), "inhouse") if "rna" in source else None
        dna={}; paths={}
        for mark,key in source["dna"].items():
            path = (histone_fragments(key,rna.barcode,work,memo,replacements)
                    if "rna" in source else input_path(key,replacements))
            paths[mark] = path
            counts=Counter(b for c,a,z,b in fragments(path) if c in NUCLEAR)
            dna[mark] = pd.DataFrame(sorted(counts.items()),columns=["source_barcode","unique_fragments"])
        rna,dna=align_barcodes(rna if rna is not None else pd.DataFrame({"barcode":pd.Series(dtype=str)}),dna)
        tables=apply_qc(rna if "rna" in source else None,dna,QC)
        for mod,frame in tables.items():
            frame["technology"]=tech; frame["modality"]=mod
            qc_frames.append(frame)
            chosen=frame[frame.pass_qc].sort_values("barcode").copy()
            chosen["comparison"]="full"
            if mod=="RNA":
                selected.append(chosen); continue
            folder=work/(tech+"_"+mod); folder.mkdir()
            bed=folder/"selected.gz"; call=folder/"calling.gz"
            selected_fragments(paths[mod],chosen,bed,call)
            peak=own_peaks(call,mod,folder/"peaks")
            metrics,hist=score_fragments(bed,chosen,peak)
            selected.append(metrics); histograms.append(hist)
    return (pd.concat(qc_frames,ignore_index=True),
            pd.concat(selected,ignore_index=True) if selected else pd.DataFrame(),
            pd.concat(histograms,ignore_index=True) if histograms else pd.DataFrame())


def build_mtr_rna(path):
    """Build Day 0 RNA counts and retain the full set of Day 0 IDs for DNA."""
    import anndata as ad

    libraries = []
    day0 = set()
    for lib, keys in MTR_RT.items():
        LOG.info("Public MTR | all-RT assignment | %s", lib)
        source = mex([keys[k] for k in ("matrix", "barcodes", "features")], unique_names=False)
        ids = parse_ids(source.obs_names, lib)
        if ids.duplicated(["physical_id", "rt_code"]).any():
            raise ValueError(f"Duplicate MTR physical ID and RT code in {lib}")
        with open_text(input_path(keys["features"])) as handle:
            features = [line.rstrip("\r\n").split("\t") for line in handle]
        rna = np.array([len(row) < 3 or bool(re.search(r"Gene Expression|RNA", row[2], re.I))
                        for row in features])
        if not rna.any():
            raise ValueError(f"No gene-expression features in MTR {lib}")
        genes = pd.DataFrame({"gene_id": [row[0] for row, keep in zip(features, rna) if keep]},
                             index=source.var_names[rna].rename("gene_symbol"))
        if not genes.index.is_unique or not genes.gene_id.is_unique:
            raise ValueError(f"Duplicate RNA gene ID or symbol in MTR {lib}")
        counts = source.X if rna.all() else source.X[:, rna]
        # GEO supplies gene counts by RT, so use all columns before RNA filtering.
        ids["rna_matrix_counts"] = np.asarray(counts.sum(axis=1, dtype=np.int64)).ravel()
        assignments = assign_timepoints(ids).set_index("physical_id")
        assigned_d0 = assignments.index[assignments.timepoint == "D0"]
        day0.update(assigned_d0)

        # Match the former Seurat preparation across all RT columns in this library.
        enough_genes = np.diff(counts.indptr) >= MTR_MIN_FEATURES
        keep_genes = counts[enough_genes].getnnz(axis=0) >= MTR_MIN_CELLS_PER_GENE
        rows = np.flatnonzero(enough_genes & ids.rt_code.eq("01") & ids.physical_id.isin(assigned_d0))
        cell_ids = pd.Index(ids.physical_id.iloc[rows], name="cell_id")
        rt = assignments.loc[cell_ids]
        obs = pd.DataFrame({
            "source_barcode": source.obs_names.to_numpy()[rows],
            "geo_accession": SOURCES[keys["matrix"]],
            "mtr_sublibrary": lib,
            "mtr_rt_code": "01",
            "mtr_timepoint": "D0",
            "mtr_rt01_gene_counts": rt.dominant_rt_counts.to_numpy(),
            "mtr_gene_counts_all_rt": rt.rna_matrix_counts_all_rt.to_numpy(),
            "mtr_rt01_fraction": rt.dominant_rt_fraction.to_numpy(),
        }, index=cell_ids)
        libraries.append(ad.AnnData(counts[rows][:, keep_genes], obs=obs,
                                    var=genes.iloc[keep_genes].copy()))
        del source, counts

    genes = pd.concat([library.var for library in libraries]).reset_index().drop_duplicates()
    if not genes.gene_symbol.is_unique or not genes.gene_id.is_unique:
        raise ValueError("Inconsistent MTR gene symbols and Ensembl IDs across libraries")
    genes = genes.set_index("gene_symbol")
    rna = ad.concat(libraries, join="outer", fill_value=0)
    if rna.n_obs == 0 or not rna.obs_names.is_unique:
        raise ValueError("Expected a nonempty MTR RNA cohort with unique Day 0 physical IDs")
    rna = rna[:, genes.index].copy()
    rna.var = genes
    rna.uns["mtr"] = {
        "geo_series": "GSE277778",
        "rt_timepoints": TIMEPOINTS,
        "rt_assignment_basis": "Gene-expression column totals across all RT codes before filtering",
        "rt_fraction_strictly_greater_than": 0.5,
        "n_day0_ids_for_dna": len(day0),
        "min_features_per_rna_column": MTR_MIN_FEATURES,
        "min_columns_per_gene": MTR_MIN_CELLS_PER_GENE,
    }
    tmp = path.with_suffix(".tmp.h5ad")
    try:
        rna.write_h5ad(tmp, compression="gzip")
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)
    LOG.info("Wrote %s | %d RNA cells, %d genes | %d Day 0 IDs for DNA",
             path, rna.n_obs, rna.n_vars, len(day0))
    return rna, day0


def public_sanitize(A,tech):
    if tech=="MTR" and not A.obs_names.is_unique:
        raise ValueError("Duplicate MTR physical IDs")
    A.obs_names=pd.Index(A.obs_names).astype(str).str.strip()
    A=A[np.flatnonzero((A.obs_names!="") & ~A.obs_names.duplicated()),:].copy()
    if tech=="TrES":
        A.obs_names=normalized_ids(A.obs_names)
    return A


# Native SnapATAC2 totals intentionally differ from the nuclear-row FRiP denominators above.
def compute_public(work, memo, out, replacements=None):
    import scanpy as sc
    import snapatac2 as snap
    mtr_rna, day0 = build_mtr_rna(out / "scRNA_D0.h5ad")
    qc=[]; frip=[]
    for tech,source in PUBLIC.items():
        LOG.info("Public | %s",tech)
        if tech == "MTR" or "rna" in source:
            if tech == "MTR":
                A = mtr_rna
                del mtr_rna
            else:
                key = source["rna"]
                A = mex(key, replacements=replacements) if isinstance(key, list) else sc.read_10x_h5(input_path(key))
            rna_ids=list(A.obs_names)
            if "feature_types" in A.var:
                keep=A.var.feature_types.astype(str).str.contains(r"Gene Expression|RNA",case=False,regex=True)
                if keep.any(): A=A[:,keep.to_numpy()].copy()
            f=rna_frame(A,"public")
            A.obs=f.set_axis(A.obs_names)
            A=public_sanitize(A[np.flatnonzero(f.pass_qc),:].copy(),tech)
            f=A.obs.copy(); f["barcode"]=A.obs_names; f["technology"]=tech; f["modality"]="RNA"
            f["n_genes"]=f.genes
            qc.append(f.reset_index(drop=True))
        for mark,key in source.get("dna",{}).items():
            path=(histone_fragments(key,rna_ids,work,memo,replacements) if tech=="TrES" else
                  assemble_mtr(mark,work,day0) if tech=="MTR" else input_path(key))
            genome=snap.genome.mm10 if tech=="NanoCNT" else snap.genome.hg38
            A=snap.pp.import_fragments(str(path),chrom_sizes=genome.chrom_sizes,min_num_fragments=0,
                sorted_by_barcode=False,chrM=["chrM","M"],tempdir=str(work),backend="hdf5",n_jobs=THREADS)
            if mark=="H3K27ac":
                snap.metrics.tsse(A,gene_anno=genome,inplace=True,n_jobs=THREADS)
            if tech=="MTR" and (not A.obs_names.is_unique or not set(A.obs_names).issubset(day0)):
                raise ValueError("MTR DNA needs unique Day 0 physical IDs")
            counts=A.obs.n_fragment
            A=A[np.flatnonzero(counts.between(QC[mark+"_min"],15000)),:].copy()
            A.obs["source_barcode"]=pd.Index(A.obs_names).astype(str)
            A=public_sanitize(A,tech)
            f=A.obs.copy();f["barcode"]=A.obs_names;f["technology"]=tech;f["modality"]=mark;f["pass_qc"]=True
            qc.append(f.reset_index(drop=True))
            folder=work/("public_"+tech+"_"+mark);folder.mkdir()
            A.obs["__group"]="all"
            files=snap.ex.export_fragments(A,groupby="__group",out_dir=str(folder),prefix="selected_",suffix=".bed.gz",compression="gzip",compression_level=1)
            bed=Path(files["all"])
            peak=own_peaks(bed,mark,folder/"peaks",public=True,mouse=tech=="NanoCNT")
            bounds={}
            for c,a,b,barcode in fragments(bed): bounds[c]=max(bounds.get(c,0),b)
            regions={"FRiP_total_fragments":[f"{c}:0-{b}" for c,b in sorted(bounds.items())]}
            if Regions(peak).bases: regions["fragments_in_own_peaks"]=peak
            snap.metrics.frip(A,regions=regions,normalized=False,count_as_insertion=False,inplace=True,n_jobs=THREADS)
            hit=A.obs.get("fragments_in_own_peaks",pd.Series(np.zeros(A.n_obs)))
            total=A.obs.FRiP_total_fragments.to_numpy()
            if np.any(total<=0) or np.any(hit.to_numpy()>total): raise ValueError("Invalid native FRiP counts")
            f["FRiP_own"]=hit.to_numpy()/total
            frip.append(f.reset_index(drop=True))
    return pd.concat(qc,ignore_index=True),pd.concat(frip,ignore_index=True)


def merged_showcase_inputs(root, work):
    """Recover the original three library groups for QC and pooled peak calling."""
    pattern = library_pattern(K422_K9_LIBRARIES)
    matrix, features, barcodes = rna_paths(root, K422_K9_DATASET)
    rna = read_mex((matrix, barcodes, features), unique_names=False)
    cells = [split_cell_id(cell, pattern) for cell in rna.obs_names]
    rna_by_library = {}
    for name, lib in K422_K9_LIBRARIES.items():
        positions = [i for i, (library, _) in enumerate(cells) if library == name]
        if not positions:
            raise ValueError(f"{K422_K9_DATASET}: no RNA cells from {name}")
        part = rna[positions].copy()
        part.obs_names = [cells[i][1] for i in positions]
        rna_by_library[lib] = part

    dna = {}
    for mark in K422_K9_DT_MARKS:
        source = fragment_path(root, K422_K9_DATASET, mark)
        paths = {lib: work / f"k422_k9_dt_{lib}_{mark}.gz" for lib in K422_K9_DT}
        originals = {lib: {} for lib in paths}
        counts = {lib: Counter() for lib in paths}
        with ExitStack() as stack:
            targets = {lib: stack.enter_context(gzip.open(path, "wt", encoding="utf-8"))
                       for lib, path in paths.items()}
            for chrom, start, end, cell in fragments(source):
                if chrom not in NUCLEAR:
                    continue
                library, original = split_cell_id(cell, pattern)
                lib = K422_K9_LIBRARIES[library]
                barcode = canonical_barcode(original)
                if originals[lib].setdefault(barcode, original) != original:
                    raise ValueError(f"{source}: barcode normalization collision in {library}: {barcode!r}")
                targets[lib].write(f"{chrom}\t{start}\t{end}\t{barcode}\n")
                counts[lib][barcode] += 1
        for lib, path in paths.items():
            frame = pd.Series(counts[lib], dtype=np.int64, name="unique_fragments")
            dna[lib, mark] = path, frame.rename_axis("barcode").sort_index().to_frame()
    return rna_by_library, dna


def k422_k9_dt_compute(work, dataset_root=None):
    if dataset_root is not None:
        merged_rna, merged_dna = merged_showcase_inputs(dataset_root, work)
    qc=[];paths={};selected={}
    for lib,source in K422_K9_DT.items():
        LOG.info("Showcase | library %s",lib)
        library_rna = mex(source["rna"],unique_names=False) if dataset_root is None else merged_rna[lib]
        rna=rna_frame(library_rna,"show")
        rna["barcode"]=canonicalize_barcodes(rna.source_barcode)
        rna["library"]=lib;rna["modality"]="RNA";qc.append(rna[rna.pass_qc].copy())
        for mark,key in source["dna"].items():
            if dataset_root is None:
                path=work/f"k422_k9_dt_{lib}_{mark}.gz"
                f, _ = extract_fragments(input_path(key),path,list(rna.barcode),threads=THREADS,tmpdir=work)
            else:
                path, f = merged_dna[lib,mark]
            f=f.reset_index();f["source_barcode"]=f.barcode;f["library"]=lib;f["modality"]=mark
            f=f[f.unique_fragments.between(300,15000)].copy()
            f["cell_id"]=lib+"::"+f.barcode
            paths[lib,mark]=path;selected[lib,mark]=f
    metrics=list(qc);histograms=[]
    # Pooled means library-tagged cells concatenated; one pooled broad call per mark.
    for mark in K422_K9_DT_MARKS:
        folder=work/("pooled_"+mark);folder.mkdir()
        calls=[];beds={}
        for lib in K422_K9_DT:
            bed=folder/(lib+".selected.gz");call=folder/(lib+".calling.gz")
            selected_fragments(paths[lib,mark],selected[lib,mark],bed,call)
            beds[lib]=bed;calls.append(call)
        combined=folder/"pooled.calling.gz"
        with combined.open("wb") as out:
            for call in calls:
                with call.open("rb") as handle: shutil.copyfileobj(handle,out)
        peak=own_peaks(combined,mark,folder/"peaks")
        for lib in K422_K9_DT:
            f,hist=score_fragments(beds[lib],selected[lib,mark],peak,library=lib)
            f["FRiP_pooled"]=f.pop("FRiP_own");metrics.append(f);histograms.append(hist)
    return pd.concat(metrics,ignore_index=True),pd.concat(histograms,ignore_index=True)


def plot_inhouse_venns(config, qc, out):
    ih_style()
    folder = Path(out)
    folder.mkdir(parents=True, exist_ok=True)
    q = config['qc']
    for tech in TECHS:
        sets = ih_qc_sets(qc, tech)
        if not sets:
            continue
        if 'RNA' not in sets:
            ih_plot_dna_venn(config, sets, tech, folder / ('QC_Venn_' + tech + '.pdf'))
            continue
        universe = set.union(*sets.values())
        fig, ax = plt.subplots(figsize=(100 * MM, 123 * MM))
        for (cx, cy), mod in zip([(-0.48, 0.3), (0.48, 0.3), (0, -0.48)], MODS):
            ax.add_patch(Circle((cx, cy), 0.9, facecolor=MOD_COLORS[mod], edgecolor=MOD_COLORS[mod], alpha=0.17, lw=1.5))
        positions = {'100': (-0.86, 0.43), '010': (0.86, 0.43), '001': (0, -0.99), '110': (0, 0.71), '101': (-0.43, -0.37), '011': (0.43, -0.37), '111': (0, -0.035)}
        for mask, xy in positions.items():
            subset = set(universe)
            for bit, mod in zip(mask, MODS):
                subset = subset & sets[mod] if bit == '1' else subset - sets[mod]
            value = len(subset)
            ax.text(*xy, f'{value:,}', ha='center', va='center', fontsize=12, fontweight='bold' if mask == '111' else 'normal', color=INK if value else '#91A0AB')
        for text, x, y, color in [('H3K27ac', -0.83, 1.39, MOD_COLORS[MARKS[0]]), ('H3K27me3', 0.83, 1.39, MOD_COLORS[MARKS[1]]), ('RNA', 0, -1.57, MOD_COLORS['RNA'])]:
            ax.text(x, y, text, ha='center', va='center', fontsize=9, color=color, fontweight='bold')
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.8, 1.65)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(experiment_label(tech) + '\nQC-passing cells', color=IH_COLORS[tech], fontweight='bold', fontsize=10, pad=14)
        loaded = set(qc.loc[qc.technology == tech, 'barcode'].astype(str))
        common = len(set.intersection(*sets.values()))
        fig.text(0.5, 0.31, f'All three: {common:,} cells | Loaded union: {len(loaded):,}', ha='center', fontsize=8, color=INK)
        threshold_text = f"H3K27ac: {q['H3K27ac_min']:,}-{q['dna_max']:,} unique fragments\nH3K27me3: {q['H3K27me3_min']:,}-{q['dna_max']:,} unique fragments\nRNA: {q['rna_min_umi']:,}-{q['rna_max_umi']:,} UMIs\nMitochondrial <= {q['max_pct_mt']:g}%; hemoglobin <= {q['max_pct_hb']:g}%\nTop 50 genes <= {q['max_pct_top50']:g}%\nCircle areas are schematic; counts are exact."
        fig.text(0.1, 0.07, threshold_text, fontsize=7.3, color=MUTED, linespacing=1.7)
        fig.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.37)
        ih_save(fig, folder / ('QC_Venn_' + tech + '.pdf'), config)


def plot_inhouse(qc, data, hist, out):
    plot_inhouse_venns(CONFIG,qc,Path(out)/GROUP_IH)
    ih_style();folder=Path(out)/GROUP_IH
    for metric,title,label,name in [("UMIs","RNA molecules","UMIs per cell (log scale)","RNA_UMIs"),
                                    ("genes","Detected genes","Genes per cell (log scale)","RNA_detected_genes")]:
        ih_single_violin(CONFIG,data,"RNA",metric,title,label,folder/(name+".pdf"))
    for mark in MARKS:
        ih_single_violin(CONFIG,data,mark,"unique_fragments",mark+" unique fragments","Unique fragments per cell (log scale)",folder/(mark+"_unique_fragments.pdf"))
        ih_fragment_curve(CONFIG,hist,mark,folder/(mark+"_fragment_lengths.pdf"))
        ih_single_violin(CONFIG,data,mark,"FRiP_own",mark+" own-peak enrichment","Fragments in own peaks (%)",folder/"supplementary"/(mark+"_own_peak_FRiP.pdf"))


def plot_public(qc, frip, out):
    folder=Path(out)/GROUP_PUBLIC
    for mark in MARKS:
        f=qc[qc.modality==mark]
        groups={tech:g.n_fragment.to_numpy() for tech,g in f.groupby("technology",sort=False)}
        pub_violin(groups,mark+" unique fragments | full","Unique fragments per cell (log scale)",folder/("Benchmark_Combined_"+mark+"_Fragments_Violin_full.pdf"),log=True)
    f=qc[qc.modality=="H3K27ac"]
    groups={tech:pd.to_numeric(g.tsse,errors="coerce").dropna().to_numpy() for tech,g in f.groupby("technology",sort=False)}
    pub_violin(groups,"H3K27ac TSS enrichment | full","TSS enrichment",folder/"Benchmark_Combined_H3K27ac_TSSE_Violin_full.pdf")
    for metric,label in [("UMIs","UMIs"),("n_genes","Detected genes")]:
        f=qc[qc.modality=="RNA"]
        groups={tech:g[metric].to_numpy() for tech,g in f.groupby("technology",sort=False)}
        pub_violin(groups,"RNA "+label.lower()+" | full",label+" per cell (log scale)",folder/("Benchmark_RNA_"+metric+"_Violin_full.pdf"),log=True)
    for tech in ("TrES","MTR","NanoCNT","NTT"):
        f=qc[qc.technology==tech]
        order=("RNA",*MARKS) if tech in ("TrES","MTR") else MARKS
        sets={mod:set(f.loc[f.modality==mod,"barcode"].astype(str)) for mod in order}
        title=tech+(" | Day 0 QC-passing cells" if tech=="MTR" else " | QC-passing cells")
        venn(sets,title,folder/("Benchmark_"+tech+"_QC_Venn.pdf"))
    pub_style();fig,axes=plt.subplots(2,1,figsize=(180*MM,170*MM))
    for ax,mark in zip(axes,MARKS):
        f=frip[frip.modality==mark]
        groups={tech:g.FRiP_own.to_numpy() for tech,g in f.groupby("technology",sort=False)}
        pub_violin(groups,mark+" | full QC cells, own peaks","Fragments in own peaks (%)",None,percent=True,ax=ax)
    fig.subplots_adjust(left=.10,right=.98,bottom=.17,top=.94,hspace=.73)
    pub_save(fig,Path(out)/GROUP_FRIP/"FRiPBench_BothMarks_full_own_FRiP_Violin.pdf")


def k422_k9_dt_plot(metrics, hist, out):
    k422_k9_dt_style();folder=Path(out)/GROUP_K422_K9_DT/"pooled"
    for metric,title,label,stem in [("UMIs","RNA molecules","UMIs per cell","RNA_UMIs"),
                                    ("genes","Detected genes","Genes per cell","RNA_detected_genes")]:
        k422_k9_dt_violin(CONFIG,metrics[metrics.modality=="RNA"],metric,title,label,folder/stem,True,True)
    for mark in K422_K9_DT_MARKS:
        f=metrics[metrics.modality==mark]
        k422_k9_dt_violin(CONFIG,f,"unique_fragments",mark+" unique fragments","Unique fragments per cell",folder/(mark+"_unique_fragments"),True,True)
        k422_k9_dt_violin(CONFIG,f,"FRiP_pooled",mark+" own-peak enrichment","Fragments in own peaks (%)",folder/(mark+"_own_peak_FRiP"),True)
        k422_k9_dt_length_plot(CONFIG,hist,mark,folder/(mark+"_fragment_lengths"),True)


def public_input_keys():
    keys = set(MTR_RNA_INPUTS)
    for tech, source in PUBLIC.items():
        for key in source.get("dna", {}).values():
            if isinstance(key, tuple):
                keys.update(key)
            else:
                keys.add(key)
                if tech == "TrES":
                    keys.add(key + "_index")
        rna = source.get("rna", [])
        keys.update([rna] if isinstance(rna, str) else rna)
    return keys


def check_inputs(required=None, dataset_root=None, replacements=None):
    if set(INPUTS) != set(SOURCES): raise ValueError("INPUTS and SOURCES keys differ")
    keys = set(INPUTS) if required is None else set(required)
    paths = []
    if dataset_root is not None:
        showcase_keys = {key for key in keys if key.startswith("K422_K9_DT_")}
        if showcase_keys:
            keys -= showcase_keys
            paths.extend(rna_paths(dataset_root, K422_K9_DATASET))
            paths.extend(fragment_path(dataset_root, K422_K9_DATASET, mark) for mark in K422_K9_DT_MARKS)
    for key in INPUTS:
        if key not in keys:
            continue
        if key.endswith("_index") and input_path(key.removesuffix("_index"), replacements).suffix != ".bam":
            continue
        paths.append(input_path(key, replacements))
    paths = [p.resolve() for p in paths]
    missing=[str(p) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input files:\n"+"\n".join(missing))
    if len(set(paths))!=len(paths):
        raise ValueError("INPUTS must list each file only once")
    return paths


def check_outputs(out, expected=EXPECTED_PDFS):
    actual={str(p.relative_to(out)) for p in Path(out).rglob("*.pdf")}
    if actual!=expected:
        raise RuntimeError(f"PDF list differs. Missing: {sorted(expected-actual)}; extra: {sorted(actual-expected)}")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--public-only",action="store_true",help="Run the public benchmark panels only")
    mode.add_argument("--mtr-rna-only",action="store_true",help="Build scRNA_D0.h5ad from the GEO RNA triplets only")
    parser.add_argument("--dataset-root", type=Path,
                        help="Root of the merged in-house RNA and histone dataset folders")
    parser.add_argument("--check-inputs",action="store_true",help="Check required paths and stop")
    args=parser.parse_args()
    logging.basicConfig(level=logging.INFO,format="%(message)s")
    required = None
    if args.mtr_rna_only:
        required = MTR_RNA_INPUTS
    elif args.public_only:
        required = public_input_keys()
    dataset_root = args.dataset_root.expanduser().resolve() if args.dataset_root is not None else None
    replacements = dataset_inputs(dataset_root)
    paths=check_inputs(required, dataset_root, replacements)
    if args.check_inputs:
        print(f"All {len(paths)} required files exist.");return
    if not args.mtr_rna_only:
        executables = ["macs3", "sort", "awk"]
        if not args.public_only and dataset_root is None:
            executables.append("samtools")
        missing=[x for x in executables if not shutil.which(x)]
        if missing: raise RuntimeError("Required executables absent from PATH: "+", ".join(missing))
    out=OUTPUT_ROOT.expanduser().resolve()
    if out.exists(): raise FileExistsError(f"Choose a new OUTPUT_ROOT; it already exists: {out}")
    if any(out==p.parent or p.parent in out.parents for p in paths):
        raise ValueError("OUTPUT_ROOT must be outside input assay directories")
    out.mkdir(parents=True)
    if args.mtr_rna_only:
        build_mtr_rna(out / "scRNA_D0.h5ad")
        return
    SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="benchmark-figures-", dir=SCRATCH_ROOT) as temp:
        work=Path(temp);memo={}
        qc,frip=compute_public(work,memo,out,replacements);plot_public(qc,frip,out)
        if not args.public_only:
            metrics,hist=k422_k9_dt_compute(work,dataset_root);k422_k9_dt_plot(metrics,hist,out)
            qc,metrics,hist=compute_inhouse(work,memo,replacements);plot_inhouse(qc,metrics,hist,out)
    expected = EXPECTED_PDFS
    if args.public_only:
        expected = {p for p in expected if p.split("/", 1)[0] in (GROUP_PUBLIC, GROUP_FRIP)}
    check_outputs(out, expected)
    print(f"Wrote {len(expected)} PDFs under {out}")


if __name__=="__main__":
    sys.exit(main())
