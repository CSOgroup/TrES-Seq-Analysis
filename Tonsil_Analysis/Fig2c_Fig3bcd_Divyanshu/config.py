"""Paths and constants for the manuscript notebooks."""

from __future__ import annotations

from pathlib import Path

PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parent.parent
CODES = REPO / "codes"
PAPER_ASSETS = PACKAGE / "paper_assets"
PAPER_RESULTS = PACKAGE / "paper_results"

# Preprocessed modality objects (counts + embeddings).
PREPROCESSED_ROOT = Path(
    "/mnt/dataNas/ahrmad/tonsil_qc_workflow/tonsil_final_reviewed_protocol/full"
)
PREPROCESSED_OBJECTS = PREPROCESSED_ROOT / "objects"
PREPROCESSED_PER_SAMPLE = PREPROCESSED_ROOT / "per_sample"

MERGED_H5AD = {
    "RNA": PREPROCESSED_OBJECTS / "RNA_merged.h5ad",
    "H3K27ac": PREPROCESSED_OBJECTS / "H3K27ac_merged.h5ad",
    "H3K27me3": PREPROCESSED_OBJECTS / "H3K27me3_merged.h5ad",
}


def tonsil1_rna_h5ad() -> Path:
    return PREPROCESSED_PER_SAMPLE / "Tonsil1" / "Tonsil1_RNA_final_counts.h5ad"


def tonsil1_dna_h5ad(modality: str) -> Path:
    return PREPROCESSED_PER_SAMPLE / "Tonsil1" / f"Tonsil1_{modality}_final_peaks.h5ad"


BCELLS = REPO / "processed" / "four_tonsils_bcellsonly" / "bcells.h5mu"
STEPHANIE_MERS = REPO / "assets" / "Sungalee_2021_MERs.liftedTo.hg38.bed"
PROMOTERS_BED = PAPER_ASSETS / "promoters_tss_window.bed"
CHROM_SIZES = PAPER_ASSETS / "hg38_chrom_sizes.bed"

MODALITIES = ("RNA", "H3K27ac", "H3K27me3")
EMBEDDING_KEY = {
    "RNA": "X_pca",
    "H3K27ac": "X_spectral_harmony",
    "H3K27me3": "X_spectral_harmony",
}
NEIGHBORS_KEY = "neighbors"
CONNECTIVITIES_KEY = "connectivities"

# Column name used in these notebooks. Preprocessed objects store the same
# labels under ``lymphoid_final``; loaders rename into CellType.
CELLTYPE_KEY = "CellType"
CELLTYPE_SOURCE_KEY = "lymphoid_final"
SAMPLE_KEY = "sample"

# Leiden grid used for the cell-type resolution-sweep panels.
LEIDEN_RES_SWEEP = tuple(i / 4 for i in range(1, 13))  # 0.25 … 3.00
RANDOM_STATE = 42

TONSIL1_N_CELLS = 2831
FOUR_TONSILS_N_CELLS = 4049

# Modality colors for resolution-sweep curves (paper legend).
MOD_COLORS = {
    "RNA": "#015D32",
    "H3K27ac": "#272B87",
    "H3K27me3": "#E91D21",
}

# Frozen Leiden-vs-cell-type tables used by notebooks 01 / 02 (exact paper curves).
TONSIL1_LEIDEN_VS_CELLTYPE_CSV = PAPER_RESULTS / "tonsil1_celltype_leiden_resolution.csv"
FOUR_TONSILS_LEIDEN_VS_CELLTYPE_CSV = PAPER_RESULTS / "four_tonsils_celltype_leiden_resolution.csv"

# Frozen stage-UMAP coordinates for notebook 03 (exact paper panels).
# Views: all_tiles, promoters, distal_enhancers. Stages stored as preGC/GC/postGC.
H3K27AC_STAGE_UMAP_NPZ = PAPER_RESULTS / "H3K27ac_stage_umap_coordinates.npz"
H3K27ME3_STAGE_UMAP_NPZ = PAPER_RESULTS / "H3K27me3_stage_umap_coordinates.npz"
STAGE_UMAP_NPZ = {
    "H3K27ac": H3K27AC_STAGE_UMAP_NPZ,
    "H3K27me3": H3K27ME3_STAGE_UMAP_NPZ,
}
STAGE_UMAP_VIEWS = ("all tiles", "promoters", "distal enhancers")
STAGE_UMAP_NPZ_KEYS = {
    "all tiles": "all_tiles",
    "promoters": "promoters",
    "distal enhancers": "distal_enhancers",
}
STAGE_UMAP_EXPECTED_N_BINS = {
    "H3K27ac": {
        "all tiles": 1_618_834,
        "promoters": 31_143,
        "distal enhancers": 77_350,
    },
    "H3K27me3": {
        "all tiles": 260_058,
        "promoters": 15_645,
        "distal enhancers": 8_030,
    },
}
BCELL_STAGE_N = 2635

# Frozen differential-test tables for notebook 05 volcanoes (exact paper panels).
# Stored stage names in filenames / pair column: preGC / GC / postGC.
DIFF_TESTED_DIR = PAPER_RESULTS
DIFF_TESTED_FILES = {
    "promoters": {
        "preGC_vs_GC": DIFF_TESTED_DIR / "promoter_tested_preGC_vs_GC.tsv",
        "GC_vs_postGC": DIFF_TESTED_DIR / "promoter_tested_GC_vs_postGC.tsv",
        "preGC_vs_postGC": DIFF_TESTED_DIR / "promoter_tested_preGC_vs_postGC.tsv",
    },
    "enhancers": {
        "preGC_vs_GC": DIFF_TESTED_DIR / "enhancer_tested_preGC_vs_GC.tsv",
        "GC_vs_postGC": DIFF_TESTED_DIR / "enhancer_tested_GC_vs_postGC.tsv",
        "preGC_vs_postGC": DIFF_TESTED_DIR / "enhancer_tested_preGC_vs_postGC.tsv",
    },
}

# Promoter annotation BED stores TSS −2 kb / +500 bp; analysis window is rebuilt.
STORED_PROMOTER_UPSTREAM = 2000
STORED_PROMOTER_DOWNSTREAM = 500
PROMOTER_UPSTREAM = 1000
PROMOTER_DOWNSTREAM = 250

# Stage nomenclature used in these notebooks. Stored objects may still use
# preGC / postGC; loaders remap into Naive / Memory.
STAGE_ORDER = ["Naive", "GC", "Memory"]
STAGE_PALETTE = {"Naive": "#39437E", "GC": "#3CB24C", "Memory": "#28A9DD"}
STAGE_SOURCE_MAP = {"preGC": "Naive", "GC": "GC", "postGC": "Memory"}
STAGE_OF = {
    "B_naive": "Naive",
    "B_Cycling": "GC",
    "B_GC_DZ": "GC",
    "B_GC_LZ": "GC",
    "B_mem": "Memory",
}
B_ORDER = ["B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive"]
B_DISPLAY = {
    "B_Cycling": "Cycling",
    "B_GC_DZ": "GC dark zone",
    "B_GC_LZ": "GC light zone",
    "B_mem": "Memory",
    "B_naive": "Naive B-cells",
}
B_PALETTE = {
    "B_Cycling": "#57BAC0",
    "B_GC_DZ": "#A8CF91",
    "B_GC_LZ": "#52A053",
    "B_mem": "#488B9D",
    "B_naive": "#3B76B0",
}
SAMPLE_ORDER = ["Tonsil1", "Tonsil2", "Tonsil3", "Tonsil4"]
SAMPLE_PALETTE = {
    "Tonsil1": "#4C78A8",
    "Tonsil2": "#F58518",
    "Tonsil3": "#54A24B",
    "Tonsil4": "#E45756",
}

VIEW_COLORS = {
    "all tiles": "#333333",
    "promoters": "#4C78A8",  # blue
    "enhancers": "#54A24B",  # green (light band in plots)
    "distal enhancers": "#72B7B2",
    "random tiles": "#333333",
}

# B-cell chromatin notebooks
BCELL_RESOLUTIONS = [round(i / 10, 1) for i in range(1, 11)]
N_COMPS = 30
N_NEIGHBORS = 15
MIN_DIST = 0.3
BCELL_RANDOM_STATE = 0
N_MATCHED_REPEATS = 20
