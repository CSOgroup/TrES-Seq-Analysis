# treseq_paper

Notebooks that reproduce the manuscript analyses.

## Notebooks

| Notebook | Main-paper figure | Analysis |
|---|---|---|
| `01_tonsil1_leiden_vs_celltype.ipynb` | Figure 2 | Tonsil1 (n=2,831) Leiden resolution sweep vs `CellType` (ARI, NMI, conditional entropy, homogeneity); curves from frozen table |
| `02_four_tonsils_leiden_vs_celltype.ipynb` | Figure 3 | Tonsil1–4 (n=4,049) same sweep; curves from frozen table |
| `03_bcells_chromatin_views.ipynb` | Figure 3 | Non-plasma B cells: H3K27ac / H3K27me3 stage UMAPs (Naive / GC / Memory) on all tiles, promoters, distal enhancers; from frozen coords |
| `04_matched_random_tiles.ipynb` | Figure 3 | H3K27ac size-matched promoters vs random enhancers vs random genomic tiles vs stage (Naive / GC / Memory; 20 repeats) |
| `05_differential_volcanoes.ipynb` | Figure 3 | Stage-pair H3K27ac volcanoes on promoters and Stephanie MERs; from frozen tested tables |

Figure and table disk writes are commented out; plots and tables stay in the notebook.

## Inputs

- Preprocessed modality objects (RNA / H3K27ac / H3K27me3 embeddings + counts); paths in `config.py`
- `processed/four_tonsils_bcellsonly/bcells.h5mu` (B-cell notebooks; see below)
- `assets/Sungalee_2021_MERs.liftedTo.hg38.bed`
- Bundled `treseq_paper/paper_assets/promoters_tss_window.bed` (TSS annotation stored as −2 kb / +500 bp; analysis windows rebuilt to −1 kb / +250 bp)
- Bundled `treseq_paper/paper_assets/hg38_chrom_sizes.bed`
- Frozen manuscript tables / coordinates in `treseq_paper/paper_results/`:
  - `tonsil1_celltype_leiden_resolution.csv`, `four_tonsils_celltype_leiden_resolution.csv` (notebooks 01 / 02)
  - `H3K27ac_stage_umap_coordinates.npz`, `H3K27me3_stage_umap_coordinates.npz` (notebook 03; all tiles / promoters / distal enhancers)
  - `promoter_tested_*.tsv`, `enhancer_tested_*.tsv` (notebook 05; preGC/GC/postGC pairs)

Cell-type labels are exposed as `CellType` in these notebooks (renamed from the stored annotation column on load).

Notebooks **01** and **02** plot frozen Leiden tables; **03** plots frozen UMAP coordinates; **05** plots frozen differential-test TSVs so panels match the manuscript figures. A commented regenerate path remains in each notebook; a full recompute can shift values slightly because neighbor graphs, Leiden, Harmony, UMAP, and differential testing are intrinsically stochastic even with a fixed seed.

## Where the calculation code lives

Frozen files in `paper_results/` are only the **numbers and coordinates** used for the paper plots. The definitions of every metric and embedding step ship in this folder:

| What | Where to read it |
|---|---|
| ARI, NMI, AMI, homogeneity / completeness, Shannon entropy, conditional entropy in bits (`H(cell type \| clusters)`, `H(clusters \| cell type)`) | `metrics.py` |
| Trimodal overlap, neighbor recompute, Leiden (`flavor=igraph`, `n_iterations=2`, `random_state=42`), resolution sweep that calls `metrics.vs_celltype`, resolution-sweep plotting | `load_resolution.py` |
| Paths, resolutions, palettes, stage maps | `config.py` |
| Notebooks **01** / **02** | Active cells load and plot the frozen tables. Optional commented cells call `load_tonsil1` / `load_four_tonsils` + `leiden_vs_celltype` (implementations in the `.py` files above). |
| Notebook **03** | Active cells load and plot frozen UMAP coordinates. Optional commented cell contains the full recompute: promoter / distal-enhancer masks, SnapATAC2 cosine spectral, Harmony, neighbors, UMAP. |
| Notebook **04** | All size-matched draws, spectral / Harmony / UMAP, Leiden, and ARI / NMI / conditional entropy vs stage are **live** in the notebook. |
| Notebook **05** | Active cells load frozen `*_tested_*.tsv` and draw the volcano grid. Optional commented cell shows interval aggregation + `snapatac2.tl.diff_test` (Naive/GC/Memory). |

Nothing required for these figure panels is deferred to another analysis tree.

## How `bcells.h5mu` is built

It is a cell-type subset of the preprocessed **merged** four-tonsil objects, then re-embedded. Plasma cells are dropped; the five retained subtypes are `B_Cycling`, `B_GC_DZ`, `B_GC_LZ`, `B_mem`, and `B_naive`. Two non-overlapping intersections are used: RNA∩H3K27ac∩H3K27me3 (2,635 cells) and RNA∩H3K27me3∩H3K9me3 (1,652 cells). RNA embeddings use the union; H3K27ac uses the first set only. Stages are `Naive` (`B_naive`), `GC` (`B_Cycling` / `B_GC_DZ` / `B_GC_LZ`), and `Memory` (`B_mem`).

Selection sketch (paths from `config.MERGED_H5AD` plus the H3K9me3 merged object):

```python
import anndata as ad
import numpy as np
import pandas as pd

B_ORDER = ["B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive"]
STAGE_OF = {
    "B_naive": "Naive",
    "B_Cycling": "GC",
    "B_GC_DZ": "GC",
    "B_GC_LZ": "GC",
    "B_mem": "Memory",
}
# Preprocessed objects store labels under this column; notebooks rename to CellType.
SRC = "lymphoid_final"

def read_obs(path):
    adata = ad.read_h5ad(path, backed="r")
    try:
        return pd.DataFrame(
            {
                "CellType": adata.obs[SRC].astype(str).to_numpy(),
                "sample": adata.obs["sample"].astype(str).to_numpy(),
            },
            index=pd.Index(adata.obs_names.astype(str), name="cell_id"),
        )
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()

def b_cell_ids(obs_by_assay, modalities):
    shared = obs_by_assay["RNA"].index
    for assay in modalities:
        if assay != "RNA":
            shared = shared.intersection(obs_by_assay[assay].index)
    rna_labels = obs_by_assay["RNA"].loc[list(shared), "CellType"]
    return [cid for cid in shared if rna_labels.loc[cid] in B_ORDER]

# obs_by_assay = {assay: read_obs(path) for assay, path in PATHS.items()}
# ac_ids = b_cell_ids(obs_by_assay, ["RNA", "H3K27ac", "H3K27me3"])   # 2635
# k9_ids = b_cell_ids(obs_by_assay, ["RNA", "H3K27me3", "H3K9me3"])   # 1652
# Then subset count matrices, drop empty features, embed (RNA: HVG+PCA;
# chromatin: SnapATAC2 spectral + Harmony on sample), and write bcells.h5mu.
```

## Library versions

Recorded from the `cnt` environment used for these analyses:

```
scanpy==1.10.2
snapatac2==2.7.1
muon==0.1.6
harmonypy==2.0.0
anndata==0.11.4
numpy==1.26.4
pandas==2.3.3
scipy==1.15.3
scikit-learn==1.7.2
matplotlib==3.10.6
seaborn==0.13.2
pybedtools==0.10.0
hdf5plugin==7.1.0
igraph==0.11.9
leidenalg==0.10.2
```
