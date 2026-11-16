# TrES-seq tonsil analysis

This workflow **starts from TrESFlow outputs** for four human tonsil samples. Run TrESFlow separately with the study's GRCh38 reference and library assignments. The workflow produces RNA and histone-mark profiles, study figures, and differential-analysis tables.

## Setup and inputs

On Linux, create the environment with micromamba:

```bash
micromamba create -f environment.yml
```

Set **`tresflow_root`** in `study.json` to the folder containing `Tonsil1/` through `Tonsil4/`. The default is `inputs/tresflow`; you can point it to your existing output folder. Set `reference_dir`, `output_root`, or `threads` in the same file. Relative roots are resolved from this repository.

| Input | Path |
|---|---|
| RNA | Under `tresflow_root`: `<sample>/rna_align/<library>.Solo.outGeneFull/filtered/` with `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` |
| Histone marks | Under `tresflow_root`: `<sample>/dna_align/<library>_<mark>_MarkedDup.bam` |
| Gene-linked enhancers | Under `reference_dir`: `Sungalee_2021_MERs.hg38.proximity_100kb.bed` (BED4; column 4 is the assigned gene) |
| Study curation | `cell_curation.tsv.gz`, included here |

Missing public references are downloaded into `reference_dir` as needed: the [GEO lymphoid H5AD](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269695/suppl/GSE269695_sc_reference.h5ad), [UCSC chromosome sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes), the [Boyle Lab blacklist](https://github.com/Boyle-Lab/Blacklist/blob/d0527cdea356bf57d30535b48783804ca8d21031/lists/hg38-blacklist.v2.bed.gz), the [UCSC knownGene GTF](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) for TSS enrichment, and the [GENCODE v46 comprehensive ALL GTF](https://www.gencodegenes.org/human/release_46.html) for effect sizes. Existing reference files take precedence. The blacklist is decompressed to BED. SnapATAC2 obtains GENCODE v41 for stages 03–04 automatically.

For the original TSS enrichment measurements, keep the exact TrESFlow `genes.gtf.gz` in `reference_dir`. The UCSC fallback is used for TSS enrichment and does not change the default cell thresholds.

**Lymphoid reference.** The default is `GSE269695_sc_reference.h5ad` from GEO. Stage 02 uses `obs["new_celltype"]`, removes `Tumor`, and then checks the retained matrix for raw UMI counts in `layers["counts"]`, `.X`, or `.raw.X`. If the download contains only processed expression, use a verified count-bearing H5AD: set `references.lymphoid.file` to its filename and remove that entry's `url`. The workflow reports an error rather than training on normalized values.

## Run

Keep both Python scripts and `differential.R` in the same directory. Run processing first, then downstream analysis:

```bash
micromamba run -n tonsil-analysis python tonsil_processing.py
micromamba run -n tonsil-analysis python tonsil_downstream.py
```

The first script reads TrESFlow outputs and writes the finalized profiles, retained fragments, and study figures. The second reads those saved outputs. If they already exist, run the downstream script directly with the same `output_root`.

To download and validate the GEO reference before starting:

```bash
micromamba run -n tonsil-analysis python tonsil_processing.py references --only lymphoid
```

| Script | Arguments | Outputs |
|---|---|---|
| `tonsil_processing.py` | `process` (optionally `--sample Tonsil1`) | `results/processing/<sample>/`: raw RNA counts, doublet scores, and native histone candidates |
| `tonsil_processing.py` | `label` | `results/final/`: per-sample and merged raw count objects, retained fragments, and peaks. Venn, QC, composition, linked UMAP, and agreement outputs are in `results/publication/`. |
| `tonsil_downstream.py` | `pairwise` | `results/pairwise/<dataset>/<pair>/`: DGE tables and volcanoes, random-subset markers, cell-type and shared gene-signal UMAPs |
| `tonsil_downstream.py` | `dz-lz` | `results/dz_lz/<scope>/{donor_grouped,pseudo_grouped}/`: DGE and volcanoes, plus `all_genes_effect_sizes.pdf` |

`label` uses all four donors' processing outputs. `pairwise` defaults to the study pairs in merged and Tonsil1; to run one pair, add `--dataset Tonsil1 --label-a B_mem --label-b B_naive`. `dz-lz` defaults to merged, Tonsil1 and Tonsil1 >600; use `--scopes Tonsil1` for one scope. Use a separate `output_root` for changed inputs or settings.

## Analysis details

Processing applies per-sample technical thresholds and exact-ID curation. Recorded DoubletDetection 4.3.0.post1 scores are reused when the eligible RNA IDs match exactly and the scores are usable. Otherwise the workflow fits DoubletDetection on the raw RNA counts. Stage 02 removes cells with a raw score greater than 70 before final RNA clustering. Signal-score cutoffs in `study.json` document the scores used to curate histone exclusions; their cell IDs are in the curation table.

Labels are transferred with scVI/scANVI, merged into study cell types and smoothed on the within-donor RNA graph. Final RNA labels are joined to histone profiles by exact physical cell ID. RNA uses PCA; histone profiles use spectral embeddings of common tiles (H3K27ac 1 kb; H3K27me3/H3K9me3 10 kb). Harmony corrects merged H3K27ac and H3K27me3 for single versus dual tagmentation. The >600 filters are applied independently by modality. Agreement uses the exact RNA/H3K27ac/H3K27me3 intersection; ARI and NMI use its labelled cells.

Pairwise RNA contrasts are A / B. DZ–LZ DESeq2 contrasts are DZ / LZ, with histone counts from gene bodies ±1 kb, alone or unioned with linked enhancers. Merged donor-paired models use biological donors as replicates. Results based on random cell subdivisions are exploratory. Effect-size panels use the corresponding RNA DGE for ranking and estimate LZ / DZ effects on matched RNA/H3K27ac/H3K27me3 cells. Histone effects use gene-body ±2 kb counts plus linked enhancer counts. Unsupported estimates are marked NA.
