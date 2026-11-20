# TrES-seq tonsil analysis

This workflow **starts from TrESFlow outputs** for four human tonsil samples. If you want to start from scratch, run TrESFlow separately with the study's GRCh38 reference and library assignments. This workflow produces RNA and histone-mark profiles, study figures, and differential-analysis tables.

## Setup and inputs

On Linux, create the environment with micromamba, mamba or conda:

```bash
micromamba create -f environment.yml
```

Set **`tresflow_root`** in `study.json` to the folder containing `Tonsil1/` through `Tonsil4/`. The default is `inputs/tresflow`; you can point it to your existing output folder. Set `reference_dir`, `output_root`, or `threads` in the same file. Relative roots are resolved from this repository.

To start from dataset-level RNA counts and histone fragments, set `dataset_root` in `study.json` to the folder containing those files.

| Input | Path |
|---|---|
| RNA | Under `tresflow_root`: `<sample>/rna_align/<library>.Solo.outGeneFull/filtered/` with `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` |
| Histone marks | Under `tresflow_root`: `<sample>/dna_align/<library>_<mark>_NoDup.bam` |
| Gene-linked enhancers | Under `reference_dir`: `Sungalee_2021_MERs.hg38.proximity_100kb.bed` (BED4; column 4 is the assigned gene) |
| Study curation | `cell_curation.tsv.gz`, included here |

The count-bearing lymphoid reference must be supplied separately. Other missing public references should be downloaded into `reference_dir` as needed. [GEO lymphoid H5AD](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269695/suppl/GSE269695_sc_reference.h5ad), [UCSC chromosome sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes), [Boyle Lab blacklist](https://github.com/Boyle-Lab/Blacklist/blob/d0527cdea356bf57d30535b48783804ca8d21031/lists/hg38-blacklist.v2.bed.gz), [UCSC knownGene GTF](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz) for TSS enrichment, and [GENCODE v46 comprehensive ALL GTF](https://www.gencodegenes.org/human/release_46.html). Existing reference files take precedence.

## Run

Keep both Python scripts and `differential.R` in the same directory. Run processing first, then downstream analysis:

```bash
micromamba run -n tonsil-analysis python tonsil_processing.py
micromamba run -n tonsil-analysis python tonsil_downstream.py
```

The first script reads TrESFlow outputs and writes the finalized profiles, retained fragments, and study figures. The second reads those saved outputs. If they already exist, run the downstream script directly with the same `output_root`.

| Script | Arguments | Outputs |
|---|---|---|
| `tonsil_processing.py` | `process` (optionally `--sample Tonsil1`) | `results/processing/<sample>/`: raw RNA counts, doublet scores, and native histone candidates |
| `tonsil_processing.py` | `label` | `results/final/`: per-sample and merged raw count objects, retained fragments, and peaks. Venn, QC, composition, linked UMAP, and agreement outputs are in `results/publication/`. |
| `tonsil_downstream.py` | `pairwise` | `results/pairwise/<dataset>/<pair>/`: DGE tables and volcanoes, random-subset markers, cell-type and shared gene-signal UMAPs, and comparison-root `gene_group_summaries.tsv` |
| `tonsil_downstream.py` | `dz-lz` | `results/dz_lz/<scope>/{donor_grouped,pseudo_grouped}/`: DGE and volcanoes, plus `all_genes_effect_sizes.pdf` |

`label` uses all four donors' processing outputs. `pairwise` defaults to the study pairs in merged and Tonsil1; to run one pair, add `--dataset Tonsil1 --label-a B_mem --label-b B_naive`. `dz-lz` defaults to merged, Tonsil1 and Tonsil1 >600; use `--scopes Tonsil1` for one scope. Use a separate `output_root` for changed inputs or settings.
