# Benchmarking

Set the input paths at the top of each script for your installation. The `SOURCES` entries record where the inputs came from. The `--check-inputs` option checks paths for the three independent workflows below.

| Workflow | Command | Output directory |
| --- | --- | --- |
| Benchmark figures | `python Make_Benchmark_Figures.py` | `benchmark_figures_regenerated/` |
| Cross-mark profiles and peak BEDs | `python Make_CrossCont_Figures.py` | `crosscont_figures/` |
| Panel k correlations | `python Make_Correlation_Heatmaps.py` | `panel_k_figures/` |

Use `--public-only` for the public benchmark.

Use `--dataset-root /path/to/data` for the benchmark and cross-mark workflows, and `--input-dir /path/to/data` for correlation heatmaps.

## Single-cell tracks

Run step 0 on three coordinate-sorted, indexed hg38 BAMs containing the cells eligible for plotting. It writes the locus BAMs to `inputs/`, where step 1 reads them. Start at step 1 if those locus BAMs are already available.

```sh
python 0_subset_single_cell_bams.py prepare \
  --h3k27ac /path/to/H3K27ac.bam \
  --h3k27me3 /path/to/H3K27me3.bam \
  --rna /path/to/RNA.bam
Rscript 1_make_single_cell_tracks.R
```

The R script reads cell IDs from the BAM tags and writes plots to `single_cell_stacks/`.