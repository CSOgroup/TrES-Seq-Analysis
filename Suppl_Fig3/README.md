# For Figures in Supplementary Fig. 3

## Set up and run

Use an analysis environment with Scanpy, AnnData, SnapATAC2, Matplotlib, matplotlib-venn, UpSetPlot, samtools, MACS3, BEDTools, deepTools, and GNU Parallel.

Set `DATA_ROOT` to your input folder. The default is `inputs/` beside `run.sh`.

```bash
cd Suppl_Fig3
bash run.sh --check-inputs
bash run.sh
```

The workflow runs stages 0–4 in order and writes to `results/`.

```bash
DATA_ROOT=/path/to/downloaded-data bash run.sh
```

## Required files

Place these files at the paths shown relative to `inputs/`.

### Histone BAMs

```text
Sc_CMD_A20_H3K27ac_NoDup.bam
Sc_CMD_A20_H3K27me3_NoDup.bam
Sc_CMD_GM12878_H3K27ac_NoDup.bam
Sc_CMD_GM12878_H3K27me3_NoDup.bam
Sc_CMD_JJN3_H3K27ac_NoDup.bam
Sc_CMD_JJN3_H3K27me3_NoDup.bam
Sc_CMD_Karpas422_H3K27ac_NoDup.bam
Sc_CMD_Karpas422_H3K27me3_NoDup.bam
Sc_CMD_WSU_H3K27ac_NoDup.bam
Sc_CMD_WSU_H3K27me3_NoDup.bam
```

Use the coordinate-sorted, deduplicated BAMs with their original `CB` tags. 

### RNA counts

Provide `CellLineMix_<line>_RNA_matrix.mtx.gz`, `CellLineMix_<line>_RNA_features.tsv.gz`, and `CellLineMix_<line>_RNA_barcodes.tsv.gz` for each cell line. STARsolo files can also be supplied under `<line>.Solo.outGeneFull/filtered/`.

### Reference files

| File | Used for |
|---|---|
| `regev_lab_cell_cycle_genes.txt` | RNA cell-cycle scoring | `https://github.com/dpcook/emt_dynamics/blob/master/regev_lab_cell_cycle_genes.txt`
| `hg38-blacklist.bed.gz` | Filtering Human DNA feature selection |`https://github.com/Boyle-Lab/Blacklist/tree/master/lists`
| `hg38-blacklist.bed` | Filtering human peak calls |`https://github.com/Boyle-Lab/Blacklist/tree/master/lists`
| `mm10_grcm39LiftOver-blacklist.v2.bed` | Filtering A20 peak calls | `https://github.com/Boyle-Lab/Blacklist/tree/master/lists`
| `human_epdnew_mjdVA.bed` | TSS coordinates for the marker genes | `https://epd.expasy.org/epd/EPDnew_select.php`


## Run stages individually

```bash
python SupplFig3_0_ProcessingCellLinesMix.py
python SupplFig3_1_HeatmapMarkerGenes.py --top-n 50 --min-logfc 1
bash SupplFig3_2_SubsetBAMs.sh
bash SupplFig3_3_CallPeaks.sh
bash SupplFig3_4_DeepToolsLaunch.sh
```