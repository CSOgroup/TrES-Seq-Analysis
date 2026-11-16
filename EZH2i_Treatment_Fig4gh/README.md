# DESeq2 DGE from RNA_DMSO_4DaysTr_processed.h5ad

The `reproduce.py` script subsets expression to **DMSO** and **4 Days Tr.**,
makes three random pseudobulk groups per condition, runs DESeq2,
and produces the results CSV and volcano plots.

## Requirements

python==3.11
numpy==1.26.4
pandas==2.1.1
matplotlib==3.10.8
adjustText==1.4.0
anndata==0.12.10
scipy==1.17.1

 + R 4.5.3, DESeq2 1.50.2.

## Run

From the directory containing `RNA_DMSO_4DaysTr_processed.h5ad`:

```bash
python deseq2_pseudobulk_RNA_4DaysTr_vs_DMSO.py --input RNA_DMSO_4DaysTr_processed.h5ad
```

The script creates these three files under `results/DESeq2_RNA_4DaysTr_vs_DMSO/`.

```text
01__4_Days_Tr_vs_DMSO__DESeq2.csv
01__4_Days_Tr_vs_DMSO__DESeq2__volcano.pdf
01__4_Days_Tr_vs_DMSO__DESeq2__volcano.png
```

Use `--outdir /path/to/selected` to change the destination.