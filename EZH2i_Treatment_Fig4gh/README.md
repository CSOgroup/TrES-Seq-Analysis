# DESeq2 DGE after 96 h of EZH2 inhibition

The script compares DMSO control cells with cells treated for 96 h. It makes three random pseudobulk groups per condition, runs DESeq2, and produces a results CSV and volcano plots.

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

From the repository root, pass the path to the downloaded H5AD.

```bash
python EZH2i_Treatment_Fig4gh/deseq2_pseudobulk_RNA_4DaysTr_vs_DMSO.py \
  --input /path/to/Karpas422_EZH2i_Study_RNA_processed.h5ad
```

The script creates these three files under `results/DESeq2_RNA_4DaysTr_vs_DMSO/`.

```text
01__4_Days_Tr_vs_DMSO__DESeq2.csv
01__4_Days_Tr_vs_DMSO__DESeq2__volcano.pdf
01__4_Days_Tr_vs_DMSO__DESeq2__volcano.png
```

Use `--outdir /path/to/selected` to change the destination.