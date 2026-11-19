# Internal preparation for data publication on GEO

This script merges library-level STARsolo RNA counts and histone fragments by dataset, and outputs single-cell matrices for each modality with library-prefixed cell barcodes to facilitate data publication on GEO and result reproduction.

```bash
python GEO_Preparation/merge_libraries.py GEO_Preparation/merge_libraries_manifest.json --outdir GEO_datasets
```
