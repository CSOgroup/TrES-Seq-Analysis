# Transcriptional and Epigenetic State by Sequencing (TrES-Seq)

Code repository related to analysing TrES-Seq data


## Overview

The analysis is organized into two stages.

First, FASTQ reads are processed into matrices and fragment files using the automated [TrESFlow pipeline](https://github.com/CSOgroup/TrESFlow). These outputs are also available from GEO under **GSE324511**.

Second, the matrices and fragment files are used as input for the Python scripts in **ProcessingScripts**, which perform QC and clustering. The resulting processed data are then used by **FiguresScripts** to reproduce the figures in the publication.

## Step-by-step

### 1. Data availability
Unprocessed (Sample-split FASTQs) and processed files (scRNA-seq matrices | [STARSolo](https://github.com/alexdobin/STAR) & H3K27ac/H3K27me3/H3K9me3 fragments files | [SnapATAC2](https://github.com/scverse/SnapATAC2)) are available on GEO under the Accession Number: **GSE324511**.

### 2. Environment
All analyses were made using the environment in **env**. To recreate that environment you can use:
```
conda env create -f env/environment.yaml
```
### 3. Modify the config
To run a script, you will need to open it to modify the path to the input files of that script.

### 4. Raw non-demultiplexed reads processing
If you'd like to run the analyses from scratch using the non-demultiplexed FASTQs, use the [TrESFlow pipeline](https://github.com/CSOgroup/TrESFlow) and ask us for the raw non-demultiplexed FASTQ files at ***ahrmad . annan @ unil . ch***
