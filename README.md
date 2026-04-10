# Joint Single Cell Profile of RNA and Histone Marks To Reveal Epigenetic and Transcriptional Dynamics (TrES-Seq)

Code repository related to analysing TrES-Seq data


## Overview

The analysis is made of 2 major parts, preprocessing (**PeprocessingScripts**) and the launchers (**ProcessingScripts**). The preprocessing scripts are needed for the launchers to demultiplex the raw reads, align them and process the result into usable matrices and fragment files. These files can then be used as input to reproduce the figures in the publication using the scripts in **FiguresScripts**. Alternatively processed files can be downloaded to directly start at the last step.

A Nextflow pipeline is currently being developed to automate the workflow through to the production of processed matrices.

## Step-by-step

### 1. Data availability
Raw files are available on GEO under the Accession Number: GSE324511.

Processed files - scRNA-seq matrices in 10x mtx format (RNA | [STARSolo](https://github.com/alexdobin/STAR)), fragments files (H3K27ac/H3K27me3 | [SnapATAC2](https://github.com/scverse/SnapATAC2)) are available as supplementary files on the same GEO repository

### 2. Environment
All analyses were made using the environment in **env**. To recreate that environment you can use:
```
conda env create -f env/environment.yaml
```
To preprocess the raw fastq files you will need to install [Codon](https://github.com/exaloop/codon) and [Seq](https://github.com/exaloop/seq), in addition to this environment. To install them, follow the instructions listed in the links.

### 3. Modify the config
To run a script, you will need to open it to modify the path to the input files of that script and the required accompanying scripts (only for preprocessing).

### 4. Raw reads processing
Before starting the preprocessing, first use the scripts in the **PeprocessingScripts/GenerateSTARSoloGenome** folder to create the proper STAR genomes that will be used during processing. Then to process a specific dataset, use *Launch.sh* inside the corresponding folder in **ProcessingScripts**. The attached whitelists *WLs* are required to demultiplex the reads. You can then use *PostProcessing.py* to obtain the relevant matrices from the aligned reads create in the previous step.
