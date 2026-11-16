# Transcriptional and Epigenetic State by Sequencing (TrES-Seq)

Analysis code for the TrES-Seq paper. Each folder's README lists the inputs, software and commands for that analysis.

| Folder | Analysis |
| --- | --- |
| [Benchmarking](Benchmarking/) | Benchmark comparisons, cross-mark profiles, correlations and single-cell tracks |
| [Tonsil_Analysis](Tonsil_Analysis/) | Processing and downstream analysis of four human tonsil samples |
| [EZH2i_Treatment_Fig4gh](EZH2i_Treatment_Fig4gh/) | Pseudobulk RNA differential expression after four days of EZH2 inhibition versus DMSO (Fig. 4g,h) |
| [Suppl_Fig3](Suppl_Fig3/) | Cell-line mixture analysis and histone-mark profiles (Supplementary Fig. 3) |

## Data and preprocessing

The GEO accession for the study's sample-split FASTQs, RNA count matrices and histone-mark fragment files is [GSE324511](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324511).

Use [TrESFlow](https://github.com/CSOgroup/TrESFlow) to process FASTQs. For the original, non-demultiplexed FASTQs, contact **ahrmad . annan @ unil . ch**.
