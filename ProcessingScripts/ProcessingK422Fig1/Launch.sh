#!/bin/bash

set -e

wl_cell="/mnt/dataFast/ahrmad/triseq_202409/triseq_202409_cell_whitelist.txt"
wl_mod="/mnt/dataFast/ahrmad/triseq_202409/triseq_202409_dna_modality_whitelist.txt"
wl_sample="/mnt/dataFast/ahrmad/triseq_202409/triseq_202409_sample_whitelist.txt"

outdir="/mnt/dataFast/ahrmad/triseq_202409/processed"
refdir="/mnt/dataFast/ahrmad"

DNA_SName="H2RD1"
RNA_SName1="H2RRO1"
RNA_SName2="H2RRE2"

I2D=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_DNA/avid01543/H2RD1_S1_I2_001.fastq.gz
I1D=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_DNA/avid01543/H2RD1_S1_I1_001.fastq.gz
R1D=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_DNA/avid01543/H2RD1_S1_R1_001.fastq.gz
R2D=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_DNA/avid01543/H2RD1_S1_R2_001.fastq.gz

R2R1=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_RNA/avid01544/H2RRO1_S2_R2_001.fastq.gz
I1R1=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_RNA/avid01544/H2RRO1_S2_I1_001.fastq.gz
R1R1=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_RNA/avid01544/H2RRO1_S2_R1_001.fastq.gz

R2R2=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_RNA/avid01545/H2RRE2_S3_R2_001.fastq.gz
I1R2=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_RNA/avid01545/H2RRE2_S3_I1_001.fastq.gz
R1R2=/mnt/svraw1/Ahrmad/ScH2R_Triseq_20241001/AVT0097/fastq_RNA/avid01545/H2RRE2_S3_R1_001.fastq.gz

# 1 QC
# Fastqc MultiQC pre trimmd
fastqc --memory 2048 --threads 24 --outdir ${outdir}/fastqc_pretrimming /mnt/dataFast/ahrmad/raw/AVT0097/fastq_*NA/avid*/*fastq.gz
multiqc --outdir ${outdir}/fastqc_pretrimming ${outdir}/fastqc_pretrimming/.

# 2 Tagging
# tags DNA: 
# SB: sample barcode: I2  -D BC_LEN=4 -D BC_START=14 REV
# MO: modality barcode I2 -D BC_LEN=8 -D BC_START=18 REV
# L1,L2,L3: ligation barcodes, CB: single-cell barcodes, RG: copy of CB -> I1 hardcoded

# Replace I1/2, R1, R2

codon run -plugin seq -release -D BC_LEN=4 -D BC_START=14 -D HD=1 /home/annan/ScH2R_TriSeq/Tag.codon \
    ${I2D} \
    ${R1D} \
    ${R2D} \
    ${wl_sample} \
    ${DNA_SName} SB ${outdir} first_pass rev

codon run -plugin seq -release -D BC_LEN=8 -D BC_START=18 -D HD=1 /home/annan/ScH2R_TriSeq/Tag.codon \
    ${I2D} \
    ${outdir}/H2RD1_S1_R1_001_SB.fastq \
    ${outdir}/H2RD1_S1_R2_001_SB.fastq \
    ${wl_mod} \
    ${DNA_SName} MO ${outdir} not_first_pass rev

codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/ScH2R_TriSeq/Tag_Lig3.codon \
    ${I1D} \
    ${outdir}/H2RD1_S1_R1_001_SB_MO.fastq \
    ${outdir}/H2RD1_S1_R2_001_SB_MO.fastq \
    ${wl_cell} \
    ${DNA_SName} CB ${outdir}

rm ${outdir}/H2RD1_S1_R1_001_SB.fastq ${outdir}/H2RD1_S1_R2_001_SB.fastq ${outdir}/H2RD1_S1_R1_001_SB_MO.fastq ${outdir}/H2RD1_S1_R2_001_SB_MO.fastq

# tags RNA: 
# SB: sample barcode R2 -D BC_LEN=4 -D BC_START=0 REV
# UM: UMI R2 -D BC_LEN=10 -D BC_START=4
# L1,L2,L3: ligation barcodes, CB: single-cell barcodes, RG: copy of CB -> I1 hardcoded

###RNA1
# SB RNA1

codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/ScH2R_TriSeq/Tag.codon \
    ${R2R1} \
    ${R1R1} \
    ${R2R1} \
    ${wl_sample} \
    ${RNA_SName1} SB ${outdir} first_pass rev

# UM RNA1
codon run -plugin seq -release -D BC_LEN=10 -D BC_START=4 /home/annan/ScH2R_TriSeq/Tag_UMI.codon \
    ${R2R1} \
    ${outdir}/H2RRO1_S2_R1_001_SB.fastq \
    ${outdir}/H2RRO1_S2_R2_001_SB.fastq \
    ${RNA_SName1} UM ${outdir}

# CB RNA1
codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/ScH2R_TriSeq/Tag_Lig3.codon \
    ${I1R1} \
    ${outdir}/H2RRO1_S2_R1_001_SB_UM.fastq \
    ${outdir}/H2RRO1_S2_R2_001_SB_UM.fastq \
    ${wl_cell} \
    ${RNA_SName1} CB ${outdir}

rm ${outdir}/H2RRO1_S2_R1_001_SB.fastq ${outdir}/H2RRO1_S2_R2_001_SB.fastq ${outdir}/H2RRO1_S2_R1_001_SB_UM.fastq ${outdir}/H2RRO1_S2_R2_001_SB_UM.fastq

###RNA2
# SB RNA2

codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/ScH2R_TriSeq/Tag.codon \
    ${R2R2} \
    ${R1R2} \
    ${R2R2} \
    ${wl_sample} \
    ${RNA_SName2} SB ${outdir} first_pass rev

# UM RNA2
codon run -plugin seq -release -D BC_LEN=10 -D BC_START=4 /home/annan/ScH2R_TriSeq/Tag_UMI.codon \
    ${R2R2} \
    ${outdir}/H2RRE2_S3_R1_001_SB.fastq \
    ${outdir}/H2RRE2_S3_R2_001_SB.fastq \
    ${RNA_SName2} UM ${outdir}

# CB RNA2
codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/ScH2R_TriSeq/Tag_Lig3.codon \
    ${I1R2} \
    ${outdir}/H2RRE2_S3_R1_001_SB_UM.fastq \
    ${outdir}/H2RRE2_S3_R2_001_SB_UM.fastq \
    ${wl_cell} \
    ${RNA_SName2} CB ${outdir}

rm ${outdir}/H2RRE2_S3_R1_001_SB.fastq ${outdir}/H2RRE2_S3_R2_001_SB.fastq ${outdir}/H2RRE2_S3_R1_001_SB_UM.fastq ${outdir}/H2RRE2_S3_R2_001_SB_UM.fastq

# 3 Sequencing adaptor trimming
# trim-galore automatic detection of seq adapters
outdir_postT=/mnt/dataFast/ahrmad/triseq_202409/processed/fastqc_posttrimming

trim_galore \
    --quality 10 \
    --cores 8 \
    --fastqc \
    --fastqc_args "--memory 2048 --threads 4 --outdir ${outdir_postT}" \
    --output_dir ${outdir} \
    --dont_gzip \
    --length 20 \
    --paired \
    ${outdir}/H2RD1_S1_R1_001_SB_MO_CB.fastq ${outdir}/H2RD1_S1_R2_001_SB_MO_CB.fastq

trim_galore \
    --quality 10 \
    --cores 8 \
    --fastqc \
    --fastqc_args "--memory 2048 --threads 4 --outdir ${outdir_postT}" \
    --output_dir ${outdir} \
    --dont_gzip \
    --length 20 \
    --paired \
    ${outdir}/H2RRO1_S2_R1_001_SB_UM_CB.fastq ${outdir}/H2RRO1_S2_R2_001_SB_UM_CB.fastq

trim_galore \
    --quality 10 \
    --cores 8 \
    --fastqc \
    --fastqc_args "--memory 2048 --threads 4 --outdir ${outdir_postT}" \
    --output_dir ${outdir} \
    --dont_gzip \
    --length 20 \
    --paired \
    ${outdir}/H2RRE2_S3_R1_001_SB_UM_CB.fastq ${outdir}/H2RRE2_S3_R2_001_SB_UM_CB.fastq

multiqc --outdir ${outdir_postT} ${outdir_postT}/.

rm H2RRE2_S3_R1_001_SB_UM_CB.fastq H2RRE2_S3_R2_001_SB_UM_CB.fastq H2RRO1_S2_R1_001_SB_UM_CB.fastq H2RRO1_S2_R2_001_SB_UM_CB.fastq
# 4 Read splitting
#codon run -plugin seq -release /home/annan/ScH2R_TriSeq/Split_Reads.codon testsample /mnt/dataFast/ahrmad/triseq_202409/test_dna dna R1.fq R2.fq
#codon run -plugin seq -release /home/annan/ScH2R_TriSeq/Split_Reads.codon testsample /mnt/dataFast/ahrmad/triseq_202409/test_rna rna R1.fq R2.fq

codon run -plugin seq -release /home/annan/ScH2R_TriSeq/Split_Reads.codon \
    ${DNA_SName} \
    ${outdir} \
    dna ${outdir}/H2RD1_S1_R1_001_SB_MO_CB_val_1.fq ${outdir}/H2RD1_S1_R2_001_SB_MO_CB_val_2.fq

codon run -plugin seq -release /home/annan/ScH2R_TriSeq/Split_Reads.codon \
    ${RNA_SName1} \
    ${outdir} \
    rna ${outdir}/H2RRO1_S2_R1_001_SB_UM_CB_val_1.fq ${outdir}/H2RRO1_S2_R2_001_SB_UM_CB_val_2.fq

codon run -plugin seq -release /home/annan/ScH2R_TriSeq/Split_Reads.codon \
    ${RNA_SName2} \
    ${outdir} \
    rna ${outdir}/H2RRE2_S3_R1_001_SB_UM_CB_val_1.fq ${outdir}/H2RRE2_S3_R2_001_SB_UM_CB_val_2.fq

rm H2RD1_S1_R1_001_SB_MO_CB_val_1.fq H2RD1_S1_R2_001_SB_MO_CB_val_2.fq
rm H2RRE2_S3_R1_001_SB_UM_CB_val_1.fq H2RRE2_S3_R2_001_SB_UM_CB_val_2.fq H2RRO1_S2_R1_001_SB_UM_CB_val_1.fq H2RRO1_S2_R2_001_SB_UM_CB_val_2.fq

# Concatenate RNA fastqs
cat H2RRE2_Mouse_R1.fq H2RRO1_Mouse_R1.fq > H2RRNA_Mouse_R1.fq
cat H2RRE2_Mouse_R2.fq H2RRO1_Mouse_R2.fq > H2RRNA_Mouse_R2.fq
cat H2RRE2_Human_R1.fq H2RRO1_Human_R1.fq > H2RRNA_Human_R1.fq
cat H2RRE2_Human_R2.fq H2RRO1_Human_R2.fq > H2RRNA_Human_R2.fq

cat SAM_RG_Header_H2RRE2_Mouse.tsv SAM_RG_Header_H2RRO1_Mouse.tsv | awk -F'\t' 'BEGIN {OFS="\t"} {$3 = "SM:H2RRNA"; print}' | sort --unique --buffer-size=5% > SAM_RG_Header_H2RRNA_Mouse.tsv
cat SAM_RG_Header_H2RRE2_Human.tsv SAM_RG_Header_H2RRO1_Human.tsv | awk -F'\t' 'BEGIN {OFS="\t"} {$3 = "SM:H2RRNA"; print}' | sort --unique --buffer-size=5% > SAM_RG_Header_H2RRNA_Human.tsv

rm H2RRE2_Mouse_R1.fq H2RRO1_Mouse_R1.fq H2RRE2_Mouse_R2.fq H2RRO1_Mouse_R2.fq SAM_RG_Header_H2RRE2_Mouse.tsv SAM_RG_Header_H2RRO1_Mouse.tsv
rm H2RRE2_Human_R1.fq H2RRO1_Human_R1.fq H2RRE2_Human_R2.fq H2RRO1_Human_R2.fq SAM_RG_Header_H2RRE2_Human.tsv SAM_RG_Header_H2RRO1_Human.tsv

# 5 Alignment

#DNA
mm10_eff_size=2652783500
hg38_eff_size=2913022398

mod=H3K27ac
org=Human
/home/annan/ScH2R_TriSeq/AlignDNA.sh ${org} ${mod} ${DNA_SName} \
    ${outdir}/${DNA_SName}_${org}_${mod}_R1.fq \
    ${outdir}/${DNA_SName}_${org}_${mod}_R2.fq \
    ${refdir}/hg38-blacklist.bed \
    ${outdir}/SAM_RG_Header_${DNA_SName}_${org}_${mod}.tsv \
    ${refdir}/hg38_bwamem2_index/hg38.fa ${hg38_eff_size}

mod=H3K27ac
org=Mouse
/home/annan/ScH2R_TriSeq/AlignDNA.sh ${org} ${mod} ${DNA_SName} \
    ${outdir}/${DNA_SName}_${org}_${mod}_R1.fq \
    ${outdir}/${DNA_SName}_${org}_${mod}_R2.fq \
    ${refdir}/mm10-blacklist.bed \
    ${outdir}/SAM_RG_Header_${DNA_SName}_${org}_${mod}.tsv \
    ${refdir}/mm10_bwamem2_index/mm10.fa ${mm10_eff_size}

mod=H3K27me3
org=Mouse
/home/annan/ScH2R_TriSeq/AlignDNA.sh ${org} ${mod} ${DNA_SName} \
    ${outdir}/${DNA_SName}_${org}_${mod}_R1.fq \
    ${outdir}/${DNA_SName}_${org}_${mod}_R2.fq \
    ${refdir}/mm10-blacklist.bed \
    ${outdir}/SAM_RG_Header_${DNA_SName}_${org}_${mod}.tsv \
    ${refdir}/mm10_bwamem2_index/mm10.fa ${mm10_eff_size}

mod=H3K27me3
org=Human
/home/annan/ScH2R_TriSeq/AlignDNA.sh ${org} ${mod} ${DNA_SName} \
    ${outdir}/${DNA_SName}_${org}_${mod}_R1.fq \
    ${outdir}/${DNA_SName}_${org}_${mod}_R2.fq \
    ${refdir}/hg38-blacklist.bed \
    ${outdir}/SAM_RG_Header_${DNA_SName}_${org}_${mod}.tsv \
    ${refdir}/hg38_bwamem2_index/hg38.fa ${hg38_eff_size}

#RNA
org=Human
#/home/annan/ScH2R_TriSeq/AlignRNA.sh ${org} H2RRNA \
#    H2RRNA_${org}_R1.fq H2RRNA_${org}_R2.fq \
#    ${refdir}/STARIndex_hg38

/home/annan/ScH2R_TriSeq/AlignRNA.sh ${org} H2RRNA \
    H2RRNA_${org}_R1.fq H2RRNA_${org}_R2.fq \
    ${refdir}/refdata-gex-GRCh38-2024-A/star

org=Mouse
#/home/annan/ScH2R_TriSeq/AlignRNA.sh ${org} H2RRNA \
#    H2RRNA_${org}_R1.fq H2RRNA_${org}_R2.fq \
#    ${refdir}/STARIndex_mm10

/home/annan/ScH2R_TriSeq/AlignRNA.sh ${org} H2RRNA \
    H2RRNA_${org}_R1.fq H2RRNA_${org}_R2.fq \
    ${refdir}/refdata-gex-GRCm39-2024-A/star

# Duplicate detection
#DNA
mod=H3K27me3
org=Mouse
/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
    -I H2RD1_${org}_${mod}.bam -O H2RD1_${org}_${mod}_MarkedDup.bam \
    -M /mnt/dataFast/ahrmad/triseq_202409/H2RD1_${org}_${mod}.DuplicateMetrics.txt \
    --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000

mod=H3K27ac
/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
    -I H2RD1_${org}_${mod}.bam -O H2RD1_${org}_${mod}_MarkedDup.bam \
    -M /mnt/dataFast/ahrmad/triseq_202409/H2RD1_${org}_${mod}.DuplicateMetrics.txt \
    --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000

mod=H3K27me3
org=Human
/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
    -I H2RD1_${org}_${mod}.bam -O H2RD1_${org}_${mod}_MarkedDup.bam \
    -M /mnt/dataFast/ahrmad/triseq_202409/H2RD1_${org}_${mod}.DuplicateMetrics.txt \
    --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000

mod=H3K27ac
/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
    -I H2RD1_${org}_${mod}.bam -O H2RD1_${org}_${mod}_MarkedDup.bam \
    -M /mnt/dataFast/ahrmad/triseq_202409/H2RD1_${org}_${mod}.DuplicateMetrics.txt \
    --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000

# Spark version does not take BARCODE_TAG parameter... 
#/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicatesSpark \
#    -I H2RD1_${org}_${mod}.bam -O H2RD1_${org}_${mod}_MarkedDup.bam \
#    -M /mnt/dataFast/ahrmad/triseq_202409/H2RD1_${org}_${mod}.DuplicateMetrics.txt \
#    --remove-all-duplicates false --create-output-bam-index --spark-master local[64]

#RNA
#gtfMM=/mnt/dataFast/ahrmad/MM10_Annot_FromEncode.gtf
#gtfHG=/mnt/dataFast/ahrmad/HG38_Annot_FromEncode.gtf
gtfMM=/mnt/dataFast/ahrmad/MM39_Annot_From10x.gtf
gtfHG=/mnt/dataFast/ahrmad/HG38_Annot_From10x.gtf

featureCounts -T 64 -p -O -g gene_name -t gene -a ${gtfMM} -o H2RRNA_Mouse_gene_assigned -R BAM H2RRNA_Mouse_NameSortedGood_Tagged.bam
samtools sort --threads 64 -m 3G H2RRNA_Mouse_NameSortedGood_Tagged.bam.featureCounts.bam -o H2RRNA_Mouse_CoordSortedGood_Tagged_Assigned.bam
samtools index --threads 64 H2RRNA_Mouse_CoordSortedGood_Tagged_Assigned.bam
rm H2RRNA_Mouse_NameSortedGood_Tagged.bam.featureCounts.bam H2RRNA_Mouse_NameSortedGood_Tagged.bam

featureCounts -T 48 -p -O -g gene_name -t gene -a ${gtfHG} -o H2RRNA_Human_gene_assigned -R BAM H2RRNA_Human_NameSortedGood_Tagged.bam
samtools sort --threads 48 -m 3G H2RRNA_Human_NameSortedGood_Tagged.bam.featureCounts.bam -o H2RRNA_Human_CoordSortedGood_Tagged_Assigned.bam
samtools index --threads 64 H2RRNA_Human_CoordSortedGood_Tagged_Assigned.bam
rm H2RRNA_Human_NameSortedGood_Tagged.bam.featureCounts.bam H2RRNA_Human_NameSortedGood_Tagged.bam

exp=("H2RRNA")
org=("Human" "Mouse")

for e in "${exp[@]}"; do
    for o in "${org[@]}"; do
        echo ${e}_${o}
        # umi_tools count
        umi_tools count --stdin=${e}_${o}_CoordSortedGood_Tagged_Assigned.bam \
            --log=${e}_${o}_UMIToolsLogCount.txt \
            --wide-format-cell-counts \
            --extract-umi-method=tag \
            --umi-tag=UM \
            --edit-distance-threshold=1 \
            --spliced-is-unique \
            --per-cell --cell-tag=CB \
            --per-gene --gene-tag=XT --assigned-status-tag=XS \
            --paired \
            --random-seed=42 \
            --stdout=${e}_${o}_NoDupUMI_count.tsv

        # umi_tools dedup
        umi_tools dedup --stdin=${e}_${o}_CoordSortedGood_Tagged_Assigned.bam \
            --log=${e}_${o}_UMIToolsLogPerCellPerGene.txt \
            --output-stats=${e}_${o}_UMIToolsStatsPerCellPerGene.txt \
            --extract-umi-method=tag \
            --umi-tag=UM \
            --edit-distance-threshold=1 \
            --spliced-is-unique \
            --multimapping-detection-method=NH \
            --per-cell --cell-tag=CB \
            --per-gene --gene-tag=XT --assigned-status-tag=XS \
            --paired \
            --buffer-whole-contig \
            --random-seed=42 \
            --stdout=${e}_${o}_NoDupUMI_dedupPerCellPerGene.bam
    done
done

samtools index --threads 64 H2RRNA_Human_NoDupUMI_dedupPerCellPerGene.bam
samtools index --threads 64 H2RRNA_Mouse_NoDupUMI_dedupPerCellPerGene.bam

# Merging Alignments for stats
#samtools merge --threads 64 -o H2RD1_Human_MarkedDup.bam H2RD1_Human_H3K27ac_MarkedDup.bam H2RD1_Human_H3K27me3_MarkedDup.bam
#samtools index --threads 64 --bai --output H2RD1_Human_MarkedDup.bam.bai H2RD1_Human_MarkedDup.bam
#
#samtools merge --threads 64 -o H2RD1_Mouse_MarkedDup.bam H2RD1_Mouse_H3K27ac_MarkedDup.bam H2RD1_Mouse_H3K27me3_MarkedDup.bam
#samtools index --threads 64 --bai --output H2RD1_Mouse_MarkedDup.bam.bai H2RD1_Mouse_MarkedDup.bam

# Getting rid of Duplicate reads
samtools view --threads 64 --bam --with-header --require-flags 0x400 --output DUP.bam --unoutput H2RD1_Mouse_H3K27ac_NoDup.bam H2RD1_Mouse_H3K27ac_MarkedDup.bam
samtools index --threads 64 --bai --output H2RD1_Mouse_H3K27ac_NoDup.bam.bai H2RD1_Mouse_H3K27ac_NoDup.bam

samtools view --threads 64 --bam --with-header --require-flags 0x400 --output DUP.bam --unoutput H2RD1_Mouse_H3K27me3_NoDup.bam H2RD1_Mouse_H3K27me3_MarkedDup.bam
samtools index --threads 64 --bai --output H2RD1_Mouse_H3K27me3_NoDup.bam.bai H2RD1_Mouse_H3K27me3_NoDup.bam

samtools view --threads 64 --bam --with-header --require-flags 0x400 --output DUP.bam --unoutput H2RD1_Human_H3K27ac_NoDup.bam H2RD1_Human_H3K27ac_MarkedDup.bam
samtools index --threads 64 --bai --output H2RD1_Human_H3K27ac_NoDup.bam.bai H2RD1_Human_H3K27ac_NoDup.bam

samtools view --threads 64 --bam --with-header --require-flags 0x400 --output DUP.bam --unoutput H2RD1_Human_H3K27me3_NoDup.bam H2RD1_Human_H3K27me3_MarkedDup.bam
samtools index --threads 64 --bai --output H2RD1_Human_H3K27me3_NoDup.bam.bai H2RD1_Human_H3K27me3_NoDup.bam
rm DUP.bam


# Plot Making
python ./Make_Plots.py

macs3 callpeak --treatment H2RD1_Human_H3K27ac_CountFiltered.bam --format BAMPE \
	--name ac --outdir ac_macs3_50kbg \
	--gsize hs --keep-dup all \
	--llocal 50000 --min-length 600 --max-gap 300 \
	--cutoff-analysis

macs3 callpeak --treatment H2RD1_Human_H3K27me3_CountFiltered.bam --format BAMPE \
	--name me3 --outdir me3_macs3_q005_broad \
	--gsize hs --keep-dup all --broad \
	--llocal 1000000 --min-length 3000 --max-gap 1000 --qvalue 0.05 \
	--cutoff-analysis


# Making BigWig tracks
bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM --extendReads --centerReads -b H2RD1_Human_H3K27ac_NoDup.bam -o H2RD1_Human_H3K27ac_NoDup_RPKM.bw -of bigwig --effectiveGenomeSize 2913022398
bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM --extendReads --centerReads -b H2RD1_Human_H3K27me3_NoDup.bam -o H2RD1_Human_H3K27me3_NoDup_RPKM.bw -of bigwig --effectiveGenomeSize 2913022398

bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM --extendReads --centerReads -b H2RD1_Human_H3K27ac_Post.bam -o H2RD1_Human_H3K27ac_Post_RPKM.bw -of bigwig --effectiveGenomeSize 2913022398
bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM --extendReads --centerReads -b H2RD1_Human_H3K27me3_Post.bam -o H2RD1_Human_H3K27me3_Post_RPKM.bw -of bigwig --effectiveGenomeSize 2913022398

bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM -b H2RRNA_Human_NoDupUMI_dedupPerCellPerGene.bam -o H2RRNA_Human_NoDupUMI_dedupPerCellPerGene_RPKM.bw -of bigwig --effectiveGenomeSize 2913022398
bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM -b H2RRNA_Human_Post.bam -o H2RRNA_Human_Post_RPKM.bw -of bigwig --effectiveGenomeSize 2913022398

bamCoverage -p 64 -bs 10 --smoothLength 100 --normalizeUsing RPKM -b /mnt/data3/evotrack/Natalaya/k422/cellranger_pooled_output/pool_K422_DMSO_exp3000/outs/possorted_genome_bam.bam -o ScRNA_Natalya.bw -of bigwig --effectiveGenomeSize 2913022398

#RNA

#################################################
################## BUILD BAMs ###################
#################################################

# BASH
#exps=("H2RD1_Human_H3K27me3" "H2RD1_Human_H3K27ac")
#for e in "${exps[@]}"; do
#    samtools view --threads 64 --with-header --read-group-file ${e}_Post_cell_names.lst --output ${e}_Post.bam ${e}_NoDup.bam
#    samtools index --threads 64 ${e}_Post.bam
#done
#
#exps=("H2RRNA_Human")
#for e in "${exps[@]}"; do
#    samtools view --threads 64 --with-header --read-group-file ${e}_Post_cell_names.lst --output ${e}_Post.bam ${e}_NoDupUMI_dedupPerCellPerGene.bam
#    samtools index --threads 64 ${e}_Post.bam
#done

 #Build Strander Signal tracks from RNA


# Define the array of BAM file names
RNA_BAM_Names=( "H2RRNA_Human_Post.bam" "H2RRNA_Human_NoDupUMI_dedupPerCellPerGene.bam" "H2RRNA_Human_CoordSortedGood_Tagged_Assigned.bam" )
RNA_BAM_Names=( "10x_ScRNA.bam" )

refdir="/mnt/dataFast/ahrmad"

# Function to process each BAM file
process_bam() {
    RName=$1
    refdir=$2

    STAR --runMode inputAlignmentsFromBAM \
    --inputBAMfile ${RName} \
    --outWigType bedGraph \
    --outWigStrand Stranded \
    --outWigReferencesPrefix chr \
    --outFileNamePrefix ${RName}

    rm ${RName}Signal.Unique.str*.out.bg
    sortBed -i ${RName}Signal.UniqueMultiple.str1.out.bg > ${RName}Signal.Str1Sorted.bg
    sortBed -i ${RName}Signal.UniqueMultiple.str2.out.bg > ${RName}Signal.Str2Sorted.bg
    bedGraphToBigWig ${RName}Signal.Str1Sorted.bg ${refdir}/hg38.chrom.sizes ${RName}_minus.bw
    bedGraphToBigWig ${RName}Signal.Str2Sorted.bg ${refdir}/hg38.chrom.sizes ${RName}_plus.bw

    mv "${RName}_plus.bw" "${RName/.bam/}_plus.bw"
    mv "${RName}_minus.bw" "${RName/.bam/}_minus.bw"
    
    rm ${RName}Signal.Str*Sorted.bg ${RName}Signal.UniqueMultiple.str*.out.bg ${RName}Log.out
}

# Export the function and variables so they can be used by GNU Parallel
export -f process_bam
export refdir

# Use GNU Parallel to run the function in parallel
parallel --jobs 6 process_bam ::: "${RNA_BAM_Names[@]}" ::: "$refdir"





#NEW RNA LAUNCH
#SCH2R_202409

#!/bin/bash

set -euo pipefail

# Globally used variables
libName=H2R_TRISEQ_202409
outdir="/mnt/dataFast/ahrmad/ScH2R_202409/RNA2"
refdir="/mnt/dataFast/ahrmad"

# 5 Alignment and Duplicate detection

RNAsamples=("H2RRNA_Human")

# 
for RName in "${RNAsamples[@]}"; do
    echo "Aligning ${RName}..."
    /home/annan/H2R/AlignRNA.sh ${RName} \
        ${outdir}/${RName}_R1.fq.gz ${outdir}/${RName}_R2.fq.gz \
        ${refdir} ${outdir} 64 28
    
    echo "Finished aligning ${RName}"

done
