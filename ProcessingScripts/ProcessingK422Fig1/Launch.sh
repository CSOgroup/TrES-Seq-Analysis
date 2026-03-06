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
