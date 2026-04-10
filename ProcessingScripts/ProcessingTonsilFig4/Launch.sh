#!/bin/bash

set -e

# Globally used variables
XDir=/mnt/dataFast/ahrmad/triseq_2025052_2
libName=H2R_TRISEQ_202505_2
outdir="${XDir}/processed"
refdir="/mnt/dataFast/ahrmad"
rawdir="${XDir}/raw"

### Tonsil Samples with 12 subsampling
DNA_Tonsil=( "Sc_VTD9_S1" "Sc_VTD11_S2")
RNA_Tonsil=( "Sc_VTr10_S2" "Sc_VTr12_S3" )


# Create output directory
mkdir -p ${outdir}

# 1 QC
# Fastqc MultiQC pre trimmd
mkdir -p ${XDir}/fastqc_pretrimming

fastqc --memory 4096 --threads 64 --outdir ${XDir}/fastqc_pretrimming ${XDir}/raw/*
multiqc --outdir ${XDir}/fastqc_pretrimming ${XDir}/fastqc_pretrimming/.



# 2 Tagging
# tags DNA: 
# SB: sample barcode: I2  -D BC_LEN=4 -D BC_START=14 REV
# MO: modality barcode I2 -D BC_LEN=8 -D BC_START=18 REV
# L1,L2,L3: ligation barcodes, CB: single-cell barcodes, RG: copy of CB -> I1 hardcoded

wl_lig="${XDir}/WLs/ligation_whitelist.txt"
wl_mod="${XDir}/WLs/dna_modality_whitelist.txt"
wl_sample8="${XDir}/WLs/sample8_whitelist.txt"
wl_sample12="${XDir}/WLs/sample12_whitelist.txt"

# Loop through each DNA name in the array
for DName in "${DNA_Tonsil[@]}"
do
    echo "Tagging ${DName}..."

    I1D1="${rawdir}/${DName}_I1_001.fastq.gz"
    I2D1="${rawdir}/${DName}_I2_001.fastq.gz"
    R1D1="${rawdir}/${DName}_R1_001.fastq.gz"
    R2D1="${rawdir}/${DName}_R2_001.fastq.gz"

    # Assign base according to DName
    case "$DName" in
        "Sc_VTD9_S1") BASE="A" ;;
        "Sc_VTD11_S2") BASE="T" ;;
        *) echo "Unknown DName: $DName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=14 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I2D1}" \
        "${R1D1}" \
        "${R2D1}" \
        "${wl_sample12}" \
        "${DName}" SB "${outdir}" first_pass_withBC_${BASE} rev

    # Second codon run
    codon run -plugin seq -release -D BC_LEN=8 -D BC_START=18 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I2D1}" \
        "${outdir}/${DName}_R1_001_SB.fastq" \
        "${outdir}/${DName}_R2_001_SB.fastq" \
        "${wl_mod}" \
        "${DName}" MO "${outdir}" not_first_pass rev

    # Third codon run
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/H2R/Tag_Lig3.codon \
        "${I1D1}" \
        "${outdir}/${DName}_R1_001_SB_MO.fastq" \
        "${outdir}/${DName}_R2_001_SB_MO.fastq" \
        "${wl_lig}" \
        "${DName}" CB "${outdir}"

    # Cleanup intermediate files
    rm "${outdir}/${DName}_R1_001_SB_MO.fastq" \
       "${outdir}/${DName}_R2_001_SB_MO.fastq" \
       "${outdir}/${DName}_R1_001_SB.fastq" \
       "${outdir}/${DName}_R2_001_SB.fastq"

    echo "Finished tagging ${DName}"

    echo "Trimming ${DName}..."
    trim_galore \
        --quality 10 \
        --cores 12 \
        --output_dir ${outdir} \
        --dont_gzip \
        --length 20 \
        --paired \
        ${outdir}/${DName}_R1_001_SB_MO_CB.fastq ${outdir}/${DName}_R2_001_SB_MO_CB.fastq

    rm ${outdir}/${DName}_R1_001_SB_MO_CB.fastq ${outdir}/${DName}_R2_001_SB_MO_CB.fastq

    echo "Finished trimming ${DName}"
done

# Loop through each RNA name in the array
for RName in "${RNA_Tonsil[@]}"
do
    echo "Tagging ${RName}..."

    I1R1="${rawdir}/${RName}_I1_001.fastq.gz"
    R1R1="${rawdir}/${RName}_R1_001.fastq.gz"
    R2R1="${rawdir}/${RName}_R2_001.fastq.gz"

    # Assign base according to RName
    case "$RName" in
        "Sc_VTr10_S2") BASE="A" ;;
        "Sc_VTr12_S3") BASE="T" ;;
        *) echo "Unknown RName: $RName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/H2R/Tag.codon \
        "${R2R1}" \
        "${R1R1}" \
        "${R2R1}" \
        "${wl_sample12}" \
        "${RName}" SB "${outdir}" first_pass_withBC_${BASE} rev

    # UM RNA1
    codon run -plugin seq -release -D BC_LEN=10 -D BC_START=4 /home/annan/H2R/Tag_UMI.codon \
        "${R2R1}" \
        "${outdir}/${RName}_R1_001_SB.fastq" \
        "${outdir}/${RName}_R2_001_SB.fastq" \
        "${RName}" UM "${outdir}"

    # CB RNA1
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/H2R/Tag_Lig3.codon \
        "${I1R1}" \
        "${outdir}/${RName}_R1_001_SB_UM.fastq" \
        "${outdir}/${RName}_R2_001_SB_UM.fastq" \
        "${wl_lig}" \
        "${RName}" CB "${outdir}"

    # Cleanup intermediate files
    rm "${outdir}/${RName}_R1_001_SB_UM.fastq" \
       "${outdir}/${RName}_R2_001_SB_UM.fastq" \
       "${outdir}/${RName}_R1_001_SB.fastq" \
       "${outdir}/${RName}_R2_001_SB.fastq"

    echo "Finished tagging ${RName}"
    
    echo "Trimming ${RName}..."
    trim_galore \
        --quality 10 \
        --cores 12 \
        --output_dir ${outdir} \
        --dont_gzip \
        --length 20 \
        --paired \
        ${outdir}/${RName}_R1_001_SB_UM_CB.fastq ${outdir}/${RName}_R2_001_SB_UM_CB.fastq
    
    rm ${outdir}/${RName}_R1_001_SB_UM_CB.fastq ${outdir}/${RName}_R2_001_SB_UM_CB.fastq

done

# 4 Read merging/splitting and SAM RG header buidling

# read / records merging
cat Sc_VTr10_S2_R1_001_SB_UM_CB_val_1.fq Sc_VTr12_S3_R1_001_SB_UM_CB_val_1.fq > Sc_VTr10_S2_Sc_VTr12_S3_R1_001_SB_UM_CB_val_1.fq
cat Sc_VTr10_S2_R2_001_SB_UM_CB_val_2.fq Sc_VTr12_S3_R2_001_SB_UM_CB_val_2.fq > Sc_VTr10_S2_Sc_VTr12_S3_R2_001_SB_UM_CB_val_2.fq
rm Sc_VTr10_S2_R1_001_SB_UM_CB_val_1.fq Sc_VTr12_S3_R1_001_SB_UM_CB_val_1.fq Sc_VTr10_S2_R2_001_SB_UM_CB_val_2.fq Sc_VTr12_S3_R2_001_SB_UM_CB_val_2.fq

cat Sc_VTD9_S1_R1_001_SB_MO_CB_val_1.fq Sc_VTD11_S2_R1_001_SB_MO_CB_val_1.fq > Sc_VTD9_S1_Sc_VTD11_S2_R1_001_SB_MO_CB_val_1.fq
cat Sc_VTD9_S1_R2_001_SB_MO_CB_val_2.fq Sc_VTD11_S2_R2_001_SB_MO_CB_val_2.fq > Sc_VTD9_S1_Sc_VTD11_S2_R2_001_SB_MO_CB_val_2.fq
rm Sc_VTD9_S1_R1_001_SB_MO_CB_val_1.fq Sc_VTD11_S2_R1_001_SB_MO_CB_val_1.fq Sc_VTD9_S1_R2_001_SB_MO_CB_val_2.fq Sc_VTD11_S2_R2_001_SB_MO_CB_val_2.fq

cat Tag_Records_Sc_VTr10_S2.tsv Tag_Records_Sc_VTr12_S3.tsv > Tag_Records_Sc_VTr10_S2_Sc_VTr12_S3.tsv
cat Tag_Records_Sc_VTD9_S1.tsv Tag_Records_Sc_VTD11_S2.tsv > Tag_Records_Sc_VTD9_S1_Sc_VTD11_S2.tsv
rm Tag_Records_Sc_VTr10_S2.tsv Tag_Records_Sc_VTr12_S3.tsv Tag_Records_Sc_VTD9_S1.tsv Tag_Records_Sc_VTD11_S2.tsv

pigz -p 64 "${outdir}"/Tag*tsv

#Tonsil DNA
DName=Sc_VTD9_S1_Sc_VTD11_S2
echo "Splitting ${DName}..."
codon run -plugin seq -release /home/annan/triseq052025_2/Split_Reads.codon \
    ${DName} \
    ${outdir} \
    ${libName} \
    exp1d ${outdir}/${DName}_R1_001_SB_MO_CB_val_1.fq ${outdir}/${DName}_R2_001_SB_MO_CB_val_2.fq
rm Sc_VTD9_S1_Sc_VTD11_S2_R1_001_SB_MO_CB_val_1.fq Sc_VTD9_S1_Sc_VTD11_S2_R2_001_SB_MO_CB_val_2.fq

#Tonsil RNA
RName=Sc_VTr10_S2_Sc_VTr12_S3
echo "Splitting ${RName}..."
codon run -plugin seq -release /home/annan/triseq052025_2/Split_Reads.codon \
    ${RName} \
    ${outdir} \
    ${libName} \
    exp1r ${outdir}/${RName}_R1_001_SB_UM_CB_val_1.fq ${outdir}/${RName}_R2_001_SB_UM_CB_val_2.fq
rm Sc_VTr10_S2_Sc_VTr12_S3_R1_001_SB_UM_CB_val_1.fq Sc_VTr10_S2_Sc_VTr12_S3_R2_001_SB_UM_CB_val_2.fq


# 5 Alignment and Duplicate detection

grcm39_eff_size=2654621783
hg38_eff_size=2913022398

gtfMM="${refdir}/MM39_Annot_From10x.gtf"
gtfHG="${refdir}/HG38_Annot_From10x.gtf"

bwaHG="${refdir}/hg38_bwamem2_index/hg38.fa"
bwaMM="${refdir}/GRCm39_bwamem2_index/GRCm39"
mods=( "H3K27ac" "H3K27me3" )


DNA_Name12_split_VTD9=( "Sc_VTD9_S1_Sc_VTD11_S2_VTHumanTonsil" "Sc_VTD9_S1_Sc_VTD11_S2_VTK422")
RNA_Name12_split_VTR10_11=( "Sc_VTr10_S2_Sc_VTr12_S3_VTHumanTonsil" "Sc_VTr10_S2_Sc_VTr12_S3_VTK422")

# Check reference files
required_refs=( "${gtfMM}" "${gtfHG}" "${bwaHG}.bwt.2bit.64" "${bwaMM}.bwt.2bit.64" "${refdir}/hg38-blacklist.v2.bed" "${refdir}/refdata-gex-GRCh38-2024-A/star" )
for ref in "${required_refs[@]}"; do
    if [ ! -e "${ref}" ]; then
        echo "ERROR: Required reference file not found: ${ref}"
        exit 1
    fi
done


## VTD
for DName in "${DNA_Name12_split_VTD9[@]}"; do
    for mod in "${mods[@]}"; do
        echo "Aligning ${DName}_${mod}..."
        /home/annan/H2R/AlignDNA.sh ${mod} ${DName} \
            ${outdir}/${DName}_${mod}_R1.fq \
            ${outdir}/${DName}_${mod}_R2.fq \
            ${refdir}/hg38-blacklist.v2.bed \
            ${outdir}/SAM_RG_Header_${DName}_${mod}.tsv \
            ${bwaHG} ${hg38_eff_size} ${outdir}
        
        /mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
            -I ${outdir}/${DName}_${mod}.bam -O ${outdir}/${DName}_${mod}_MarkedDup.bam \
            -M ${outdir}/${DName}_${mod}.DuplicateMetrics.txt \
            --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000
        
        samtools view --threads 64 --bam --with-header --require-flags 0x400 --output ${outdir}/DUP.bam --unoutput ${outdir}/${DName}_${mod}_NoDup.bam ${outdir}/${DName}_${mod}_MarkedDup.bam
        samtools index --threads 64 --bai --output ${outdir}/${DName}_${mod}_NoDup.bam.bai ${outdir}/${DName}_${mod}_NoDup.bam
        rm -f ${outdir}/DUP.bam ${outdir}/${DName}_${mod}.bam

        bamCoverage -p 64 -bs 100 --extendReads --centerReads -b ${outdir}/${DName}_${mod}_NoDup.bam -o ${outdir}/${DName}_${mod}_NoDup.bw -of bigwig --effectiveGenomeSize ${hg38_eff_size}

        echo "Finished aligning ${DName}"
    done
done

# REALIGNEMNT WITH STARSOLO

#!/bin/bash

set -euo pipefail

# Globally used variables
libName=H2R_TRISEQ_202505_2
outdir="/mnt/dataFast/ahrmad/triseq_2025052_2/RNA_reproc"
refdir="/mnt/dataFast/ahrmad"

# 5 Alignment and Duplicate detection

RNA_Name12_split_VTR10_11=("Sc_VTr10_S2_Sc_VTr12_S3_VTHumanTonsil")

# VTR
for RName in "${RNA_Name12_split_VTR10_11[@]}"; do
    echo "Aligning ${RName}..."
    /home/annan/H2R/AlignRNA.sh ${RName} \
        ${outdir}/${RName}_R1.fq ${outdir}/${RName}_R2.fq \
        ${refdir} ${outdir} 64 29
    
    echo "Finished aligning ${RName}"

done





mkdir -p ${outdir}/stats
mv ${outdir}/Reads_Per_Barcode_Sc_* ${outdir}/stats
mv ${outdir}/Barcode_Statistics_Sc_* ${outdir}/stats




#!/bin/bash

set -e

# Get seq data from EPFL GECF facility
#echo "ls AVT0320/*" | sftp -q -P 22 sftpgecf@svsftp.epfl.ch
#scp -r -P 22 sftpgecf@svsftp.epfl.ch:AVT0320/FilterOnR1/* .

#mv avid067*/* .


# Globally used variables
XDir=/mnt/dataFast/ahrmad/triseq_202508
libName=H2R_TRISEQ_202508
outdir="${XDir}/processed"
refdir="/mnt/dataFast/ahrmad"
rawdir="${XDir}/Raw_Reads"

### Tonsil Samples with 12 subsampling
DNA_Tonsil=( "Sc_VTD4T_S13" "Sc_VTD1_S3" )
RNA_Tonsil=( "Sc_VTR4F_S14" "Sc_VTr1_S4" )

bash ./merge2.sh

# 1 QC
# Fastqc MultiQC pre trimmd

mkdir -p ${XDir}/fastqc_pretrimming_unfiltered

fastqc --memory 4096 --threads 64 --outdir ${XDir}/fastqc_pretrimming_unfiltered ${rawdir}/*.fastq.gz
multiqc --outdir ${XDir}/fastqc_pretrimming_unfiltered ${XDir}/fastqc_pretrimming_unfiltered/.


# 2 Tagging
# tags DNA: 
# SB: sample barcode: I2  -D BC_LEN=4 -D BC_START=14 REV
# MO: modality barcode I2 -D BC_LEN=8 -D BC_START=18 REV
# L1,L2,L3: ligation barcodes, CB: single-cell barcodes, RG: copy of CB -> I1 hardcoded

wl_lig="${XDir}/WLs/ligation_whitelist.txt"
wl_mod="${XDir}/WLs/dna_modality_whitelist.txt"
wl_sample12="${XDir}/WLs/sample12_whitelist.txt"
wl_sample9="${XDir}/WLs/sample9_whitelist.txt"


# Loop through each DNA name in the array
for DName in "${DNA_Tonsil[@]}"
do
    echo "Tagging ${DName}..."

    I1D1="${rawdir}/${DName}_I1_001.fastq.gz"
    I2D1="${rawdir}/${DName}_I2_001.fastq.gz"
    R1D1="${rawdir}/${DName}_R1_001.fastq.gz"
    R2D1="${rawdir}/${DName}_R2_001.fastq.gz"

    # Assign base according to DName
    case "$DName" in
        "Sc_VTD4T_S13") BASE="A" ;;
        "Sc_VTD1_S3") BASE="T" ;;
        *) echo "Unknown DName: $DName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=14 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I2D1}" \
        "${R1D1}" \
        "${R2D1}" \
        "${wl_sample12}" \
        "${DName}" SB "${outdir}" first_pass_withBC_${BASE} rev

    # Second codon run
    codon run -plugin seq -release -D BC_LEN=8 -D BC_START=18 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I2D1}" \
        "${outdir}/${DName}_R1_001_SB.fastq" \
        "${outdir}/${DName}_R2_001_SB.fastq" \
        "${wl_mod}" \
        "${DName}" MO "${outdir}" not_first_pass rev

    # Third codon run
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/H2R/Tag_Lig3.codon \
        "${I1D1}" \
        "${outdir}/${DName}_R1_001_SB_MO.fastq" \
        "${outdir}/${DName}_R2_001_SB_MO.fastq" \
        "${wl_lig}" \
        "${DName}" CB "${outdir}"

    # Cleanup intermediate files
    rm "${outdir}/${DName}_R1_001_SB_MO.fastq" \
       "${outdir}/${DName}_R2_001_SB_MO.fastq" \
       "${outdir}/${DName}_R1_001_SB.fastq" \
       "${outdir}/${DName}_R2_001_SB.fastq"

    echo "Finished tagging ${DName}"

    echo "Trimming ${DName}..."
    trim_galore \
        --quality 10 \
        --cores 12 \
        --output_dir ${outdir} \
        --dont_gzip \
        --length 20 \
        --paired \
        ${outdir}/${DName}_R1_001_SB_MO_CB.fastq ${outdir}/${DName}_R2_001_SB_MO_CB.fastq

    rm ${outdir}/${DName}_R1_001_SB_MO_CB.fastq ${outdir}/${DName}_R2_001_SB_MO_CB.fastq

    echo "Finished trimming ${DName}"
done


# Loop through each RNA name in the array
for RName in "${RNA_Tonsil[@]}"
do
    echo "Tagging ${RName}..."

    I1R1="${rawdir}/${RName}_I1_001.fastq.gz"
    R1R1="${rawdir}/${RName}_R1_001.fastq.gz"
    R2R1="${rawdir}/${RName}_R2_001.fastq.gz"

    # Assign base according to RName
    case "$RName" in
        "Sc_VTR4F_S14") BASE="A" ;;
        "Sc_VTr1_S4") BASE="T" ;;
        *) echo "Unknown RNA RName: $RName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/H2R/Tag.codon \
        "${R2R1}" \
        "${R1R1}" \
        "${R2R1}" \
        "${wl_sample12}" \
        "${RName}" SB "${outdir}" first_pass_withBC_${BASE} rev

    # UM RNA1
    codon run -plugin seq -release -D BC_LEN=10 -D BC_START=4 /home/annan/H2R/Tag_UMI.codon \
        "${R2R1}" \
        "${outdir}/${RName}_R1_001_SB.fastq" \
        "${outdir}/${RName}_R2_001_SB.fastq" \
        "${RName}" UM "${outdir}"

    # CB RNA1
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/H2R/Tag_Lig3.codon \
        "${I1R1}" \
        "${outdir}/${RName}_R1_001_SB_UM.fastq" \
        "${outdir}/${RName}_R2_001_SB_UM.fastq" \
        "${wl_lig}" \
        "${RName}" CB "${outdir}"

    # Cleanup intermediate files
    rm "${outdir}/${RName}_R1_001_SB_UM.fastq" \
       "${outdir}/${RName}_R2_001_SB_UM.fastq" \
       "${outdir}/${RName}_R1_001_SB.fastq" \
       "${outdir}/${RName}_R2_001_SB.fastq"

    echo "Finished tagging ${RName}"
    
    echo "Trimming ${RName}..."
    trim_galore \
        --quality 10 \
        --cores 12 \
        --output_dir ${outdir} \
        --dont_gzip \
        --length 20 \
        --paired \
        ${outdir}/${RName}_R1_001_SB_UM_CB.fastq ${outdir}/${RName}_R2_001_SB_UM_CB.fastq
    
    rm ${outdir}/${RName}_R1_001_SB_UM_CB.fastq ${outdir}/${RName}_R2_001_SB_UM_CB.fastq

done


# read / records merging Tonsil
cat ${outdir}/Sc_VTD4T_S13_R1_001_SB_MO_CB_val_1.fq ${outdir}/Sc_VTD1_S3_R1_001_SB_MO_CB_val_1.fq > ${outdir}/Sc_VTD_R1_001_SB_MO_CB_val_1.fq
rm ${outdir}/Sc_VTD4T_S13_R1_001_SB_MO_CB_val_1.fq ${outdir}/Sc_VTD1_S3_R1_001_SB_MO_CB_val_1.fq
cat ${outdir}/Sc_VTD4T_S13_R2_001_SB_MO_CB_val_2.fq ${outdir}/Sc_VTD1_S3_R2_001_SB_MO_CB_val_2.fq > ${outdir}/Sc_VTD_R2_001_SB_MO_CB_val_2.fq
rm ${outdir}/Sc_VTD4T_S13_R2_001_SB_MO_CB_val_2.fq ${outdir}/Sc_VTD1_S3_R2_001_SB_MO_CB_val_2.fq

cat ${outdir}/Sc_VTR4F_S14_R1_001_SB_UM_CB_val_1.fq ${outdir}/Sc_VTr1_S4_R1_001_SB_UM_CB_val_1.fq > ${outdir}/Sc_VTR_R1_001_SB_UM_CB_val_1.fq
rm ${outdir}/Sc_VTR4F_S14_R1_001_SB_UM_CB_val_1.fq ${outdir}/Sc_VTr1_S4_R1_001_SB_UM_CB_val_1.fq
cat ${outdir}/Sc_VTR4F_S14_R2_001_SB_UM_CB_val_2.fq ${outdir}/Sc_VTr1_S4_R2_001_SB_UM_CB_val_2.fq > ${outdir}/Sc_VTR_R2_001_SB_UM_CB_val_2.fq
rm ${outdir}/Sc_VTR4F_S14_R2_001_SB_UM_CB_val_2.fq ${outdir}/Sc_VTr1_S4_R2_001_SB_UM_CB_val_2.fq

cat ${outdir}/Tag_Records_Sc_VTD4T_S13.tsv ${outdir}/Tag_Records_Sc_VTD1_S3.tsv > ${outdir}/Tag_Records_Sc_VTD.tsv
rm ${outdir}/Tag_Records_Sc_VTD4T_S13.tsv ${outdir}/Tag_Records_Sc_VTD1_S3.tsv

cat ${outdir}/Tag_Records_Sc_VTR4F_S14.tsv ${outdir}/Tag_Records_Sc_VTr1_S4.tsv > ${outdir}/Tag_Records_Sc_VTR.tsv
rm ${outdir}/Tag_Records_Sc_VTR4F_S14.tsv ${outdir}/Tag_Records_Sc_VTr1_S4.tsv

#Tonsil DNA
DName=Sc_VTD
echo "Splitting ${DName}..."
codon run -plugin seq -release /home/annan/triseq082025/Split_Reads2.codon \
    ${DName} \
    ${outdir} \
    ${libName} \
    exp1d ${outdir}/${DName}_R1_001_SB_MO_CB_val_1.fq ${outdir}/${DName}_R2_001_SB_MO_CB_val_2.fq
rm ${outdir}/${DName}_R1_001_SB_MO_CB_val_1.fq ${outdir}/${DName}_R2_001_SB_MO_CB_val_2.fq

#Tonsil RNA
RName=Sc_VTR
echo "Splitting ${RName}..."
codon run -plugin seq -release /home/annan/triseq082025/Split_Reads2.codon \
    ${RName} \
    ${outdir} \
    ${libName} \
    exp1r ${outdir}/${RName}_R1_001_SB_UM_CB_val_1.fq ${outdir}/${RName}_R2_001_SB_UM_CB_val_2.fq
rm ${outdir}/${RName}_R1_001_SB_UM_CB_val_1.fq ${outdir}/${RName}_R2_001_SB_UM_CB_val_2.fq



# 5 Alignment and Duplicate detection

grcm39_eff_size=2654621783
hg38_eff_size=2913022398

gtfMM="${refdir}/MM39_Annot_From10x.gtf"
gtfHG="${refdir}/HG38_Annot_From10x.gtf"

bwaHG="${refdir}/hg38_bwamem2_index/hg38.fa"
bwaMM="${refdir}/GRCm39_bwamem2_index/GRCm39"
mods=( "H3K27ac" "H3K27me3" )

DNA_Name12_split_VTD=( "Sc_VTD_HumanTonsil2" )
RNA_Name12_split_VTR=( "Sc_VTR_HumanTonsil2" )


# Check reference files
required_refs=( "${gtfMM}" "${gtfHG}" "${bwaHG}.bwt.2bit.64" "${bwaMM}.bwt.2bit.64" "${refdir}/hg38-blacklist.v2.bed" "${refdir}/refdata-gex-GRCh38-2024-A/star" )
for ref in "${required_refs[@]}"; do
    if [ ! -e "${ref}" ]; then
        echo "ERROR: Required reference file not found: ${ref}"
        exit 1
    fi
done



## VTD
for DName in "${DNA_Name12_split_VTD[@]}"; do
    for mod in "${mods[@]}"; do
        echo "Aligning ${DName}_${mod}..."
        /home/annan/H2R/AlignDNA.sh ${mod} ${DName} \
            ${outdir}/${DName}_${mod}_R1.fq \
            ${outdir}/${DName}_${mod}_R2.fq \
            ${refdir}/hg38-blacklist.v2.bed \
            ${outdir}/SAM_RG_Header_${DName}_${mod}.tsv \
            ${bwaHG} ${hg38_eff_size} ${outdir}
        
        pigz -p 40 ${outdir}/${DName}_${mod}_R1.fq ${outdir}/${DName}_${mod}_R2.fq

        /mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
            -I ${outdir}/${DName}_${mod}.bam -O ${outdir}/${DName}_${mod}_MarkedDup.bam \
            -M ${outdir}/${DName}_${mod}.DuplicateMetrics.txt \
            --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000
        
        samtools view --threads 64 --bam --with-header --require-flags 0x400 --output ${outdir}/DUP.bam --unoutput ${outdir}/${DName}_${mod}_NoDup.bam ${outdir}/${DName}_${mod}_MarkedDup.bam
        samtools index --threads 64 --bai --output ${outdir}/${DName}_${mod}_NoDup.bam.bai ${outdir}/${DName}_${mod}_NoDup.bam
        rm -f ${outdir}/DUP.bam ${outdir}/${DName}_${mod}.bam

        bamCoverage -p 64 -bs 100 --extendReads --centerReads -b ${outdir}/${DName}_${mod}_NoDup.bam -o ${outdir}/${DName}_${mod}_NoDup.bw -of bigwig --effectiveGenomeSize ${hg38_eff_size}

        echo "Finished aligning ${DName}"
    done
done


# VTR
for RName in "${RNA_Name12_split_VTR[@]}"; do
    echo "Aligning ${RName}..."
    /home/annan/H2R/AlignRNA.sh ${RName} \
        ${outdir}/${RName}_R1.fq ${outdir}/${RName}_R2.fq \
        ${refdir}/refdata-gex-GRCh38-2024-A/star ${outdir} 64

    pigz -p 40 ${outdir}/${RName}_R1.fq ${outdir}/${RName}_R2.fq
    
    featureCounts -T 64 -p -O -g gene_name -t gene -a ${gtfHG} -o ${outdir}/${RName}_gene_assigned -R BAM ${outdir}/${RName}_NameSortedGood_Tagged.bam
    samtools sort --threads 64 -m 3G ${outdir}/${RName}_NameSortedGood_Tagged.bam.featureCounts.bam -o ${outdir}/${RName}_ALL.bam
    samtools index --threads 64 --bai --output ${outdir}/${RName}_ALL.bam.bai ${outdir}/${RName}_ALL.bam
    rm ${outdir}/${RName}_NameSortedGood_Tagged.bam.featureCounts.bam ${outdir}/${RName}_NameSortedGood_Tagged.bam

    umi_tools count --stdin=${outdir}/${RName}_ALL.bam \
    --log=${outdir}/${RName}_UMIToolsLogCount.txt \
    --wide-format-cell-counts \
    --extract-umi-method=tag \
    --umi-tag=UM \
    --edit-distance-threshold=1 \
    --spliced-is-unique \
    --per-cell --cell-tag=CB \
    --per-gene --gene-tag=XT --assigned-status-tag=XS \
    --paired \
    --random-seed=42 \
    --stdout=${outdir}/${RName}_NoDupUMI_count.tsv

    # umi_tools dedup
    umi_tools dedup --stdin=${outdir}/${RName}_ALL.bam \
        --log=${outdir}/${RName}_UMIToolsLogPerCellPerGene.txt \
        --output-stats=${outdir}/${RName}_UMIToolsStatsPerCellPerGene.txt \
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
        --stdout=${outdir}/${RName}_UMI.bam
    
    samtools index --threads 64 --bai --output ${outdir}/${RName}_UMI.bam.bai ${outdir}/${RName}_UMI.bam

    #bamCoverage -p 64 -bs 100 -b ${outdir}/${RName}_UMI.bam -o ${outdir}/${RName}_UMI.bw -of bigwig --effectiveGenomeSize ${hg38_eff_size}
    #bamCoverage -p 64 -bs 100 -b ${outdir}/${RName}_ALL.bam -o ${outdir}/${RName}_ALL.bw -of bigwig --effectiveGenomeSize ${hg38_eff_size}

    STAR --runMode inputAlignmentsFromBAM \
    --inputBAMfile ${outdir}/${RName}_ALL.bam \
    --outWigType bedGraph \
    --outWigStrand Stranded \
    --outWigReferencesPrefix chr \
    --outFileNamePrefix ${outdir}/${RName}

    rm ${outdir}/${RName}Signal.Unique.str*.out.bg
    sortBed -i ${outdir}/${RName}Signal.UniqueMultiple.str1.out.bg > ${outdir}/${RName}Signal.Str1Sorted.bg
    sortBed -i ${outdir}/${RName}Signal.UniqueMultiple.str2.out.bg > ${outdir}/${RName}Signal.Str2Sorted.bg
    bedGraphToBigWig ${outdir}/${RName}Signal.Str1Sorted.bg ${refdir}/hg38.chrom.sizes ${outdir}/${RName}_minus.bw
    bedGraphToBigWig ${outdir}/${RName}Signal.Str2Sorted.bg ${refdir}/hg38.chrom.sizes ${outdir}/${RName}_plus.bw
    
    #mv "${outdir}/${RName}_plus.bw" "${outdir}/${RName/.bam/}_plus.bw"
    #mv "${outdir}/${RName}_minus.bw" "${outdir}/${RName/.bam/}_minus.bw"
    echo "Finished aligning ${RName}"

done



mkdir -p ${outdir}/stats
mv ${outdir}/Reads_Per_Barcode_Sc_* ${outdir}/stats
mv ${outdir}/Barcode_Statistics_Sc_* ${outdir}/stats
python ./mm.py

pigz -p 48 "${outdir}"/Tag*tsv

#RNA2

#!/bin/bash
set -euo pipefail

threads=72
refdir="/mnt/dataFast/ahrmad"
outdir="/mnt/dataFast/ahrmad/Tonsil2/RNA2"

RNAsamples=(
  "Sc_VTR_HumanTonsil2"
)
species="human"

for RName in "${RNAsamples[@]}"; do
  echo "Aligning ${RName}..."

  R1="${outdir}/${RName}_R1.fq.gz"
  R2="${outdir}/${RName}_R2.fq.gz"
  USAM="${outdir}/${RName}_tagged.usam"
  UBAM="${outdir}/${RName}_tagged.ubam"

  if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "WARNING: Missing merged RNA FASTQs for ${RName}, skipping (${R1}, ${R2})" >&2
    continue
  fi

  # FASTQ → tagged unmapped SAM 
  codon run -plugin seq -release /home/annan/H2R/FqToSAM.codon "${R1}" "${R2}" "${USAM}"
  rm "${R1}" "${R2}"

  # Do NOT rm -f unless you want to delete merged FASTQs
  # rm -f "${R1}" "${R2}"

  /home/annan/H2R/AlignRNA.sh "${RName}" "${USAM}" "${refdir}" "${outdir}" "${threads}" "${species}"
  rm -f "${USAM}"

  samtools view --threads "${threads}" -b -o "${UBAM}" "${USAM}"
  samtools index --threads "${threads}" "${outdir}/${RName}.filtered_cells.bam"

  echo "Finished aligning ${RName}"
done



#!/bin/bash
set -euo pipefail

# Get seq data from EPFL GECF facility
# AVT0419  AVT0433

# Globally used variables
XDir=/mnt/dataFast/ahrmad/triseq_202601
libName=H2R_TRISEQ_202601
outdir="${XDir}/processed"
refdir="/mnt/dataFast/ahrmad"
rawdir="/mnt/data3/ahrmad/NewSeq/raw"

# DT
DNA=( "Sc_TP3D_1_S5" "Sc_TP3D_3_S6" "Sc_TP4D_1_S9" "Sc_TP4D_2_S10" )
RNA=( "Sc_TP3r_1_S7" "Sc_TP3r_3_S8" "Sc_TP4r_1_S11" "Sc_TP4r_2_S12" )


# 1 QC
# Fastqc MultiQC pre trimmd
# Create output directory
mkdir -p ${outdir}

#mkdir -p ${XDir}/fastqc_pretrimming
#fastqc --memory 4096 --threads 72 --outdir ${XDir}/fastqc_pretrimming ${rawdir}/*.fastq.gz
#multiqc --outdir ${XDir}/fastqc_pretrimming ${XDir}/fastqc_pretrimming/.


# 2 Tagging
# tags DNA: 
# SB: sample barcode: I2  -D BC_LEN=4 -D BC_START=14 REV
# MO: modality barcode I2 -D BC_LEN=8 -D BC_START=18 REV
# L1,L2,L3: ligation barcodes, CB: single-cell barcodes, RG: copy of CB -> I1 hardcoded

wl_lig="${XDir}/WLs/ligation_whitelist.txt"
wl_sb_TP3="${XDir}/WLs/dna_sample_whitelist_TP3.txt"
wl_sb_TP4="${XDir}/WLs/dna_sample_whitelist_TP4.txt"
WL_MOD_USE="${XDir}/WLs/dna_modality_whitelist.txt"

wl_sb_TP3r="${XDir}/WLs/rna_sample_whitelist_TP3.txt"
wl_sb_TP4r="${XDir}/WLs/rna_sample_whitelist_TP4.txt"

# Loop through each DNA name in the array
for DName in "${DNA[@]}"
do
    echo "Tagging ${DName}..."

    I1D1="${rawdir}/${DName}_I1_001.fastq.gz"
    R1D1="${rawdir}/${DName}_R1_001.fastq.gz"
    R2D1="${rawdir}/${DName}_R2_001.fastq.gz"

    # Assign base according to DName
    case "$DName" in
        "Sc_TP3D_1_S5")  BASE="A" ; WL_SB_USE="${wl_sb_TP3}" ;;
        "Sc_TP3D_3_S6") BASE="T" ; WL_SB_USE="${wl_sb_TP3}" ;;
        "Sc_TP4D_1_S9") BASE="A" ; WL_SB_USE="${wl_sb_TP4}" ;;
        "Sc_TP4D_2_S10") BASE="T" ; WL_SB_USE="${wl_sb_TP4}" ;;
        *) echo "Unknown DName: $DName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=3 -D BC_START=0 -D HD=0 /home/annan/H2R/Tag.codon \
        "${I1D1}" \
        "${R1D1}" \
        "${R2D1}" \
        "${WL_SB_USE}" \
        "${DName}" SB "${outdir}" first_pass_withBC_${BASE} fw

    # Second codon run (use wl_mod_atac only for Sc_ATD_1_S3)
    codon run -plugin seq -release -D BC_LEN=8 -D BC_START=3 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I1D1}" \
        "${outdir}/${DName}_R1_001_SB.fastq" \
        "${outdir}/${DName}_R2_001_SB.fastq" \
        "${WL_MOD_USE}" \
        "${DName}" MO "${outdir}" not_first_pass fw

    rm "${outdir}/${DName}_R1_001_SB.fastq" \
       "${outdir}/${DName}_R2_001_SB.fastq"

    # Third codon run
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/H2R/Tag_Lig3_DT.codon \
        "${I1D1}" \
        "${outdir}/${DName}_R1_001_SB_MO.fastq" \
        "${outdir}/${DName}_R2_001_SB_MO.fastq" \
        "${wl_lig}" \
        "${DName}" CB "${outdir}"

    # Cleanup intermediate files
    rm "${outdir}/${DName}_R1_001_SB_MO.fastq" \
       "${outdir}/${DName}_R2_001_SB_MO.fastq"

    echo "Finished tagging ${DName}"

    echo "Trimming ${DName}..."
    trim_galore \
        --quality 10 \
        --cores 12 \
        --output_dir ${outdir} \
        --dont_gzip \
        --length 20 \
        --paired \
        ${outdir}/${DName}_R1_001_SB_MO_CB.fastq ${outdir}/${DName}_R2_001_SB_MO_CB.fastq

    rm ${outdir}/${DName}_R1_001_SB_MO_CB.fastq ${outdir}/${DName}_R2_001_SB_MO_CB.fastq

    echo "Finished trimming ${DName}"
done


for RName in "${RNA[@]}"
do
    echo "Tagging ${RName}..."

    I1R1="${rawdir}/${RName}_I1_001.fastq.gz"
    R1R1="${rawdir}/${RName}_R1_001.fastq.gz"
    R2R1="${rawdir}/${RName}_R2_001.fastq.gz"

    # Assign base according to DName
    case "$RName" in
        "Sc_TP3r_1_S7")  BASE="A" ; WL_SB_USE="${wl_sb_TP3r}" ;;
        "Sc_TP3r_3_S8") BASE="T" ; WL_SB_USE="${wl_sb_TP3r}" ;;
        "Sc_TP4r_1_S11") BASE="A" ; WL_SB_USE="${wl_sb_TP4r}" ;;
        "Sc_TP4r_2_S12") BASE="T" ; WL_SB_USE="${wl_sb_TP4r}" ;;
        *) echo "Unknown DName: $DName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/H2R/Tag.codon \
        "${R2R1}" \
        "${R1R1}" \
        "${R2R1}" \
        "${WL_SB_USE}" \
        "${RName}" SB "${outdir}" first_pass_withBC_${BASE} rev

    # UM RNA1
    codon run -plugin seq -release -D BC_LEN=10 -D BC_START=4 /home/annan/H2R/Tag_UMI.codon \
        "${R2R1}" \
        "${outdir}/${RName}_R1_001_SB.fastq" \
        "${outdir}/${RName}_R2_001_SB.fastq" \
        "${RName}" UM "${outdir}"

    rm "${outdir}/${RName}_R1_001_SB.fastq" \
       "${outdir}/${RName}_R2_001_SB.fastq"

    # CB RNA1
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/H2R/Tag_Lig3.codon \
        "${I1R1}" \
        "${outdir}/${RName}_R1_001_SB_UM.fastq" \
        "${outdir}/${RName}_R2_001_SB_UM.fastq" \
        "${wl_lig}" \
        "${RName}" CB "${outdir}"

    # Cleanup intermediate files
    rm "${outdir}/${RName}_R1_001_SB_UM.fastq" \
       "${outdir}/${RName}_R2_001_SB_UM.fastq"

    echo "Finished tagging ${RName}"

    echo "Trimming ${RName}..."
    trim_galore \
        --quality 10 \
        --cores 12 \
        --output_dir ${outdir} \
        --dont_gzip \
        --length 20 \
        --paired \
        ${outdir}/${RName}_R1_001_SB_UM_CB.fastq ${outdir}/${RName}_R2_001_SB_UM_CB.fastq
    
    rm ${outdir}/${RName}_R1_001_SB_UM_CB.fastq ${outdir}/${RName}_R2_001_SB_UM_CB.fastq

    echo "Finished trimming ${RName}"

done



for s in "${RNA[@]}"; do
  codon run -plugin seq -release /home/annan/H2R/Split_ReadsV2.codon \
    "$s" "${outdir}" "${libName}" rna \
    - \
    "${outdir}/${s}_R1_001_SB_UM_CB_val_1.fq.gz" \
    "${outdir}/${s}_R2_001_SB_UM_CB_val_2.fq.gz" \
    "${XDir}/WLs/sb_map_tonsil_RNA.tsv"
done


for s in "${DNA[@]}"; do
  codon run -plugin seq -release /home/annan/H2R/Split_ReadsV2.codon \
    "$s" "${outdir}" "${libName}" dna \
    "${XDir}/WLs/mo_map_tonsil.tsv" \
    "${outdir}/${s}_R1_001_SB_MO_CB_val_1.fq.gz" \
    "${outdir}/${s}_R2_001_SB_MO_CB_val_2.fq.gz" \
    "${XDir}/WLs/sb_map_tonsil.tsv"
done



# ---- helper: concat gz fastqs safely (gzip streams concatenate) ----
cat_gz_pair() {
  local in1="$1"
  local in2="$2"
  local out="$3"
  # require both to exist (hard fail to avoid silent partial merges)
  [[ -e "$in1" ]] || { echo "ERROR: missing $in1" >&2; exit 1; }
  [[ -e "$in2" ]] || { echo "ERROR: missing $in2" >&2; exit 1; }
  cat "$in1" "$in2" > "$out"
}

# ---- helper: merge RG headers (unique lines) ----
merge_rg() {
  local in1="$1"
  local in2="$2"
  local out="$3"
  [[ -e "$in1" ]] || { echo "ERROR: missing $in1" >&2; exit 1; }
  [[ -e "$in2" ]] || { echo "ERROR: missing $in2" >&2; exit 1; }
  cat "$in1" "$in2" | sort -u > "$out"
}

# ------------------------------------------------------------
# RNA merges
# ------------------------------------------------------------
merge_rna_group() {
  local s1="$1" s2="$2" merged_prefix="$3" group="$4"
  echo "[RNA] ${merged_prefix}_${group}"

  cat_gz_pair \
    "${outdir}/${s1}_${group}_R1.fq.gz" \
    "${outdir}/${s2}_${group}_R1.fq.gz" \
    "${outdir}/${merged_prefix}_${group}_R1.fq.gz"

  cat_gz_pair \
    "${outdir}/${s1}_${group}_R2.fq.gz" \
    "${outdir}/${s2}_${group}_R2.fq.gz" \
    "${outdir}/${merged_prefix}_${group}_R2.fq.gz"

  merge_rg \
    "${outdir}/SAM_RG_Header_${s1}_${group}.tsv" \
    "${outdir}/SAM_RG_Header_${s2}_${group}.tsv" \
    "${outdir}/SAM_RG_Header_${merged_prefix}_${group}.tsv"
}

# TP3 RNA
merge_rna_group "Sc_TP3r_1_S7" "Sc_TP3r_3_S8" "Sc_TP3r" "TP3_A"
merge_rna_group "Sc_TP3r_1_S7" "Sc_TP3r_3_S8" "Sc_TP3r" "TP3_B"

# TP4 RNA
merge_rna_group "Sc_TP4r_1_S11" "Sc_TP4r_2_S12" "Sc_TP4r" "TP4_A"
merge_rna_group "Sc_TP4r_1_S11" "Sc_TP4r_2_S12" "Sc_TP4r" "TP4_B"


# ------------------------------------------------------------
# DNA merges (per sb_group + mark)
# Marks are derived from your mo_map:
# TP3_A: H3K27me3, H3K9me3
# TP3_B: H3K27me3, H3K27ac
# TP4_A: H3K27me3, H3K27ac
# TP4_B: H3K27me3, H3K9me3
# ------------------------------------------------------------
merge_dna_group_mark() {
  local s1="$1" s2="$2" merged_prefix="$3" group="$4" mark="$5"
  echo "[DNA] ${merged_prefix}_${group}_${mark}"

  cat_gz_pair \
    "${outdir}/${s1}_${group}_${mark}_R1.fq.gz" \
    "${outdir}/${s2}_${group}_${mark}_R1.fq.gz" \
    "${outdir}/${merged_prefix}_${group}_${mark}_R1.fq.gz"

  cat_gz_pair \
    "${outdir}/${s1}_${group}_${mark}_R2.fq.gz" \
    "${outdir}/${s2}_${group}_${mark}_R2.fq.gz" \
    "${outdir}/${merged_prefix}_${group}_${mark}_R2.fq.gz"

  merge_rg \
    "${outdir}/SAM_RG_Header_${s1}_${group}_${mark}.tsv" \
    "${outdir}/SAM_RG_Header_${s2}_${group}_${mark}.tsv" \
    "${outdir}/SAM_RG_Header_${merged_prefix}_${group}_${mark}.tsv"
}

# TP3 DNA
# Inputs: Sc_TP3D_1_S5 + Sc_TP3D_3_S6  -> merged prefix "Sc_TP3D"
merge_dna_group_mark "Sc_TP3D_1_S5" "Sc_TP3D_3_S6" "Sc_TP3D" "TP3_A" "H3K27me3"
merge_dna_group_mark "Sc_TP3D_1_S5" "Sc_TP3D_3_S6" "Sc_TP3D" "TP3_A" "H3K9me3"
merge_dna_group_mark "Sc_TP3D_1_S5" "Sc_TP3D_3_S6" "Sc_TP3D" "TP3_B" "H3K27me3"
merge_dna_group_mark "Sc_TP3D_1_S5" "Sc_TP3D_3_S6" "Sc_TP3D" "TP3_B" "H3K27ac"

# TP4 DNA
# Inputs: Sc_TP4D_1_S9 + Sc_TP4D_2_S10 -> merged prefix "Sc_TP4D"
merge_dna_group_mark "Sc_TP4D_1_S9" "Sc_TP4D_2_S10" "Sc_TP4D" "TP4_A" "H3K27me3"
merge_dna_group_mark "Sc_TP4D_1_S9" "Sc_TP4D_2_S10" "Sc_TP4D" "TP4_A" "H3K27ac"
merge_dna_group_mark "Sc_TP4D_1_S9" "Sc_TP4D_2_S10" "Sc_TP4D" "TP4_B" "H3K27me3"
merge_dna_group_mark "Sc_TP4D_1_S9" "Sc_TP4D_2_S10" "Sc_TP4D" "TP4_B" "H3K9me3"

echo "DONE: merged FASTQs and SAM_RG_Header TSVs into ${outdir}"




# 5 Alignment and Duplicate detection
grcm39_eff_size=2654621783
hg38_eff_size=2913022398

bwaHG="${refdir}/hg38_bwamem2_index/hg38.fa"
bwaMM="${refdir}/GRCm39_bwamem2_index/GRCm39"

blacklistHG="${refdir}/hg38-blacklist.v2.bed"
blacklistMM="${refdir}/mm10_grcm39LiftOver-blacklist.v2.bed"

MO_MAP="${XDir}/WLs/mo_map_tonsil.tsv"
SB_MAP="${XDir}/WLs/sb_map_tonsil.tsv"
threads=72

# ---- merged DNA targets (tonsil) ----
# Species: set once per merged cohort (edit if needed)
declare -A DNA_MERGED_SPECIES=(
  [Sc_TP3D]=human
  [Sc_TP4D]=human
)

# Groups per merged cohort
declare -A DNA_GROUPS=(
  [Sc_TP3D]="TP3_A TP3_B"
  [Sc_TP4D]="TP4_A TP4_B"
)

# Marks per (cohort,group) from your 4-col mo_map_tonsil.tsv
declare -A DNA_MARKS=(
  [Sc_TP3D:TP3_A]="H3K27me3 H3K9me3"
  [Sc_TP3D:TP3_B]="H3K27me3 H3K27ac"
  [Sc_TP4D:TP4_A]="H3K27me3 H3K27ac"
  [Sc_TP4D:TP4_B]="H3K27me3 H3K9me3"
)

for cohort in Sc_TP3D Sc_TP4D; do
  sp="${DNA_MERGED_SPECIES[${cohort}]:-}"
  if [[ -z "${sp}" ]]; then
    echo "ERROR: No species set for merged cohort ${cohort} in DNA_MERGED_SPECIES[]" >&2
    exit 1
  fi

  if [[ "${sp}" == "mouse" ]]; then
    bwaRef="${bwaMM}"
    effSize="${grcm39_eff_size}"
    blacklist="${blacklistMM}"
  elif [[ "${sp}" == "human" ]]; then
    bwaRef="${bwaHG}"
    effSize="${hg38_eff_size}"
    blacklist="${blacklistHG}"
  else
    echo "ERROR: species for ${cohort} must be human|mouse (got: ${sp})" >&2
    exit 1
  fi

  for grp in ${DNA_GROUPS[${cohort}]}; do
    key="${cohort}:${grp}"
    for mod in ${DNA_MARKS[${key}]}; do
      # merged FASTQs are gz
      R1="${outdir}/${cohort}_${grp}_${mod}_R1.fq.gz"
      R2="${outdir}/${cohort}_${grp}_${mod}_R2.fq.gz"
      RG="${outdir}/SAM_RG_Header_${cohort}_${grp}_${mod}.tsv"

      if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "WARNING: Missing FASTQs for ${cohort}_${grp}_${mod}, skipping (${R1}, ${R2})" >&2
        continue
      fi
      if [[ ! -s "${RG}" ]]; then
        echo "WARNING: Missing RG header for ${cohort}_${grp}_${mod}, skipping (${RG})" >&2
        continue
      fi

      echo "Aligning ${cohort}_${grp}_${mod} (${sp})..."
      /home/annan/H2R/AlignDNA.sh "${mod}" "${cohort}_${grp}" \
        "${R1}" \
        "${R2}" \
        "${blacklist}" \
        "${RG}" \
        "${bwaRef}" "${effSize}" "${outdir}"

      /mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
        -I "${outdir}/${cohort}_${grp}_${mod}.bam" -O "${outdir}/${cohort}_${grp}_${mod}_MarkedDup.bam" \
        -M "${outdir}/${cohort}_${grp}_${mod}.DuplicateMetrics.txt" \
        --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000

      samtools view --threads "${threads}" --bam --with-header \
        --require-flags 0x400 \
        --output "${outdir}/DUP.bam" \
        --unoutput "${outdir}/${cohort}_${grp}_${mod}_NoDup.bam" \
        "${outdir}/${cohort}_${grp}_${mod}_MarkedDup.bam"

      samtools index --threads "${threads}" --bai \
        --output "${outdir}/${cohort}_${grp}_${mod}_NoDup.bam.bai" \
        "${outdir}/${cohort}_${grp}_${mod}_NoDup.bam"

      rm -f "${outdir}/DUP.bam" "${outdir}/${cohort}_${grp}_${mod}.bam"

      bamCoverage -p "${threads}" -bs 100 --extendReads --centerReads \
        -b "${outdir}/${cohort}_${grp}_${mod}_NoDup.bam" \
        -o "${outdir}/${cohort}_${grp}_${mod}_NoDup.bw" \
        -of bigwig --effectiveGenomeSize "${effSize}"

      echo "Finished aligning ${cohort}_${grp}_${mod}"
    done
  done
done


RNAsamples=( "Sc_TP3r_TP3_A" "Sc_TP3r_TP3_B" "Sc_TP4r_TP4_A" "Sc_TP4r_TP4_B" )
species="human"

for RName in "${RNAsamples[@]}"; do
  echo "Aligning ${RName}..."

  R1="${outdir}/${RName}_R1.fq.gz"
  R2="${outdir}/${RName}_R2.fq.gz"
  USAM="${outdir}/${RName}_tagged.usam"
  UBAM="${outdir}/${RName}_tagged.ubam"

  if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "WARNING: Missing merged RNA FASTQs for ${RName}, skipping (${R1}, ${R2})" >&2
    continue
  fi

  # FASTQ → tagged unmapped SAM (must accept gz; if your FqToSAM.codon doesn't, tell me and I'll patch it)
  codon run -plugin seq -release /home/annan/H2R/FqToSAM.codon "${R1}" "${R2}" "${USAM}"

  # Do NOT rm -f unless you want to delete merged FASTQs
  # rm -f "${R1}" "${R2}"

  /home/annan/H2R/AlignRNA.sh "${RName}" "${USAM}" "${refdir}" "${outdir}" "${threads}" "${species}"

  samtools view --threads "${threads}" -b -o "${UBAM}" "${USAM}"
  samtools index --threads "${threads}" "${outdir}/${RName}.filtered_cells.bam"
  rm -f "${USAM}"

  echo "Finished aligning ${RName}"
done


mkdir -p ${outdir}/stats
mv ${outdir}/Reads_Per_Barcode_Sc_* ${outdir}/stats
mv ${outdir}/Barcode_Statistics_Sc_* ${outdir}/stats

