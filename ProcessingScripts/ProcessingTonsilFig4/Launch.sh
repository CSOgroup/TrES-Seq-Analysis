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





