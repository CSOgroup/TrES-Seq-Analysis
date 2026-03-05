#!/bin/bash

set -e

# Get seq data from EPFL GECF facility
#echo "ls AVT0266/*" | sftp -q -P 22 sftpgecf@svsftp.epfl.ch
#scp -r -P 22 sftpgecf@svsftp.epfl.ch:AVT0266/* .

# Globally used variables
XDir=/mnt/dataFast/ahrmad/triseq_2025052
libName=H2R_TRISEQ_202505_2
outdir="${XDir}/processed"
refdir="/mnt/dataFast/ahrmad"
rawdir="${XDir}/raw"

### Tonsil Samples with 12 subsampling
DNA_Tonsil=( "Sc_VTD9_S1" )
RNA_Tonsil=( "Sc_VTr10_S2" "Sc_VTr11_S3" )

### Cell line Samples with 8 subsampling
DNA_CellLine=( "Sc_CMD1_S4" "Sc_CMD2_S6" "Sc_CMD3_S8" )
RNA_CellLine=( "Sc_CMr1_S5" "Sc_CMr2_S7" "Sc_CMr3_S9" )



# Create output directory
mkdir -p ${outdir}

# 0 Env setup
#mm env create -n t2
# Launch
#mm install samtools screen bwa-mem2 pandas polars ipython pysam pybedtools numpy matplotlib seaborn scipy pyarrow star==2.7.1a umi_tools==1.1.5 fastqc multiqc trim-galore deeptools subread coreutils parallel upsetplot ucsc-bedGraphToBigWig
# analysis
#mm install snapatac2 anndata scanpy matplotlib-venn leidenalg bbknn


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
for DName in "${DNA_CellLine[@]}"
do
    echo "Tagging ${DName}..."

    I1D1="${rawdir}/${DName}_I1_001.fastq.gz"
    I2D1="${rawdir}/${DName}_I2_001.fastq.gz"
    R1D1="${rawdir}/${DName}_R1_001.fastq.gz"
    R2D1="${rawdir}/${DName}_R2_001.fastq.gz"

    # Assign base according to DName
    case "$DName" in
        "Sc_CMD1_S4") BASE="A" ;;
        "Sc_CMD2_S6") BASE="T" ;;
        "Sc_CMD3_S8") BASE="C" ;;
        *) echo "Unknown DName: $DName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=14 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I2D1}" \
        "${R1D1}" \
        "${R2D1}" \
        "${wl_sample8}" \
        "${DName}" SB "${outdir}" first_pass_withBC_${BASE}  rev

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
for RName in "${RNA_CellLine[@]}"
do
    echo "Tagging ${RName}..."

    I1R1="${rawdir}/${RName}_I1_001.fastq.gz"
    R1R1="${rawdir}/${RName}_R1_001.fastq.gz"
    R2R1="${rawdir}/${RName}_R2_001.fastq.gz"

    # Assign base according to RName
    case "$RName" in
        "Sc_CMr1_S5") BASE="A" ;;
        "Sc_CMr2_S7") BASE="T" ;;
        "Sc_CMr3_S9") BASE="C" ;;
        *) echo "Unknown RNA RName: $RName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/H2R/Tag.codon \
        "${R2R1}" \
        "${R1R1}" \
        "${R2R1}" \
        "${wl_sample8}" \
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

    echo "Finished trimming ${RName}"

done


# 4 Read merging/splitting and SAM RG header buidling


### Cell line Samples with 8 subsampling
DNA_CellLine=( "Sc_CMD1_S4" "Sc_CMD2_S6" "Sc_CMD3_S8" )
RNA_CellLine=( "Sc_CMr1_S5" "Sc_CMr2_S7" "Sc_CMr3_S9" )

# Cell line 
cat Sc_CMD1_S4_R1_001_SB_MO_CB_val_1.fq Sc_CMD2_S6_R1_001_SB_MO_CB_val_1.fq Sc_CMD3_S8_R1_001_SB_MO_CB_val_1.fq > Sc_CMD_R1_001_SB_MO_CB_val_1.fq
cat Sc_CMD1_S4_R2_001_SB_MO_CB_val_2.fq Sc_CMD2_S6_R2_001_SB_MO_CB_val_2.fq Sc_CMD3_S8_R2_001_SB_MO_CB_val_2.fq > Sc_CMD_R2_001_SB_MO_CB_val_2.fq
rm Sc_CMD1_S4_R1_001_SB_MO_CB_val_1.fq Sc_CMD2_S6_R1_001_SB_MO_CB_val_1.fq Sc_CMD3_S8_R1_001_SB_MO_CB_val_1.fq
rm Sc_CMD1_S4_R2_001_SB_MO_CB_val_2.fq Sc_CMD2_S6_R2_001_SB_MO_CB_val_2.fq Sc_CMD3_S8_R2_001_SB_MO_CB_val_2.fq

cat Sc_CMr1_S5_R1_001_SB_UM_CB_val_1.fq Sc_CMr2_S7_R1_001_SB_UM_CB_val_1.fq Sc_CMr3_S9_R1_001_SB_UM_CB_val_1.fq > Sc_CMr_R1_001_SB_UM_CB_val_1.fq
cat Sc_CMr1_S5_R2_001_SB_UM_CB_val_2.fq Sc_CMr2_S7_R2_001_SB_UM_CB_val_2.fq Sc_CMr3_S9_R2_001_SB_UM_CB_val_2.fq > Sc_CMr_R2_001_SB_UM_CB_val_2.fq
rm Sc_CMr1_S5_R1_001_SB_UM_CB_val_1.fq Sc_CMr2_S7_R1_001_SB_UM_CB_val_1.fq Sc_CMr3_S9_R1_001_SB_UM_CB_val_1.fq
rm Sc_CMr1_S5_R2_001_SB_UM_CB_val_2.fq Sc_CMr2_S7_R2_001_SB_UM_CB_val_2.fq Sc_CMr3_S9_R2_001_SB_UM_CB_val_2.fq

cat Tag_Records_Sc_CMD1_S4.tsv Tag_Records_Sc_CMD2_S6.tsv Tag_Records_Sc_CMD3_S8.tsv > Tag_Records_Sc_CMD.tsv
cat Tag_Records_Sc_CMr1_S5.tsv Tag_Records_Sc_CMr2_S7.tsv Tag_Records_Sc_CMr3_S9.tsv > Tag_Records_Sc_CMr.tsv
rm Tag_Records_Sc_CMD1_S4.tsv Tag_Records_Sc_CMD2_S6.tsv Tag_Records_Sc_CMD3_S8.tsv
rm Tag_Records_Sc_CMr1_S5.tsv Tag_Records_Sc_CMr2_S7.tsv Tag_Records_Sc_CMr3_S9.tsv

#Cell Line DNA
DName=Sc_CMD
echo "Splitting ${DName}..."
codon run -plugin seq -release /home/annan/triseq052025_2/Split_Reads.codon \
    ${DName} \
    ${outdir} \
    ${libName} \
    exp2d ${outdir}/${DName}_R1_001_SB_MO_CB_val_1.fq ${outdir}/${DName}_R2_001_SB_MO_CB_val_2.fq
rm ${DName}_R1_001_SB_MO_CB_val_1.fq ${DName}_R2_001_SB_MO_CB_val_2.fq

#Cell Line RNA
RName=Sc_CMr
echo "Splitting ${RName}..."
codon run -plugin seq -release /home/annan/triseq052025_2/Split_Reads.codon \
    ${RName} \
    ${outdir} \
    ${libName} \
    exp2r ${outdir}/${RName}_R1_001_SB_UM_CB_val_1.fq ${outdir}/${RName}_R2_001_SB_UM_CB_val_2.fq
rm ${RName}_R1_001_SB_UM_CB_val_1.fq ${RName}_R2_001_SB_UM_CB_val_2.fq


# 5 Alignment and Duplicate detection

grcm39_eff_size=2654621783
hg38_eff_size=2913022398

gtfMM="${refdir}/MM39_Annot_From10x.gtf"
gtfHG="${refdir}/HG38_Annot_From10x.gtf"

bwaHG="${refdir}/hg38_bwamem2_index/hg38.fa"
bwaMM="${refdir}/GRCm39_bwamem2_index/GRCm39"
mods=( "H3K27ac" "H3K27me3" )

DNA_CellLine_split_HG38=( "Sc_CMD_A549" "Sc_CMD_GM12878" "Sc_CMD_HCC827" "Sc_CMD_JJN2" "Sc_CMD_Karpas422" "Sc_CMD_WSU")
RNA_CellLine_split_HG38=( "Sc_CMr_A549" "Sc_CMr_GM12878" "Sc_CMr_HCC827" "Sc_CMr_JJN2" "Sc_CMr_Karpas422" "Sc_CMr_WSU")

DNA_CellLine_split_MM39=( "Sc_CMD_A20" "Sc_CMD_KPAR" )
RNA_CellLine_split_MM39=( "Sc_CMr_A20" "Sc_CMr_KPAR" )

# Check reference files
required_refs=( "${gtfMM}" "${gtfHG}" "${bwaHG}.bwt.2bit.64" "${bwaMM}.bwt.2bit.64" "${refdir}/hg38-blacklist.v2.bed" "${refdir}/refdata-gex-GRCh38-2024-A/star" )
for ref in "${required_refs[@]}"; do
    if [ ! -e "${ref}" ]; then
        echo "ERROR: Required reference file not found: ${ref}"
        exit 1
    fi
done


#Cell Line Human
## DNA
for DName in "${DNA_CellLine_split_HG38[@]}"; do
    for mod in "${mods[@]}"; do
        echo "Aligning ${DName}_${mod}..."
        /home/annan/H2R/AlignDNA.sh ${mod} ${DName} \
            ${outdir}/${DName}_${mod}_R1.fq \
            ${outdir}/${DName}_${mod}_R2.fq \
            ${refdir}/hg38-blacklist.v2.bed \
            ${outdir}/SAM_RG_Header_${DName}_${mod}.tsv \
            ${bwaHG} ${hg38_eff_size} ${outdir}
        
        nreads=$(samtools view -c ${outdir}/${DName}_${mod}.bam)
        if [ "$nreads" -gt 10000 ]; then
            /mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
                -I ${outdir}/${DName}_${mod}.bam -O ${outdir}/${DName}_${mod}_MarkedDup.bam \
                -M ${outdir}/${DName}_${mod}.DuplicateMetrics.txt \
                --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000
            
            samtools view --threads 64 --bam --with-header --require-flags 0x400 --output ${outdir}/DUP.bam --unoutput ${outdir}/${DName}_${mod}_NoDup.bam ${outdir}/${DName}_${mod}_MarkedDup.bam
            samtools index --threads 64 --bai --output ${outdir}/${DName}_${mod}_NoDup.bam.bai ${outdir}/${DName}_${mod}_NoDup.bam
            rm -f ${outdir}/DUP.bam ${outdir}/${DName}_${mod}.bam

            bamCoverage -p 64 -bs 100 --extendReads --centerReads -b ${outdir}/${DName}_${mod}_NoDup.bam -o ${outdir}/${DName}_${mod}_NoDup.bw -of bigwig --effectiveGenomeSize ${hg38_eff_size}
        else
            echo "Skipping duplicate marking and coverage for ${DName}_${mod} - insufficient reads (${nreads} < 10000)"
            mv ${outdir}/${DName}_${mod}.bam ${outdir}/${DName}_${mod}_FAILED.bam
        fi
        echo "Finished aligning ${DName}"
    done
done



#Cell Line Mouse
## DNA
for DName in "${DNA_CellLine_split_MM39[@]}"; do
    for mod in "${mods[@]}"; do
        echo "Aligning ${DName}_${mod}..."
        /home/annan/H2R/AlignDNA.sh ${mod} ${DName} \
            ${outdir}/${DName}_${mod}_R1.fq \
            ${outdir}/${DName}_${mod}_R2.fq \
            ${refdir}/mm10_grcm39LiftOver-blacklist.v2.bed \
            ${outdir}/SAM_RG_Header_${DName}_${mod}.tsv \
            ${bwaMM} ${grcm39_eff_size} ${outdir}
        
        nreads=$(samtools view -c ${outdir}/${DName}_${mod}.bam)
        if [ "$nreads" -gt 10000 ]; then
            /mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk MarkDuplicates \
                -I ${outdir}/${DName}_${mod}.bam -O ${outdir}/${DName}_${mod}_MarkedDup.bam \
                -M ${outdir}/${DName}_${mod}.DuplicateMetrics.txt \
                --REMOVE_DUPLICATES false --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000
            
            samtools view --threads 64 --bam --with-header --require-flags 0x400 --output ${outdir}/DUP.bam --unoutput ${outdir}/${DName}_${mod}_NoDup.bam ${outdir}/${DName}_${mod}_MarkedDup.bam
            samtools index --threads 64 --bai --output ${outdir}/${DName}_${mod}_NoDup.bam.bai ${outdir}/${DName}_${mod}_NoDup.bam
            rm -f ${outdir}/DUP.bam ${outdir}/${DName}_${mod}.bam

            bamCoverage -p 64 -bs 100 --extendReads --centerReads -b ${outdir}/${DName}_${mod}_NoDup.bam -o ${outdir}/${DName}_${mod}_NoDup.bw -of bigwig --effectiveGenomeSize ${grcm39_eff_size}
        else
            echo "Skipping duplicate marking and coverage for ${DName}_${mod} - insufficient reads (${nreads} < 10000)"
            mv ${outdir}/${DName}_${mod}.bam ${outdir}/${DName}_${mod}_FAILED.bam
        fi
        echo "Finished aligning ${DName}"
    done
done

#!/bin/bash

set -euo pipefail

# Globally used variables
indir="."
outdir="/mnt/dataFast/ahrmad/triseq_2025052/RNA_RE"
refdir="/mnt/dataFast/ahrmad"

# 5 Alignment and Duplicate detection
RNAsamples=("GM12878" "JJN2" "Karpas422" "WSU" "A20")
# 
for RName in "${RNAsamples[@]}"; do
    echo "Aligning ${RName}..."
    /home/annan/H2R/AlignRNA.sh ${RName} \
        ${indir}/Sc_CMr_${RName}_R1.fq.gz ${indir}/Sc_CMr_${RName}_R2.fq.gz \
        ${refdir} ${outdir} 64 29
    
    echo "Finished aligning ${RName}"

done


samtools merge -@ 24 -o Sc_CMD_H3K27ac_MarkedDup.tmp.bam Sc_CMD_*H3K27ac*_MarkedDup.bam
samtools sort -@ 64 -o Sc_CMD_H3K27ac_MarkedDup.bam Sc_CMD_H3K27ac_MarkedDup.tmp.bam
rm Sc_CMD_H3K27ac_MarkedDup.tmp.bam
samtools index -@ 48 Sc_CMD_H3K27ac_MarkedDup.bam

samtools merge -@ 24 -o Sc_CMD_H3K27me3_MarkedDup.tmp.bam Sc_CMD_*H3K27me3*_MarkedDup.bam
samtools sort -@ 64 -o Sc_CMD_H3K27me3_MarkedDup.bam Sc_CMD_H3K27me3_MarkedDup.tmp.bam
rm Sc_CMD_H3K27me3_MarkedDup.tmp.bam
samtools index -@ 48 Sc_CMD_H3K27me3_MarkedDup.bam

samtools merge -@ 24 -o Sc_CMr_ALL.tmp.bam Sc_CMr_*ALL.bam
samtools sort -@ 64 -o Sc_CMr_ALL.bam Sc_CMr_ALL.tmp.bam
rm Sc_CMr_ALL.tmp.bam
samtools index -@ 48 Sc_CMr_ALL.bam

samtools merge -@ 24 -o Sc_CMr_UMI.tmp.bam Sc_CMr_*UMI.bam
samtools sort -@ 64 -o Sc_CMr_UMI.bam Sc_CMr_UMI.tmp.bam
rm Sc_CMr_UMI.tmp.bam
samtools index -@ 48 Sc_CMr_UMI.bam

cat *CMD*_H3K27ac_ProperPairedMapped_reads_per_barcode.tsv > Sc_CMD_H3K27ac_ProperPairedMapped_reads_per_barcode.tsv
cat *CMD*_H3K27me3_ProperPairedMapped_reads_per_barcode.tsv > Sc_CMD_H3K27me3_ProperPairedMapped_reads_per_barcode.tsv
cat *CMr*_ProperPairedMapped_reads_per_barcode.tsv > Sc_CMr_ProperPairedMapped_reads_per_barcode.tsv

pigz -p 24 *fq Tag*tsv