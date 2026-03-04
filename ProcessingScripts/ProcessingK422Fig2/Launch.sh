#!/bin/bash

set -e

# Get seq data from EPFL GECF facility
#echo "ls AVT0369/Fastq/FilterOnR1" | sftp -q -P 22 sftpgecf@svsftp.epfl.ch
#sftp -P 22 -B 262144 -R 64 -r #sftpgecf@svsftp.epfl.ch:AVT0369/Fastq/FilterOnR1/ .

# Globally used variables
XDir=/mnt/dataFast/ahrmad/triseq_202510
libName=H2R_TRISEQ_202510
outdir="${XDir}/processed"
refdir="/mnt/dataFast/ahrmad"
rawdir="${XDir}/Raw_Reads"

# K422_A20
DNA_K422_DT=( "Sc_K9D1_S1" "Sc_K9D2_S2" "Sc_K9D3_S3" )
RNA_K422_DT=( "Sc_K9r1_S4" "Sc_K9r2_S5" "Sc_K9r3_S6" )


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
wl_mod="${XDir}/WLs/dna_modality_whitelist.txt"
wl_sample5="${XDir}/WLs/sample5_whitelist.txt"
wl_sample5_dna="${XDir}/WLs/sample5_whitelist_dna.txt"


# Loop through each DNA name in the array
for DName in "${DNA_K422_DT[@]}"
do
    echo "Tagging ${DName}..."

    I1D1="${rawdir}/${DName}_I1_001.fastq.gz"
    I2D1="${rawdir}/${DName}_I2_001.fastq.gz"
    R1D1="${rawdir}/${DName}_R1_001.fastq.gz"
    R2D1="${rawdir}/${DName}_R2_001.fastq.gz"

    # Assign base according to DName
    case "$DName" in
        "Sc_K9D1_S1") BASE="A" ;;
        "Sc_K9D2_S2") BASE="T" ;;
        "Sc_K9D3_S3") BASE="C" ;;
        *) echo "Unknown DName: $DName" && continue ;;
    esac

    # First codon run
    codon run -plugin seq -release -D BC_LEN=3 -D BC_START=0 -D HD=0 /home/annan/H2R/Tag.codon \
        "${I1D1}" \
        "${R1D1}" \
        "${R2D1}" \
        "${wl_sample5_dna}" \
        "${DName}" SB "${outdir}" first_pass_withBC_${BASE} fw

    # Second codon run
    codon run -plugin seq -release -D BC_LEN=8 -D BC_START=3 -D HD=1 /home/annan/H2R/Tag.codon \
        "${I1D1}" \
        "${outdir}/${DName}_R1_001_SB.fastq" \
        "${outdir}/${DName}_R2_001_SB.fastq" \
        "${wl_mod}" \
        "${DName}" MO "${outdir}" not_first_pass fw

    rm "${outdir}/${DName}_R1_001_SB.fastq" \
       "${outdir}/${DName}_R2_001_SB.fastq"

    # Third codon run
    codon run -plugin seq -release -D BC_LEN=8 -D HD=1 /home/annan/triseq102025/Tag_L3_DT.codon \
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


for RName in "${RNA_K422_DT[@]}"
do
    echo "Tagging ${RName}..."

    I1R1="${rawdir}/${RName}_I1_001.fastq.gz"
    R1R1="${rawdir}/${RName}_R1_001.fastq.gz"
    R2R1="${rawdir}/${RName}_R2_001.fastq.gz"

    # Assign base according to RName
    case "$RName" in
        "Sc_K9r1_S4") BASE="A" ;;
        "Sc_K9r2_S5") BASE="T" ;;
        "Sc_K9r3_S6") BASE="C" ;;
        *) echo "Unknown RName: $RName" && continue ;;
    esac


    # First codon run
    codon run -plugin seq -release -D BC_LEN=4 -D BC_START=0 -D HD=1 /home/annan/H2R/Tag.codon \
        "${R2R1}" \
        "${R1R1}" \
        "${R2R1}" \
        "${wl_sample5}" \
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

# 4 Read merging/splitting and SAM RG header buidling
# K422_A20
DNA_K422_DT=( "Sc_K9D1_S1" "Sc_K9D2_S2" "Sc_K9D3_S3" )
RNA_K422_DT=( "Sc_K9r1_S4" "Sc_K9r2_S5" "Sc_K9r3_S6" )

# read / records merging
cat ${outdir}/Tag_Records_Sc_K9D1_S1.tsv ${outdir}/Tag_Records_Sc_K9D2_S2.tsv ${outdir}/Tag_Records_Sc_K9D3_S3.tsv > ${outdir}/Tag_Records_Sc_K9D.tsv
rm ${outdir}/Tag_Records_Sc_K9D1_S1.tsv ${outdir}/Tag_Records_Sc_K9D2_S2.tsv ${outdir}/Tag_Records_Sc_K9D3_S3.tsv

cat ${outdir}/Tag_Records_Sc_K9r1_S4.tsv ${outdir}/Tag_Records_Sc_K9r2_S5.tsv ${outdir}/Tag_Records_Sc_K9r3_S6.tsv > ${outdir}/Tag_Records_Sc_K9r.tsv
rm ${outdir}/Tag_Records_Sc_K9r1_S4.tsv ${outdir}/Tag_Records_Sc_K9r2_S5.tsv ${outdir}/Tag_Records_Sc_K9r3_S6.tsv
pigz -p 72 Tag_Rec*

cat ${outdir}/Sc_K9D1_S1_R1_001_SB_MO_CB_val_1.fq ${outdir}/Sc_K9D2_S2_R1_001_SB_MO_CB_val_1.fq ${outdir}/Sc_K9D3_S3_R1_001_SB_MO_CB_val_1.fq > ${outdir}/Sc_K9D_R1_001_SB_MO_CB_val_1.fq
rm ${outdir}/Sc_K9D1_S1_R1_001_SB_MO_CB_val_1.fq ${outdir}/Sc_K9D2_S2_R1_001_SB_MO_CB_val_1.fq ${outdir}/Sc_K9D3_S3_R1_001_SB_MO_CB_val_1.fq
cat ${outdir}/Sc_K9D1_S1_R2_001_SB_MO_CB_val_2.fq ${outdir}/Sc_K9D2_S2_R2_001_SB_MO_CB_val_2.fq ${outdir}/Sc_K9D3_S3_R2_001_SB_MO_CB_val_2.fq > ${outdir}/Sc_K9D_R2_001_SB_MO_CB_val_2.fq
rm ${outdir}/Sc_K9D1_S1_R2_001_SB_MO_CB_val_2.fq ${outdir}/Sc_K9D2_S2_R2_001_SB_MO_CB_val_2.fq ${outdir}/Sc_K9D3_S3_R2_001_SB_MO_CB_val_2.fq

cat ${outdir}/Sc_K9r1_S4_R1_001_SB_UM_CB_val_1.fq ${outdir}/Sc_K9r2_S5_R1_001_SB_UM_CB_val_1.fq ${outdir}/Sc_K9r3_S6_R1_001_SB_UM_CB_val_1.fq > ${outdir}/Sc_K9r_R1_001_SB_UM_CB_val_1.fq
cat ${outdir}/Sc_K9r1_S4_R2_001_SB_UM_CB_val_2.fq ${outdir}/Sc_K9r2_S5_R2_001_SB_UM_CB_val_2.fq ${outdir}/Sc_K9r3_S6_R2_001_SB_UM_CB_val_2.fq > ${outdir}/Sc_K9r_R2_001_SB_UM_CB_val_2.fq
rm ${outdir}/Sc_K9r1_S4_R1_001_SB_UM_CB_val_1.fq ${outdir}/Sc_K9r2_S5_R1_001_SB_UM_CB_val_1.fq ${outdir}/Sc_K9r3_S6_R1_001_SB_UM_CB_val_1.fq
rm ${outdir}/Sc_K9r1_S4_R2_001_SB_UM_CB_val_2.fq ${outdir}/Sc_K9r2_S5_R2_001_SB_UM_CB_val_2.fq ${outdir}/Sc_K9r3_S6_R2_001_SB_UM_CB_val_2.fq

#K422
DName=Sc_K9D
echo "Splitting ${DName}..."
codon run -plugin seq -release /home/annan/triseq102025/Split_Reads.codon \
    ${DName} \
    ${outdir} \
    ${libName} \
    exp2d ${outdir}/${DName}_R1_001_SB_MO_CB_val_1.fq ${outdir}/${DName}_R2_001_SB_MO_CB_val_2.fq
#rm ${outdir}/${DName}_R1_001_SB_MO_CB_val_1.fq ${outdir}/${DName}_R2_001_SB_MO_CB_val_2.fq


RName=Sc_K9r
echo "Splitting ${RName}..."
codon run -plugin seq -release /home/annan/triseq102025/Split_Reads.codon \
    ${RName} \
    ${outdir} \
    ${libName} \
    exp2r ${outdir}/${RName}_R1_001_SB_UM_CB_val_1.fq ${outdir}/${RName}_R2_001_SB_UM_CB_val_2.fq
#rm ${outdir}/${RName}_R1_001_SB_UM_CB_val_1.fq ${outdir}/${RName}_R2_001_SB_UM_CB_val_2.fq


# 5 Alignment and Duplicate detection

grcm39_eff_size=2654621783
hg38_eff_size=2913022398

gtfMM="${refdir}/MM39_Annot_From10x.gtf"
gtfHG="${refdir}/HG38_Annot_From10x.gtf"

bwaHG="${refdir}/hg38_bwamem2_index/hg38.fa"
bwaMM="${refdir}/GRCm39_bwamem2_index/GRCm39"

mods2=( "H3K9me3" "H3K27me3" )

DNA_K422_Split=( "Sc_K9D" )
RNA_K422_Split=( "Sc_K9r" )

### HUMAN

## K422
for DName in "${DNA_K422_Split[@]}"; do
    for mod in "${mods2[@]}"; do
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

for RName in "${RNA_K422_Split[@]}"; do
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
python ./Make_Plots.py


