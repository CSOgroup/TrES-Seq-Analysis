#!/bin/bash

# Make sure hg19-blacklist.bed is in the assets dir (from https://github.com/nf-core/cutandrun/blob/master/assets/blacklists/hg19-blacklist.bed)

# ./AlignDNA.sh Human H3K27ac H2RD1 H2RD1_Human_H3K27ac_R1.fastq H2RD1_Human_H3K27ac_R2.fastq hg38-blacklist.bed SAM_RG_Header_H2RD1_Human_H3K27ac.tsv hg38_bwamem2_index 2913022398

set -e

threads=72
view_threads=12
sort_threads=64
min_good_reads_in_cells=100

# Input Variables
modality="${1}"
sample_name="${2}"
R1="${3}"
R2="${4}"
blacklist_bed="${5}"
PathSam_header="${6}"
path_bwarefDB="${7}"
effsize="${8}"
outdir="${9}"

PathOutputBam=${outdir}/${sample_name}_${modality}.bam

RGID="${sample_name}_${modality}"

# Alignement
echo "Aligning"
bwa-mem2 mem -t ${threads} -C -o ${RGID}_TEMP.sam ${path_bwarefDB} ${R1} ${R2}

echo "Changing header and sorting alignment..."
# Prepare the header
samtools view -H ${RGID}_TEMP.sam > ${RGID}_TEMPHEADER.sam
# Find the line number of the last occurrence of "@SQ"
last_sq_line=$(grep -n "@SQ" "${RGID}_TEMPHEADER.sam" | tail -n1 | cut -d':' -f1)
# Split the header based on the line number of the last "@SQ" occurrence
head -n "$last_sq_line" "${RGID}_TEMPHEADER.sam" > "${RGID}_TEMPHEADER1.sam"
tail -n +"$((last_sq_line + 1))" "${RGID}_TEMPHEADER.sam" > "${RGID}_TEMPHEADER2.sam"

# Append the new header to the RG-replaced SAM file and convert to BAM
{ cat ${RGID}_TEMPHEADER1.sam ${PathSam_header} ${RGID}_TEMPHEADER2.sam; samtools view --threads ${view_threads} ${RGID}_TEMP.sam; } | samtools sort --threads ${sort_threads} -m 2G -n -o ${RGID}_TEMP.bam -

echo "Get numbers of properly paired and mapped reads per barcode..."
# Remove reads overlapping with the blacklist
samtools view --threads ${view_threads} --bam --with-header --output ${RGID}_TEMP_inBLRegions.bam --unoutput ${RGID}_TEMP_outBLRegions.bam -L ${blacklist_bed} ${RGID}_TEMP.bam

# Get the number of Properly paired and mapped reads per barcode
samtools view --require-flags 0x2 ${RGID}_TEMP_outBLRegions.bam | awk '/RG:Z:/{match($0, /RG:Z:([^\t]*)/); cnts[substr($0, RSTART+5, RLENGTH-5)]++} END {for (bc in cnts) print cnts[bc] "\t" bc}' > ${outdir}/${RGID}_ProperPairedMapped_reads_per_barcode.tsv

echo "Filtering and sorting..."
# Get the barcodes with more than min_good_reads_in_cells reads
awk -v threshold="$min_good_reads_in_cells" -F'\t' '$1 >= threshold {print $2}' ${outdir}/${RGID}_ProperPairedMapped_reads_per_barcode.tsv > ${RGID}_GoodBarcodes.txt
# Keep only the Properly paired and mapped reads from the barcodes with less than min_good_reads_in_cells reads
samtools view --threads ${view_threads} --with-header --require-flags 0x2 --read-group-file ${RGID}_GoodBarcodes.txt --output ${RGID}_TEMP_Good.bam ${RGID}_TEMP_outBLRegions.bam

## Sort and output final bam
samtools sort -@ ${threads} -m 2G -o ${PathOutputBam} ${RGID}_TEMP_Good.bam

# Index the BAM
samtools index --threads ${threads} --bai --output ${PathOutputBam}.bai ${PathOutputBam}

#Remove temp files
echo "Removing temp files..."
rm ${RGID}_TEMP.sam ${RGID}_TEMPHEADER.sam ${RGID}_TEMPHEADER1.sam ${RGID}_TEMPHEADER2.sam ${RGID}_TEMP.bam ${RGID}_TEMP_Good.bam ${RGID}_GoodBarcodes.txt ${RGID}_TEMP_inBLRegions.bam ${RGID}_TEMP_outBLRegions.bam
