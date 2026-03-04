#!/bin/bash
# Usage:
#   ./AlignRNA.sh SAMPLE_NAME TAGGED.usam REF_BASE_DIR OUTDIR THREADS SPECIES
# SPECIES: human | mouse

set -euo pipefail

if [[ $# -lt 6 ]]; then
    echo "Usage: $0 SAMPLE_NAME TAGGED.usam REF_BASE_DIR OUTDIR THREADS SPECIES(human|mouse)" >&2
    exit 1
fi

# Input Variables
sample_name="${1}"
USAM_IN="${2}"
path_ref="${3}"
outdir="${4}"
threads="${5}"
species="${6}"

# References
path_refDB="${path_ref}/GRCh38_TrES/star"
path_refCHROMSIZES="${path_ref}/hg38.chrom.sizes"

if [[ "${species}" == "mouse" ]]; then
    path_refDB="${path_ref}/GRCm39_TrES/star"
    path_refCHROMSIZES="${path_ref}/mm39.chrom.sizes"
elif [[ "${species}" != "human" ]]; then
    echo "ERROR: SPECIES must be 'human' or 'mouse' (got: ${species})" >&2
    exit 1
fi

if [[ ! -s "${USAM_IN}" ]]; then
    echo "ERROR: Input SAM missing or empty: ${USAM_IN}" >&2
    exit 1
fi

# Fixed UMI length
UMIlen=10

# Auto-detect CBlen from CB tag or from CR tag (fallback: len(CR)-UMIlen)
CBlen="$(
  awk -v UMIlen="${UMIlen}" '
    BEGIN { FS="\t" }
    $0 ~ /^@/ { next }
    {
      cb = ""; cr = "";
      for (i=12; i<=NF; i++) {
        if ($i ~ /^CB:Z:/) { cb = substr($i,6); break }
      }
      if (cb != "") { print length(cb); exit }

      for (i=12; i<=NF; i++) {
        if ($i ~ /^CR:Z:/) { cr = substr($i,6); break }
      }
      if (cr != "") { print length(cr) - UMIlen; exit }
    }
  ' "${USAM_IN}"
)"

if [[ -z "${CBlen}" ]] || [[ "${CBlen}" -le 0 ]]; then
  echo "ERROR: Could not detect CB length from ${USAM_IN} (no CB:Z: or CR:Z: found)." >&2
  exit 1
fi

# UMI starts right after CB
UMIstart=$(( CBlen + 1 ))

echo "Using species=${species}"
echo "Detected CBlen=${CBlen} => UMIstart=${UMIstart} UMIlen=${UMIlen}"

ulimit -n 32000 2>/dev/null || echo "WARNING: using ulimit -n = $(ulimit -n)" >&2

# STAR run: mapping + STARsolo from SAM
STAR \
  --genomeDir "${path_refDB}" \
  --runThreadN "${threads}" \
  --readFilesIn "${USAM_IN}" \
  --readFilesType SAM PE \
  --twopassMode Basic \
  --outFilterType BySJout \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --sjdbScore 2 \
  --limitSjdbInsertNsj 4000000 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --outFilterScoreMinOverLread 0.33 \
  --outFilterMatchNminOverLread 0.33 \
  --outFilterMultimapNmax 20 \
  --outSAMstrandField intronMotif \
  --outFilterIntronMotifs RemoveNoncanonical \
  --soloType CB_UMI_Simple \
  --soloInputSAMattrBarcodeSeq CR \
  --soloInputSAMattrBarcodeQual - \
  --soloCBwhitelist None \
  --soloCBstart 1 --soloCBlen "${CBlen}" \
  --soloUMIstart "${UMIstart}" --soloUMIlen "${UMIlen}" \
  --soloBarcodeReadLength 0 \
  --soloFeatures GeneFull \
  --soloStrand Forward \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR 15000 0.99 30 20000 90000 50 0.01 20000 0.05 10000 \
  --soloCellReadStats Standard \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI nM AS CB UB GX GN NM MD jM jI MC  \
  --outSAMunmapped None \
  --soloOutFileNames "Solo.out" "features.tsv" "barcodes.tsv" "matrix.mtx" \
  --outFileNamePrefix "${outdir}/${sample_name}."

BARCODES="${outdir}/${sample_name}.Solo.outGeneFull/filtered/barcodes.tsv"
INBAM="${outdir}/${sample_name}.Aligned.sortedByCoord.out.bam"
OUTBAM="${outdir}/${sample_name}.filtered_cells.bam"

samtools view --threads ${threads} --with-header --exclude-flags 0x100 --require-flags 0x1,0x2 --tag-file RG:${BARCODES} --output ${OUTBAM} ${INBAM}
samtools index --threads ${threads} ${INBAM}

rm "${outdir}/${sample_name}.Aligned.sortedByCoord.out.bam"*

STAR \
  --runMode inputAlignmentsFromBAM \
  --runThreadN "${threads}" \
  --genomeDir "${path_refDB}" \
  --inputBAMfile "${OUTBAM}" \
  --outWigType bedGraph \
  --outWigStrand Stranded \
  --outWigNorm RPM \
  --outWigReferencesPrefix chr \
  --outFileNamePrefix "${outdir}/${sample_name}.stranded_"

STAR \
  --runMode inputAlignmentsFromBAM \
  --runThreadN "${threads}" \
  --genomeDir "${path_refDB}" \
  --inputBAMfile "${OUTBAM}" \
  --outWigType bedGraph \
  --outWigStrand Unstranded \
  --outWigNorm RPM \
  --outWigReferencesPrefix chr \
  --outFileNamePrefix "${outdir}/${sample_name}.unstranded_"

# Stranded:
for f in "${outdir}/${sample_name}.stranded_Signal.Unique.str"*.bg; do
  sort -k1,1 -k2,2n "$f" > "${f%.bedGraph}.sorted.bg"
  bedGraphToBigWig "${f%.bedGraph}.sorted.bg" "${path_refCHROMSIZES}" "${f%.bg}.bw"
done

# Unstranded:
for f in "${outdir}/${sample_name}.unstranded_Signal.Unique.str"*.bg; do
  sort -k1,1 -k2,2n "$f" > "${f%.bedGraph}.sorted.bg"
  bedGraphToBigWig "${f%.bedGraph}.sorted.bg" "${path_refCHROMSIZES}" "${f%.bg}.bw"
done

# Keep cleanup close to original, but avoid deleting unrelated bg files
rm -f "${outdir}/${sample_name}"*.bg \
      "${outdir}/${sample_name}"*.bg