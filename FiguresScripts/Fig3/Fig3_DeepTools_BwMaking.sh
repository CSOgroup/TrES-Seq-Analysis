#!/usr/bin/env bash
set -euo pipefail

# ------------------------
# Tunables (override via env)
# ------------------------
TOTAL_CORES="${TOTAL_CORES:-80}"

# Output BW directory (override per config)
BW_DIR="${BW_DIR:-v3}"

# bamCoverage knobs (override per config)
BAMC_BS="${BAMC_BS:-20}"            # --binSize
BAMC_SMOOTH="${BAMC_SMOOTH:-100}"  # --smoothLength; set "" or 0 to disable
BAMC_CENTER="${BAMC_CENTER:-1}"     # 1=>--centerReads
BAMC_EXTEND="${BAMC_EXTEND:-1}"     # 1=>--extendReads
BAMC_NORM="${BAMC_NORM:-RPGC}"      # RPGC or CPM
BAMC_MINMQ="${BAMC_MINMQ:-}"        # e.g. 30; empty to skip
BAMC_IGNDUPS="${BAMC_IGNDUPS:-0}"   # 1=>--ignoreDuplicates

# ------------------------
# BAM files (unchanged)
# ------------------------
BAM_FILES=(
  GM12878_H3K27ac_FiltShared.bam
  Karpas422_H3K27ac_FiltShared.bam
  JJN2_H3K27ac_FiltShared.bam
  WSU_H3K27ac_FiltShared.bam
  GM12878_H3K27me3_FiltShared.bam
  Karpas422_H3K27me3_FiltShared.bam
  JJN2_H3K27me3_FiltShared.bam
  WSU_H3K27me3_FiltShared.bam
  A20_H3K27ac_FiltShared.bam
  A20_H3K27me3_FiltShared.bam
)

# Divide cores equally among jobs (≥1 per job)
CORES_PER_JOB=$(( TOTAL_CORES / ${#BAM_FILES[@]} ))
(( CORES_PER_JOB < 1 )) && CORES_PER_JOB=1

# Effective genome sizes
EGS_DEFAULT=2913022398   # human
EGS_A20=2654621783       # A20 (mouse)

export CORES_PER_JOB EGS_DEFAULT EGS_A20 BW_DIR \
       BAMC_BS BAMC_SMOOTH BAMC_CENTER BAMC_EXTEND BAMC_NORM \
       BAMC_MINMQ BAMC_IGNDUPS

mkdir -p "${BW_DIR}"

echo "== aa.sh config =="
echo "BW_DIR=${BW_DIR}"
echo "binSize=${BAMC_BS} smooth=${BAMC_SMOOTH:-none} center=${BAMC_CENTER} extend=${BAMC_EXTEND}"
echo "norm=${BAMC_NORM} minMQ=${BAMC_MINMQ:-none} ignoreDup=${BAMC_IGNDUPS}"
echo "TOTAL_CORES=${TOTAL_CORES} -> per-job=${CORES_PER_JOB}"
echo

# Run bamCoverage in parallel with per-file genome size selection
parallel -j ${#BAM_FILES[@]} '
  bam=$(basename {})
  out_prefix=${BW_DIR}/${bam%.bam}_RPGC.bw

  # Choose effective genome size based on filename
  egs=$EGS_DEFAULT
  [[ "$bam" == *A20* ]] && egs=$EGS_A20

  # Build extra flags
  extra=()
  [[ "${BAMC_EXTEND}" == "1" ]] && extra+=(--extendReads)
  [[ "${BAMC_CENTER}" == "1" ]] && extra+=(--centerReads)
  if [[ -n "${BAMC_SMOOTH}" && "${BAMC_SMOOTH}" != "0" ]]; then
    extra+=(--smoothLength "${BAMC_SMOOTH}")
  fi
  [[ -n "${BAMC_MINMQ}" ]] && extra+=(--minMappingQuality "${BAMC_MINMQ}")
  [[ "${BAMC_IGNDUPS}" == "1" ]] && extra+=(--ignoreDuplicates)

  echo "Processing $bam → $out_prefix using $CORES_PER_JOB threads (EGS=$egs)"
  bamCoverage \
    -p "$CORES_PER_JOB" \
    -bs "${BAMC_BS}" \
    --normalizeUsing "${BAMC_NORM}" \
    -b "$bam" \
    -o "$out_prefix" \
    -of bigwig \
    --effectiveGenomeSize "$egs" \
    "${extra[@]}"
' ::: "${BAM_FILES[@]}"

echo "✅ All bamCoverage jobs completed successfully."
