#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
BW_DIR="${BW_DIR:-v3}"
TOTAL_CORES="${TOTAL_CORES:-80}"
BAMC_BS="${BAMC_BS:-20}"
BAMC_SMOOTH="${BAMC_SMOOTH:-100}"  # Set to 0 to disable smoothing.
BAMC_CENTER="${BAMC_CENTER:-1}"
BAMC_EXTEND="${BAMC_EXTEND:-1}"
BAMC_NORM="${BAMC_NORM:-RPGC}"
BAMC_MINMQ="${BAMC_MINMQ:-}"
BAMC_IGNDUPS="${BAMC_IGNDUPS:-0}"
EGS_DEFAULT=2913022398
EGS_A20=2654621783
BAM_FILES=(
  GM12878_H3K27ac_FiltShared.bam
  Karpas422_H3K27ac_FiltShared.bam
  JJN3_H3K27ac_FiltShared.bam
  WSU_H3K27ac_FiltShared.bam
  GM12878_H3K27me3_FiltShared.bam
  Karpas422_H3K27me3_FiltShared.bam
  JJN3_H3K27me3_FiltShared.bam
  WSU_H3K27me3_FiltShared.bam
  A20_H3K27ac_FiltShared.bam
  A20_H3K27me3_FiltShared.bam
)

coverage() {
  local bam="$1" name="${1##*/}" genome_size="$EGS_DEFAULT"
  local -a flags=()
  [[ "$name" == *A20* ]] && genome_size="$EGS_A20"
  [[ "$BAMC_EXTEND" == 1 ]] && flags+=(--extendReads)
  [[ "$BAMC_CENTER" == 1 ]] && flags+=(--centerReads)
  [[ -n "$BAMC_SMOOTH" && "$BAMC_SMOOTH" != 0 ]] && flags+=(--smoothLength "$BAMC_SMOOTH")
  [[ -n "$BAMC_MINMQ" ]] && flags+=(--minMappingQuality "$BAMC_MINMQ")
  [[ "$BAMC_IGNDUPS" == 1 ]] && flags+=(--ignoreDuplicates)

  bamCoverage -p "$CORES_PER_JOB" -bs "$BAMC_BS" --normalizeUsing "$BAMC_NORM" \
    -b "$bam" -o "$BW_DIR/${name%.bam}_RPGC.bw" -of bigwig \
    --effectiveGenomeSize "$genome_size" "${flags[@]}"
}

cd -- "$WORK_DIR"
mkdir -p "$BW_DIR"
CORES_PER_JOB=$(( TOTAL_CORES / ${#BAM_FILES[@]} ))
(( CORES_PER_JOB >= 1 )) || CORES_PER_JOB=1
export -f coverage
export BW_DIR CORES_PER_JOB EGS_DEFAULT EGS_A20
export BAMC_BS BAMC_SMOOTH BAMC_CENTER BAMC_EXTEND BAMC_NORM BAMC_MINMQ BAMC_IGNDUPS
PARALLEL_SHELL="$BASH" parallel -j "${#BAM_FILES[@]}" coverage ::: "${BAM_FILES[@]}"
