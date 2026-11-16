#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
DATA_ROOT="${DATA_ROOT:-$SCRIPT_DIR/inputs}"
OUTDIR="${OUTDIR:-macs3_out}"
NCORES="${NCORES:-72}"
BLACKLIST_BED_HS="${BLACKLIST_BED_HS:-$DATA_ROOT/hg38-blacklist.bed}"
BLACKLIST_BED_MM="${BLACKLIST_BED_MM:-$DATA_ROOT/mm10_grcm39LiftOver-blacklist.v2.bed}"
AC_MODE="${AC_MODE:-narrow}"
ME3_MODE="${ME3_MODE:-broad}"
QVAL_AC="${QVAL_AC:-1e-3}"
QVAL_ME3="${QVAL_ME3:-1e-3}"
BROAD_CUTOFF_AC="${BROAD_CUTOFF_AC:-0.1}"
BROAD_CUTOFF_ME3="${BROAD_CUTOFF_ME3:-0.1}"
KEEP_DUP_AC="${KEEP_DUP_AC:-0}"
KEEP_DUP_ME3="${KEEP_DUP_ME3:-0}"
FRAGMENTS=(
  GM12878_H3K27ac_FiltShared.bed.gz
  Karpas422_H3K27ac_FiltShared.bed.gz
  JJN3_H3K27ac_FiltShared.bed.gz
  WSU_H3K27ac_FiltShared.bed.gz
  GM12878_H3K27me3_FiltShared.bed.gz
  Karpas422_H3K27me3_FiltShared.bed.gz
  JJN3_H3K27me3_FiltShared.bed.gz
  WSU_H3K27me3_FiltShared.bed.gz
  A20_H3K27ac_FiltShared.bed.gz
  A20_H3K27me3_FiltShared.bed.gz
)

call_peaks() {
  set -e
  local fragments="$1" sample="${1##*/}" mark mode qvalue broad keep_dup
  local genome=hs blacklist="$BLACKLIST_BED_HS" ext peaks
  sample="${sample%.bed.gz}"
  case "$sample" in
    *A20*|*mm*|*mouse*) genome=mm; blacklist="$BLACKLIST_BED_MM" ;;
  esac
  case "$fragments" in
    *H3K27ac*) mark=H3K27ac; mode="$AC_MODE"; qvalue="$QVAL_AC"; broad="$BROAD_CUTOFF_AC"; keep_dup="$KEEP_DUP_AC" ;;
    *H3K27me3*) mark=H3K27me3; mode="$ME3_MODE"; qvalue="$QVAL_ME3"; broad="$BROAD_CUTOFF_ME3"; keep_dup="$KEEP_DUP_ME3" ;;
    *) echo "Skipping $fragments (unknown mark)" >&2; return 0 ;;
  esac

  local -a peak_options
  if [[ "$mode" == narrow ]]; then
    peak_options=(-q "$qvalue" --call-summits)
  else
    peak_options=(--broad --broad-cutoff "$broad")
  fi
  echo "Calling $mark peaks for $sample ($genome)"
  macs3 callpeak -t "$fragments" -f FRAG -g "$genome" --keep-dup "$keep_dup" \
    --nolambda --nomodel --name "$sample" --outdir "$OUTDIR" --seed 99 --verbose 1 "${peak_options[@]}"

  if [[ -f "$blacklist" ]]; then
    for ext in narrowPeak broadPeak gappedPeak; do
      peaks="$OUTDIR/${sample}_peaks.$ext"
      if [[ -f "$peaks" ]]; then
        bedtools intersect -v -a "$peaks" -b "$blacklist" > "$peaks.blacklistFiltered"
      fi
    done
  fi
}

[[ "$BLACKLIST_BED_HS" = /* ]] || BLACKLIST_BED_HS="$PWD/$BLACKLIST_BED_HS"
[[ "$BLACKLIST_BED_MM" = /* ]] || BLACKLIST_BED_MM="$PWD/$BLACKLIST_BED_MM"
cd -- "$WORK_DIR"
mkdir -p "$OUTDIR"
export -f call_peaks
export OUTDIR AC_MODE ME3_MODE QVAL_AC QVAL_ME3 BROAD_CUTOFF_AC BROAD_CUTOFF_ME3
export KEEP_DUP_AC KEEP_DUP_ME3 BLACKLIST_BED_HS BLACKLIST_BED_MM
PARALLEL_SHELL="$BASH" parallel -j "$NCORES" call_peaks ::: "${FRAGMENTS[@]}"
