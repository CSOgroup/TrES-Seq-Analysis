#!/usr/bin/env bash
set -euo pipefail

# Run after the marker BEDs, filtered BAMs, and peaks are available.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
export TOTAL_CORES="${TOTAL_CORES:-80}"
export BAMC_NORM=RPGC BAMC_CENTER=1 BAMC_EXTEND=1 BAMC_MINMQ=30 BAMC_IGNDUPS=1

cd -- "$WORK_DIR"
export WORK_DIR="$PWD"

# TSS coverage and plots.
BW_DIR=v3_tss_s100_bs80 BAMC_SMOOTH=100 BAMC_BS=80 \
  bash "$SCRIPT_DIR/SupplFig3_4a_MakeBigWigs.sh"
BW_DIR=v3_tss_s100_bs80 MAT_DIR=matrices_tss_s100_bs80 FIG_DIR=figs_tss_s100_bs80 \
  YMAX_TSS_PROFILE=40 YMAX_TSS_HEAT=40 TSS_BIN=80 topn=50 BW_SUFFIX="$BAMC_NORM" \
  bash "$SCRIPT_DIR/SupplFig3_4b_DeepToolsTSS.sh"

# Peak-centered coverage and plots.
MODE=sum_s0_bs10
BW_DIR="v3_${MODE}" BAMC_SMOOTH=0 BAMC_BS=10 \
  bash "$SCRIPT_DIR/SupplFig3_4a_MakeBigWigs.sh"
export YMAX_SUMMITS_PROFILE_GM12878=40 YMAX_SUMMITS_HEAT_GM12878=25
export YMAX_SUMMITS_PROFILE_WSU=20 YMAX_SUMMITS_HEAT_WSU=10
export YMAX_SUMMITS_PROFILE_JJN3=20 YMAX_SUMMITS_HEAT_JJN3=20
export YMAX_SUMMITS_PROFILE_Karpas422=20 YMAX_SUMMITS_HEAT_Karpas422=10
export YMAX_SUMMITS_PROFILE_A20=8 YMAX_SUMMITS_HEAT_A20=10
BW_DIR="v3_${MODE}" MAT_DIR="matrices_${MODE}" PROF_DIR="profiles_${MODE}" HEAT_DIR="heatmaps_${MODE}" \
  BIN=10 BW_SUFFIX="$BAMC_NORM" bash "$SCRIPT_DIR/SupplFig3_4c_DeepToolsPEAKS.sh"
