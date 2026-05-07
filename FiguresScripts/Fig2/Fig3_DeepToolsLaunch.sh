#!/usr/bin/env bash
set -euo pipefail

# ===== Shared bamCoverage settings =====
export TOTAL_CORES="${TOTAL_CORES:-80}"
export BAMC_NORM=RPGC      
export BAMC_CENTER=1
export BAMC_EXTEND=1
export BAMC_MINMQ=30
export BAMC_IGNDUPS=1

# ==============================
# TSS block
# ==============================
(
  export BW_DIR="v3_tss_s3000_bs80"
  export BAMC_SMOOTH=3000
  export BAMC_BS=80
  ./Fig3_DeepTools_BwMaking.sh

  BW_DIR="${BW_DIR}" \
  MAT_DIR="matrices_tss_s3000_bs80" \
  FIG_DIR="figs_tss_s3000_bs80" \
  YMAX_TSS_PROFILE=40 \
  YMAX_TSS_HEAT=40 \
  TSS_BIN=80 \
  topn=50 \
  BW_SUFFIX="${BAMC_NORM}" \
  ./Fig3_DeepToolsTSS.sh
)

# ============================================
# Summits block
# ============================================
(
  MODE=sum_s0_bs10
  export BW_DIR="v3_${MODE}"
  export BAMC_SMOOTH=0
  export BAMC_BS=10
  ./Fig3_DeepTools_BwMaking.sh

  # Per-cell caps
  export YMAX_SUMMITS_PROFILE_GM12878=40
  export YMAX_SUMMITS_HEAT_GM12878=25
  export YMAX_SUMMITS_PROFILE_WSU=20
  export YMAX_SUMMITS_HEAT_WSU=10
  export YMAX_SUMMITS_PROFILE_JJN2=20
  export YMAX_SUMMITS_HEAT_JJN2=20
  export YMAX_SUMMITS_PROFILE_Karpas422=20
  export YMAX_SUMMITS_HEAT_Karpas422=10
  export YMAX_SUMMITS_PROFILE_A20=8
  export YMAX_SUMMITS_HEAT_A20=10

  BW_DIR="${BW_DIR}" \
  MAT_DIR="matrices_${MODE}" \
  PROF_DIR="profiles_${MODE}" \
  HEAT_DIR="heatmaps_${MODE}" \
  BIN=10 \
  BW_SUFFIX="${BAMC_NORM}" \
  ./Fig3_DeepToolsPEAKS.sh
)
