#!/usr/bin/env bash
set -euo pipefail

# ---------- USER CONFIG ----------
OUTDIR="macs3_TonsilProcessed"

# Genome ("hs" or "mm") – used for MACS -g and for picking blacklist.
GENOME="hs"

# blacklists
BLACKLIST_BED_HS="../../../hg38-blacklist.bed"

# input FRAGSs (paired-end)
FRAGS_AC="./Fragments_Processed/Sc_VTD9_S1_Sc_VTD11_S2_VTHumanTonsil_H3K27ac.bed.gz"
FRAGS_ME3="./Fragments_Processed/Sc_VTD9_S1_Sc_VTD11_S2_VTHumanTonsil_H3K27me3.bed.gz"

# mark modes (as in script 2)
AC_MODE="narrow"   # "narrow" or "broad"
ME3_MODE="broad"  # set to "broad" if you want broad H3K27me3

# cutoffs
QVAL_AC="1e-3"
QVAL_ME3="1e-3"
BROAD_CUTOFF_AC="0.1"
BROAD_CUTOFF_ME3="0.1"

# duplicates
KEEP_DUP_AC="all"
KEEP_DUP_ME3="all"

mkdir -p "$OUTDIR"

# ---------- HELPERS ----------

pick_mode_and_cutoff() {
  local mode="$1" qval="$2" broad="$3"
  if [[ "$mode" == "narrow" ]]; then
    printf -- "-q %s --call-summits" "$qval"
  else
    printf -- "--broad --broad-cutoff %s" "$broad"
  fi
}

mode_args_for_mark() {
  local mark="$1"
  if [[ "$mark" == "H3K27ac" ]]; then
    pick_mode_and_cutoff "$AC_MODE" "$QVAL_AC" "$BROAD_CUTOFF_AC"
  else
    pick_mode_and_cutoff "$ME3_MODE" "$QVAL_ME3" "$BROAD_CUTOFF_ME3"
  fi
}

filter_blacklist_all_outputs() {
  # expects MACS --name "${sample}"  → ${OUTDIR}/${sample}_peaks.{narrowPeak,broadPeak,gappedPeak}
  local sample="$1" bl="$2"
  [[ -z "$bl" || ! -f "$bl" ]] && return 0
  shopt -s nullglob
  for ext in narrowPeak broadPeak gappedPeak; do
    local infile="${OUTDIR}/${sample}_peaks.${ext}"
    if [[ -f "$infile" ]]; then
      echo ">>> Blacklist-filtering ${infile}"
      bedtools intersect -v -a "$infile" -b "$bl" > "${infile}.blacklistFiltered"
    fi
  done
}

run_macs3_frags() {
  local frags="$1" mark="$2"

  local sample genome_flag bl keepdup modeargs
  sample=$(basename "$frags" .bed.gz)

  genome_flag="-g hs"
  bl=$BLACKLIST_BED_HS

  if [[ "$mark" == "H3K27ac" ]]; then
    keepdup="$KEEP_DUP_AC"
  else
    keepdup="$KEEP_DUP_ME3"
  fi

  modeargs=$(mode_args_for_mark "$mark")

  echo ">>> Running MACS3 on ${sample} (${mark}) | genome=${GENOME}"
  echo "    MACS3: ${modeargs}  --keep-dup ${keepdup} --nolambda --nomodel -f FRAG"

  # shellcheck disable=SC2086
  macs3 callpeak \
    -t "$frags" \
    -f FRAG \
    ${genome_flag} \
    --keep-dup "$keepdup" \
    --nolambda \
    --nomodel \
    --name "$sample" \
    --outdir "$OUTDIR" \
    --seed 99 \
    --verbose 1 \
    ${modeargs}

  filter_blacklist_all_outputs "$sample" "$bl"
}

# ---------- RUN ----------
run_macs3_frags "$FRAGS_AC"  "H3K27ac"
run_macs3_frags "$FRAGS_ME3" "H3K27me3"

echo "✅ Done. Outputs (including *.blacklistFiltered) are in: $OUTDIR"
