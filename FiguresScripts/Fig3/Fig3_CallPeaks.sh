#!/usr/bin/env bash
set -euo pipefail

NCORES=72
OUTDIR="macs3_out"              

# mark modes
AC_MODE="narrow"
ME3_MODE="broad" 

# cutoffs
QVAL_AC="1e-3"
QVAL_ME3="1e-3"
BROAD_CUTOFF_AC="0.1"
BROAD_CUTOFF_ME3="0.1"

# duplicates
KEEP_DUP_AC="0"
KEEP_DUP_ME3="0"

# blacklists
BLACKLIST_BED_HS="../../hg38-blacklist.bed"
BLACKLIST_BED_MM="../../mm10_grcm39LiftOver-blacklist.v2.bed"

# species hints (first match wins)
MOUSE_HINTS="A20;mm39;mm;mouse"
HUMAN_HINTS="GM12878;Karpas422;JJN2;WSU;hg38;hs;human"

# input Fragments
FiltSharedFrags=(
  GM12878_H3K27ac_FiltShared.bed.gz
  Karpas422_H3K27ac_FiltShared.bed.gz
  JJN2_H3K27ac_FiltShared.bed.gz
  WSU_H3K27ac_FiltShared.bed.gz
  GM12878_H3K27me3_FiltShared.bed.gz
  Karpas422_H3K27me3_FiltShared.bed.gz
  JJN2_H3K27me3_FiltShared.bed.gz
  WSU_H3K27me3_FiltShared.bed.gz
  A20_H3K27ac_FiltShared.bed.gz
  A20_H3K27me3_FiltShared.bed.gz
)

mkdir -p "$OUTDIR"

# ---------- HELPERS ----------
contains_any_hint() {
  local hay="$1" hints="$2" IFS=';'
  read -r -a arr <<< "$hints"
  for h in "${arr[@]}"; do [[ -n "$h" && "$hay" == *"$h"* ]] && return 0; done
  return 1
}

infer_species() {
  local base; base=$(basename "$1")
  if contains_any_hint "$base" "$MOUSE_HINTS"; then echo mm; elif contains_any_hint "$base" "$HUMAN_HINTS"; then echo hs; else echo hs; fi
}

blacklist_for_species() { [[ "$1" == "mm" ]] && echo "$BLACKLIST_BED_MM" || echo "$BLACKLIST_BED_HS"; }
keepdup_for_mark()      { [[ "$1" == "H3K27ac" ]] && echo "$KEEP_DUP_AC"   || echo "$KEEP_DUP_ME3"; }

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
  if [[ "$mark" == "H3K27ac" ]]; then pick_mode_and_cutoff "$AC_MODE" "$QVAL_AC" "$BROAD_CUTOFF_AC"
  else                                   pick_mode_and_cutoff "$ME3_MODE" "$QVAL_ME3" "$BROAD_CUTOFF_ME3"
  fi
}

filter_blacklist_all_outputs() {
  # expects MACS --name "${sample}"  → files: ${OUTDIR}/${sample}_peaks.{narrowPeak|broadPeak|gappedPeak}
  local sample="$1" bl="$2"
  [[ -z "$bl" || ! -f "$bl" ]] && return 0
  shopt -s nullglob
  for ext in narrowPeak broadPeak gappedPeak; do
    local infile="${OUTDIR}/${sample}_peaks.${ext}"
    [[ -f "$infile" ]] && bedtools intersect -v -a "$infile" -b "$bl" > "${infile}.blacklistFiltered"
  done
}

run_macs3_one() {
  local fragFile="$1" sample mark species genome_flag bl keepdup modeargs
  sample=$(basename "$fragFile" .bed.gz)

  if   [[ "$fragFile" == *"H3K27ac"* ]];  then mark="H3K27ac"
  elif [[ "$fragFile" == *"H3K27me3"* ]]; then mark="H3K27me3"
  else echo "Skipping $fragFile (cannot infer mark)"; return 0; fi

  species=$(infer_species "$fragFile")
  genome_flag=$([[ "$species" == "mm" ]] && echo "-g mm" || echo "-g hs")
  bl=$(blacklist_for_species "$species")
  keepdup=$(keepdup_for_mark "$mark")
  modeargs=$(mode_args_for_mark "$mark")

  echo ">>> ${sample} | ${mark} | species=${species} | ${genome_flag}"
  echo "    MACS3: ${modeargs}  --keep-dup ${keepdup} --nolambda --nomodel -f FRAG"

  # Use --name "${sample}" to avoid _peaks_peaks.* and make filtering consistent
  # shellcheck disable=SC2086
  macs3 callpeak \
    -t "$fragFile" \
    -f FRAG \
    ${genome_flag} \
    --keep-dup "$keepdup" \
    --nolambda \
    --nomodel \
    --name "${sample}" \
    --outdir "$OUTDIR" \
    --seed 99 \
    --verbose 1 \
    ${modeargs}

  filter_blacklist_all_outputs "$sample" "$bl"
}

export -f run_macs3_one contains_any_hint infer_species blacklist_for_species \
          keepdup_for_mark pick_mode_and_cutoff mode_args_for_mark \
          filter_blacklist_all_outputs
export NCORES OUTDIR AC_MODE ME3_MODE QVAL_AC QVAL_ME3 BROAD_CUTOFF_AC BROAD_CUTOFF_ME3 \
       KEEP_DUP_AC KEEP_DUP_ME3 BLACKLIST_BED_HS BLACKLIST_BED_MM MOUSE_HINTS HUMAN_HINTS

# ---------- RUN ----------
printf "%s\n" "${FiltSharedFrags[@]}" | parallel -j "$NCORES" run_macs3_one {}
echo "✅ Done. Outputs (including *.blacklistFiltered) are in: $OUTDIR"
