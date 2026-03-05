#!/usr/bin/env bash
set -euo pipefail

# =======================
# Config (override via env)
# =======================
NPROC="${NPROC:-88}"
FLANK_SUMMITS="${FLANK_SUMMITS:-5000}"         # ± bp around peak centers
BIN="${BIN:-20}"                              

# Global plot caps (per-cell overrides can replace these)
YMAX_SUMMITS_PROFILE="${YMAX_SUMMITS_PROFILE:-60}"   # y-max for profiles
YMAX_SUMMITS_HEAT="${YMAX_SUMMITS_HEAT:-40}"         # z-max for heatmaps

# I/O
BW_DIR="${BW_DIR:-v3}"
MAT_DIR="${MAT_DIR:-matrices_summits}"
PROF_DIR="${PROF_DIR:-profiles}"
HEAT_DIR="${HEAT_DIR:-heatmaps}"
BW_SUFFIX="${BW_SUFFIX:-RPGC}"                       # expects *_Filt_${BW_SUFFIX}.bw

mkdir -p "$MAT_DIR" "$PROF_DIR" "$HEAT_DIR"

# Colors
COL_AC="${COL_AC:-#D55E00}"   # H3K27ac line
COL_ME3="${COL_ME3:-#5E3D99}" # H3K27me3 line

# Tools check
for exe in computeMatrix plotProfile plotHeatmap parallel; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: '$exe' not found in PATH" >&2; exit 1; }
done

# =======================
# Per-cell cap helpers (env override)
# Set env like: YMAX_SUMMITS_PROFILE_GM12878=55 or YMAX_SUMMITS_HEAT_A20=35
# =======================
get_profile_cap() {
  local cell="$1"
  local vname="YMAX_SUMMITS_PROFILE_${cell}"
  local val="${!vname-}"
  if [[ -n "${val:-}" ]]; then echo "$val"; else echo "$YMAX_SUMMITS_PROFILE"; fi
}
get_heat_cap() {
  local cell="$1"
  local vname="YMAX_SUMMITS_HEAT_${cell}"
  local val="${!vname-}"
  if [[ -n "${val:-}" ]]; then echo "$val"; else echo "$YMAX_SUMMITS_HEAT"; fi
}

# =======================
# Discover eligible cells for Peak Centers block
# =======================
CELLS_ALL=(GM12878 JJN2 Karpas422 WSU A20)
CELLS_SUMMITS=()
for C in "${CELLS_ALL[@]}"; do
  ac_bw="${BW_DIR}/ac_${C}_H3K27ac_Filt_${BW_SUFFIX}.bw"
  me3_bw="${BW_DIR}/me3_${C}_H3K27me3_Filt_${BW_SUFFIX}.bw"
  ac_summits="./macs3_out/ac_${C}_H3K27ac_Filt_peaks.narrowPeak.blacklistFiltered"
  me3_summits="./macs3_out/me3_${C}_H3K27me3_Filt_peaks.broadPeak.blacklistFiltered"
  if [[ -s "$ac_bw" && -s "$me3_bw" && -s "$ac_summits" && -s "$me3_summits" ]]; then
    CELLS_SUMMITS+=("$C")
  else
    echo "Summits: skipping ${C} (missing one of: $ac_bw $me3_bw $ac_summits $me3_summits)" >&2 || true
  fi
done
(( ${#CELLS_SUMMITS[@]} > 0 )) || { echo "No cells available for Peak Centers analyses; exiting." >&2; exit 0; }

# Thread planning
CORES_PER_CELL=$(( NPROC / ${#CELLS_SUMMITS[@]} ))
(( CORES_PER_CELL < 1 )) && CORES_PER_CELL=1

# =======================
# Per-cell work
# =======================
summit_block_per_cell() {
  local C="$1"
  local THREADS="$2"
  local ac_bw="${BW_DIR}/ac_${C}_H3K27ac_Filt_${BW_SUFFIX}.bw"
  local me3_bw="${BW_DIR}/me3_${C}_H3K27me3_Filt_${BW_SUFFIX}.bw"
  local ac_summits="./macs3_out/ac_${C}_H3K27ac_Filt_peaks.narrowPeak.blacklistFiltered"
  local me3_summits="./macs3_out/me3_${C}_H3K27me3_Filt_peaks.broadPeak.blacklistFiltered"

  # Per-cell caps (env overrides fallback to global)
  local YP; YP="$(get_profile_cap "$C")"
  local YH; YH="$(get_heat_cap "$C")"

  echo "== ${C} (threads/cell=${THREADS})  caps: profile=${YP}  heat=${YH} =="

  # ----- Profiles over AC Peak Centers (AC + ME3 overlaid) -----
  local mat_ac="${MAT_DIR}/${C}.acSummits.ac+me3.mat.gz"
  computeMatrix reference-point \
    --referencePoint center \
    -R "$ac_summits" \
    -S "$ac_bw" "$me3_bw" \
    --sortRegions descend \
    --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" \
    --averageTypeBins mean \
    --binSize "$BIN" \
    --missingDataAsZero \
    --samplesLabel "H3K27ac" "H3K27me3" \
    --numberOfProcessors "$THREADS" \
    --outFileName "$mat_ac"

  plotProfile -m "$mat_ac" \
    -o "${PROF_DIR}/${C}.profile_over_AC_PeakCenter.pdf" \
    --plotFileFormat pdf \
    --plotType lines \
    --perGroup \
    --plotTitle "${C}: ±$((FLANK_SUMMITS/1000)) kb around H3K27ac Peak Centers" \
    --refPointLabel "Peak" \
    --yAxisLabel "Mean RPGC" \
    --colors "$COL_AC" "$COL_ME3" \
    --legendLocation upper-right \
    --yMin 0 --yMax "$YP"

  # ----- Profiles over ME3 Peak Centers (AC + ME3 overlaid) -----
  local mat_me3="${MAT_DIR}/${C}.me3Summits.ac+me3.mat.gz"
  computeMatrix reference-point \
    --referencePoint center \
    -R "$me3_summits" \
    -S "$ac_bw" "$me3_bw" \
    --sortRegions descend \
    --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" \
    --binSize "$BIN" \
    --averageTypeBins mean \
    --missingDataAsZero \
    --samplesLabel "H3K27ac" "H3K27me3" \
    --numberOfProcessors "$THREADS" \
    --outFileName "$mat_me3"

  plotProfile -m "$mat_me3" \
    -o "${PROF_DIR}/${C}.profile_over_ME3_PeakCenter.pdf" \
    --plotFileFormat pdf \
    --plotType lines \
    --perGroup \
    --plotTitle "${C}: ±$((FLANK_SUMMITS/1000)) kb around H3K27me3 Peak Centers" \
    --refPointLabel "Peak" \
    --yAxisLabel "Mean RPGC" \
    --colors "$COL_AC" "$COL_ME3" \
    --legendLocation upper-right \
    --yMin 0 --yMax "$YP"

  # ----- Heatmaps (single-track each) -----
  # 1) AC signal over AC Peak Centers
  local mat_ac_over_ac="${MAT_DIR}/${C}.acSummits.acOnly.mat.gz"
  computeMatrix reference-point \
    --referencePoint center \
    --sortRegions descend \
    -R "$ac_summits" \
    -S "$ac_bw" \
    --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" \
    --binSize "$BIN" \
    --missingDataAsZero \
    -p "$THREADS" \
    -o "$mat_ac_over_ac"

  plotHeatmap -m "$mat_ac_over_ac" \
    -o "${HEAT_DIR}/${C}.heatmap_ACsignal_over_ACPeaks.pdf" \
    --plotFileFormat pdf \
    --whatToShow "heatmap and colorbar" \
    --plotTitle "${C} • H3K27ac signal over H3K27ac Peak Centers (±$((FLANK_SUMMITS/1000)) kb)" \
    --refPointLabel "Peak" \
    --regionsLabel "H3K27ac regions" \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges \
    --sortRegions descend --sortUsing mean \
    --zMin 0 --zMax "$YH"

  # 2) ME3 signal over AC Peak Centers
  local mat_ac_over_me3sig="${MAT_DIR}/${C}.acSummits.me3Only.mat.gz"
  computeMatrix reference-point \
    --referencePoint center \
    --sortRegions descend \
    -R "$ac_summits" \
    -S "$me3_bw" \
    --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" \
    --binSize "$BIN" \
    --missingDataAsZero \
    -p "$THREADS" \
    -o "$mat_ac_over_me3sig"

  plotHeatmap -m "$mat_ac_over_me3sig" \
    -o "${HEAT_DIR}/${C}.heatmap_ME3signal_over_ACPeaks.pdf" \
    --plotFileFormat pdf \
    --whatToShow "heatmap and colorbar" \
    --plotTitle "${C} • H3K27me3 signal over H3K27ac Peak Centers (±$((FLANK_SUMMITS/1000)) kb)" \
    --refPointLabel "Peak" \
    --regionsLabel "H3K27ac regions" \
    --xAxisLabel "Distance (bp)" \
    --colorMap Purples \
    --sortRegions descend --sortUsing mean \
    --zMin 0 --zMax "$YH"

  # 3) ME3 signal over ME3 Peak Centers
  local mat_me3_over_me3="${MAT_DIR}/${C}.me3Summits.me3Only.mat.gz"
  computeMatrix reference-point \
    --referencePoint center \
    --sortRegions descend \
    -R "$me3_summits" \
    -S "$me3_bw" \
    --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" \
    --binSize "$BIN" \
    --missingDataAsZero \
    -p "$THREADS" \
    -o "$mat_me3_over_me3"

  plotHeatmap -m "$mat_me3_over_me3" \
    -o "${HEAT_DIR}/${C}.heatmap_ME3signal_over_ME3Peaks.pdf" \
    --plotFileFormat pdf \
    --whatToShow "heatmap and colorbar" \
    --plotTitle "${C} • H3K27me3 signal over H3K27me3 Peak Centers (±$((FLANK_SUMMITS/1000)) kb)" \
    --refPointLabel "Peak" \
    --regionsLabel "H3K27me3 regions" \
    --xAxisLabel "Distance (bp)" \
    --colorMap Purples \
    --sortRegions descend --sortUsing mean \
    --zMin 0 --zMax "$YH"

  # 4) AC signal over ME3 Peak Centers
  local mat_me3_over_acsig="${MAT_DIR}/${C}.me3Summits.acOnly.mat.gz"
  computeMatrix reference-point \
    --referencePoint center \
    --sortRegions descend \
    -R "$me3_summits" \
    -S "$ac_bw" \
    --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" \
    --binSize "$BIN" \
    --missingDataAsZero \
    -p "$THREADS" \
    -o "$mat_me3_over_acsig"

  plotHeatmap -m "$mat_me3_over_acsig" \
    -o "${HEAT_DIR}/${C}.heatmap_ACsignal_over_ME3Peaks.pdf" \
    --plotFileFormat pdf \
    --whatToShow "heatmap and colorbar" \
    --plotTitle "${C} • H3K27ac signal over H3K27me3 Peak Centers (±$((FLANK_SUMMITS/1000)) kb)" \
    --refPointLabel "Peak" \
    --regionsLabel "H3K27me3 regions" \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges \
    --sortRegions descend --sortUsing mean \
    --zMin 0 --zMax "$YH"
}

export -f summit_block_per_cell get_profile_cap get_heat_cap
export BW_DIR MAT_DIR PROF_DIR HEAT_DIR COL_AC COL_ME3 FLANK_SUMMITS BIN BW_SUFFIX \
       YMAX_SUMMITS_PROFILE YMAX_SUMMITS_HEAT

echo "== Script 2 (summits) across: ${CELLS_SUMMITS[*]} =="
parallel -j "${#CELLS_SUMMITS[@]}" --line-buffer \
  'summit_block_per_cell {1} "'"$CORES_PER_CELL"'"' ::: "${CELLS_SUMMITS[@]}"

echo "✅ All jobs completed."

