#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
BW_DIR="${BW_DIR:-v3}"
PEAK_DIR="${PEAK_DIR:-${OUTDIR:-./macs3_out}}"
MAT_DIR="${MAT_DIR:-matrices_summits}"
PROF_DIR="${PROF_DIR:-profiles}"
HEAT_DIR="${HEAT_DIR:-heatmaps}"
BW_SUFFIX="${BW_SUFFIX:-RPGC}"
NPROC="${NPROC:-88}"
FLANK_SUMMITS="${FLANK_SUMMITS:-5000}"
BIN="${BIN:-20}"
YMAX_SUMMITS_PROFILE="${YMAX_SUMMITS_PROFILE:-60}"
YMAX_SUMMITS_HEAT="${YMAX_SUMMITS_HEAT:-40}"
COL_AC="${COL_AC:-#D55E00}"
COL_ME3="${COL_ME3:-#5E3D99}"
CELLS=(GM12878 JJN3 Karpas422 WSU A20)

load_paths() {
  local cell="$1"
  ac_bw="$BW_DIR/${cell}_H3K27ac_FiltShared_${BW_SUFFIX}.bw"
  me3_bw="$BW_DIR/${cell}_H3K27me3_FiltShared_${BW_SUFFIX}.bw"
  ac_peaks="$PEAK_DIR/${cell}_H3K27ac_FiltShared_peaks.narrowPeak.blacklistFiltered"
  me3_peaks="$PEAK_DIR/${cell}_H3K27me3_FiltShared_peaks.broadPeak.blacklistFiltered"
}

plot_over_peaks() {
  local peaks="$1" mark="$2" short="$3" stem="$4"
  local matrix="$MAT_DIR/$cell.$stem.ac+me3.mat.gz"
  computeMatrix reference-point --referencePoint center -R "$peaks" -S "$ac_bw" "$me3_bw" \
    --sortRegions descend --beforeRegionStartLength "$FLANK_SUMMITS" \
    --afterRegionStartLength "$FLANK_SUMMITS" --binSize "$BIN" --averageTypeBins mean \
    --missingDataAsZero --samplesLabel H3K27ac H3K27me3 --numberOfProcessors "$threads" \
    --outFileName "$matrix"
  plotProfile -m "$matrix" -o "$PROF_DIR/$cell.profile_over_${short}_PeakCenter.pdf" \
    --plotFileFormat pdf --plotType lines --perGroup \
    --plotTitle "$cell: ±$(( FLANK_SUMMITS / 1000 )) kb around $mark Peak Centers" \
    --refPointLabel Peak --yAxisLabel "Mean RPGC" --colors "$COL_AC" "$COL_ME3" \
    --legendLocation upper-right --yMin 0 --yMax "$profile_cap"
}

plot_peak_heatmap() {
  local peaks="$1" peak_mark="$2" peak_short="$3" peak_stem="$4"
  local bw="$5" signal_mark="$6" signal_short="$7" signal_stem="$8" cmap="$9"
  local matrix="$MAT_DIR/$cell.$peak_stem.$signal_stem.mat.gz"
  computeMatrix reference-point --referencePoint center --sortRegions descend -R "$peaks" -S "$bw" \
    --beforeRegionStartLength "$FLANK_SUMMITS" --afterRegionStartLength "$FLANK_SUMMITS" \
    --binSize "$BIN" --missingDataAsZero -p "$threads" -o "$matrix"
  plotHeatmap -m "$matrix" -o "$HEAT_DIR/$cell.heatmap_${signal_short}signal_over_${peak_short}Peaks.pdf" \
    --plotFileFormat pdf --whatToShow "heatmap and colorbar" \
    --plotTitle "$cell • $signal_mark signal over $peak_mark Peak Centers (±$(( FLANK_SUMMITS / 1000 )) kb)" \
    --refPointLabel Peak --regionsLabel "$peak_mark regions" --xAxisLabel "Distance (bp)" \
    --colorMap "$cmap" --sortRegions descend --sortUsing mean --zMin 0 --zMax "$heat_cap"
}

peak_plots() {
  set -e
  local cell="$1" threads="$2" ac_bw me3_bw ac_peaks me3_peaks
  local profile_var="YMAX_SUMMITS_PROFILE_$1" heat_var="YMAX_SUMMITS_HEAT_$1"
  local profile_cap="${!profile_var:-$YMAX_SUMMITS_PROFILE}" heat_cap="${!heat_var:-$YMAX_SUMMITS_HEAT}"
  load_paths "$cell"
  plot_over_peaks "$ac_peaks" H3K27ac AC acSummits
  plot_over_peaks "$me3_peaks" H3K27me3 ME3 me3Summits
  plot_peak_heatmap "$ac_peaks" H3K27ac AC acSummits "$ac_bw" H3K27ac AC acOnly Oranges
  plot_peak_heatmap "$ac_peaks" H3K27ac AC acSummits "$me3_bw" H3K27me3 ME3 me3Only Purples
  plot_peak_heatmap "$me3_peaks" H3K27me3 ME3 me3Summits "$me3_bw" H3K27me3 ME3 me3Only Purples
  plot_peak_heatmap "$me3_peaks" H3K27me3 ME3 me3Summits "$ac_bw" H3K27ac AC acOnly Oranges
}

cd -- "$WORK_DIR"
mkdir -p "$MAT_DIR" "$PROF_DIR" "$HEAT_DIR"
for tool in computeMatrix plotProfile plotHeatmap parallel; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
eligible=()
for cell in "${CELLS[@]}"; do
  load_paths "$cell"
  if [[ -s "$ac_bw" && -s "$me3_bw" && -s "$ac_peaks" && -s "$me3_peaks" ]]; then
    eligible+=("$cell")
  else
    echo "Skipping $cell (need $ac_bw, $me3_bw, $ac_peaks and $me3_peaks)" >&2
  fi
done
(( ${#eligible[@]} )) || { echo "No eligible cells for peak-centered plots." >&2; exit 0; }
threads=$(( NPROC / ${#eligible[@]} ))
(( threads >= 1 )) || threads=1
export -f peak_plots load_paths plot_over_peaks plot_peak_heatmap
export BW_DIR PEAK_DIR MAT_DIR PROF_DIR HEAT_DIR COL_AC COL_ME3 FLANK_SUMMITS BIN BW_SUFFIX
export YMAX_SUMMITS_PROFILE YMAX_SUMMITS_HEAT
PARALLEL_SHELL="$BASH" parallel --line-buffer -j "${#eligible[@]}" peak_plots '{1}' "$threads" ::: "${eligible[@]}"
