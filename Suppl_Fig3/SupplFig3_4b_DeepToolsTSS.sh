#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
BW_DIR="${BW_DIR:-v3}"
TSS_BED_DIR="${TSS_BED_DIR:-tss_beds_topmarkers}"
MAT_DIR="${MAT_DIR:-matrices_tss}"
FIG_DIR="${FIG_DIR:-figs}"
BW_SUFFIX="${BW_SUFFIX:-RPGC}"
topn="${topn:-30}"
NPROC="${NPROC:-88}"
FLANK_TSS="${FLANK_TSS:-20000}"
TSS_BIN="${TSS_BIN:-100}"
YMAX_TSS_PROFILE="${YMAX_TSS_PROFILE:-75}"
YMAX_TSS_HEAT="${YMAX_TSS_HEAT:-40}"
COL_AC="${COL_AC:-#D55E00}"
COL_ME3="${COL_ME3:-#5E3D99}"
CELLS=(GM12878 JJN3 Karpas422 WSU)
PALETTE=("#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b" "#e377c2" "#7f7f7f" "#bcbd22" "#17becf")

cd -- "$WORK_DIR"
mkdir -p "$MAT_DIR" "$FIG_DIR"
for tool in computeMatrix plotProfile plotHeatmap; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done

samples=() regions=() labels=() ac_signals=() me3_signals=()
for cell in "${CELLS[@]}"; do
  bed="$TSS_BED_DIR/tss_top${topn}_${cell}.bed"
  ac="$BW_DIR/${cell}_H3K27ac_FiltShared_${BW_SUFFIX}.bw"
  me3="$BW_DIR/${cell}_H3K27me3_FiltShared_${BW_SUFFIX}.bw"
  if [[ -s "$bed" && -s "$ac" && -s "$me3" ]]; then
    samples+=("$cell")
    regions+=("$bed")
    labels+=("Top${topn}_${cell}")
    ac_signals+=("$ac")
    me3_signals+=("$me3")
  else
    echo "Skipping $cell (need $bed, $ac and $me3)" >&2
  fi
done
(( ${#samples[@]} )) || { echo "No eligible cells for TSS." >&2; exit 0; }
colors=("${PALETTE[@]:0:${#samples[@]}}")

make_matrix() {
  local matrix="$1"
  shift
  computeMatrix reference-point -S "$@" -R "${regions[@]}" --sortRegions keep \
    --beforeRegionStartLength "$FLANK_TSS" --afterRegionStartLength "$FLANK_TSS" \
    --binSize "$TSS_BIN" --averageTypeBins mean --skipZeros --numberOfProcessors "$NPROC" \
    --outFileName "$matrix" --samplesLabel "${samples[@]}"
}

plot_heatmap() {
  local short="$1" mark="$2" color="$3"
  plotHeatmap --matrixFile "$MAT_DIR/TSS_$short.mat.gz" \
    --outFileName "$FIG_DIR/Top${topn}_CL_TSS_${short}_heatmap.pdf" \
    --whatToShow "plot, heatmap and colorbar" --plotType lines --averageTypeSummaryPlot mean \
    --legendLocation upper-right --sortRegions keep --perGroup --boxAroundHeatmaps no \
    --colorList "white,$color" --heatmapHeight 16 --dpi 300 \
    --zMin 0 --zMax "$YMAX_TSS_HEAT" --yMin 0 --yMax "$YMAX_TSS_PROFILE" \
    --refPointLabel TSS --samplesLabel "${samples[@]}" --regionsLabel "${labels[@]}" \
    --plotTitle "$mark signal around TSS (Top${topn} markers per cell line)"
}

plot_profile() {
  local short="$1" mark="$2"
  plotProfile --matrixFile "$MAT_DIR/TSS_$short.mat.gz" \
    --outFileName "$FIG_DIR/Top${topn}_CL_TSS_${short}_profile.pdf" \
    --numPlotsPerRow 4 --perGroup --legendLocation upper-right --yAxisLabel "mean RPGC" \
    --yMin 0 --yMax "$YMAX_TSS_PROFILE" --samplesLabel "${samples[@]}" \
    --regionsLabel "${labels[@]}" --colors "${colors[@]}" \
    --plotTitle "$mark: aggregate signal (Top${topn} markers per cell line)"
}

make_matrix "$MAT_DIR/TSS_AC.mat.gz" "${ac_signals[@]}" & ac_pid=$!
make_matrix "$MAT_DIR/TSS_ME3.mat.gz" "${me3_signals[@]}" & me3_pid=$!
wait "$ac_pid"
wait "$me3_pid"

plot_heatmap AC H3K27ac "$COL_AC" & ac_pid=$!
plot_heatmap ME3 H3K27me3 "$COL_ME3" & me3_pid=$!
wait "$ac_pid"
wait "$me3_pid"

plot_profile AC H3K27ac
plot_profile ME3 H3K27me3
