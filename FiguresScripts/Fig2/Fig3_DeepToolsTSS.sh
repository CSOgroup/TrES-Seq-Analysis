#!/usr/bin/env bash
set -euo pipefail

# =======================
# Config (override via env)
# =======================
topn="${topn:-30}"
NPROC="${NPROC:-88}"
FLANK_TSS="${FLANK_TSS:-20000}"            # ± bp around TSS
TSS_BIN="${TSS_BIN:-100}"                  # binSize for computeMatrix at TSS

# Plot caps (defaults per your latest request)
YMAX_TSS_PROFILE="${YMAX_TSS_PROFILE:-75}" # y-max for summary/profile
YMAX_TSS_HEAT="${YMAX_TSS_HEAT:-40}"       # z-max for heatmap

# I/O
BW_DIR="${BW_DIR:-v3}"
MAT_DIR="${MAT_DIR:-matrices_tss}"
FIG_DIR="${FIG_DIR:-figs}"
TSS_BED_DIR="${TSS_BED_DIR:-tss_beds_topmarkers}"
BW_SUFFIX="${BW_SUFFIX:-RPGC}"             # expects files like *_Filt_${BW_SUFFIX}.bw

mkdir -p "$MAT_DIR" "$FIG_DIR"

# Colors for marks (can override via env)
COL_AC="${COL_AC:-#D55E00}"   # H3K27ac heatmap high color
COL_ME3="${COL_ME3:-#5E3D99}" # H3K27me3 heatmap high color

# Colors (Tab10) — quote hex to avoid '#' as comment
PALETTE=("#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b" "#e377c2" "#7f7f7f" "#bcbd22" "#17becf")

# Tools check
for exe in computeMatrix plotProfile plotHeatmap; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: '$exe' not found in PATH" >&2; exit 1; }
done

# =======================
# Discover eligible cells for TSS block
# =======================
TSS_BASE_CELLS=(GM12878 JJN2 Karpas422 WSU)
CELLS_TSS=() AC_SIGNALS=() ME3_SIGNALS=() REGIONS=() REGION_LABELS=()

maybe_add_tss_cell() {
  local C="$1"
  local bed="${TSS_BED_DIR}/tss_top${topn}_${C}.bed"
  local ac_bw="${BW_DIR}/${C}_H3K27ac_FiltShared_${BW_SUFFIX}.bw"
  local me3_bw="${BW_DIR}/${C}_H3K27me3_FiltShared_${BW_SUFFIX}.bw"
  if [[ -s "$bed" && -s "$ac_bw" && -s "$me3_bw" ]]; then
    CELLS_TSS+=("$C")
    REGIONS+=("$bed")
    REGION_LABELS+=("Top${topn}_${C}")
    AC_SIGNALS+=("$ac_bw")
    ME3_SIGNALS+=("$me3_bw")
  else
    echo "TSS: skipping ${C} (need bed + both BWs). Missing one of: $bed $ac_bw $me3_bw" >&2 || true
  fi
}

for C in "${TSS_BASE_CELLS[@]}"; do maybe_add_tss_cell "$C"; done
(( ${#CELLS_TSS[@]} > 0 )) || { echo "No eligible cells for TSS; exiting." >&2; exit 0; }

# Labels/colors
SAMPLES_LABEL=( "${CELLS_TSS[@]}" )
COLORS_PROFILE=( "${PALETTE[@]:0:${#CELLS_TSS[@]}}" )

# =======================
# Compute matrices (AC & ME3 in background)
# =======================
(
  computeMatrix reference-point \
    -S "${AC_SIGNALS[@]}" \
    -R "${REGIONS[@]}" \
    --sortRegions keep \
    --beforeRegionStartLength "$FLANK_TSS" \
    --afterRegionStartLength "$FLANK_TSS" \
    --binSize "$TSS_BIN" \
    --averageTypeBins mean \
    --skipZeros \
    --numberOfProcessors "$NPROC" \
    --outFileName "$MAT_DIR/TSS_AC.mat.gz" \
    --samplesLabel "${SAMPLES_LABEL[@]}"
) &

(
  computeMatrix reference-point \
    -S "${ME3_SIGNALS[@]}" \
    -R "${REGIONS[@]}" \
    --sortRegions keep \
    --beforeRegionStartLength "$FLANK_TSS" \
    --afterRegionStartLength "$FLANK_TSS" \
    --binSize "$TSS_BIN" \
    --averageTypeBins mean \
    --skipZeros \
    --numberOfProcessors "$NPROC" \
    --outFileName "$MAT_DIR/TSS_ME3.mat.gz" \
    --samplesLabel "${SAMPLES_LABEL[@]}"
) &
wait

# =======================
# Heatmaps + summary plot (AC & ME3)
# =======================
(
  plotHeatmap \
    --matrixFile "$MAT_DIR/TSS_AC.mat.gz" \
    --outFileName "$FIG_DIR/Top${topn}_CL_TSS_AC_heatmap.pdf" \
    --whatToShow "plot, heatmap and colorbar" \
    --plotType lines \
    --averageTypeSummaryPlot mean \
    --legendLocation upper-right \
    --sortRegions keep \
    --perGroup \
    --boxAroundHeatmaps no \
    --colorList "white,${COL_AC}" \
    --heatmapHeight 16 \
    --dpi 300 \
    --zMin 0 --zMax "${YMAX_TSS_HEAT}" \
    --yMin 0 --yMax "${YMAX_TSS_PROFILE}" \
    --refPointLabel "TSS" \
    --samplesLabel "${SAMPLES_LABEL[@]}" \
    --regionsLabel "${REGION_LABELS[@]}" \
    --plotTitle "H3K27ac signal around TSS (Top${topn} markers per cell line)"
) &

(
  plotHeatmap \
    --matrixFile "$MAT_DIR/TSS_ME3.mat.gz" \
    --outFileName "$FIG_DIR/Top${topn}_CL_TSS_ME3_heatmap.pdf" \
    --whatToShow "plot, heatmap and colorbar" \
    --plotType lines \
    --averageTypeSummaryPlot mean \
    --legendLocation upper-right \
    --sortRegions keep \
    --perGroup \
    --boxAroundHeatmaps no \
    --colorList "white,${COL_ME3}" \
    --heatmapHeight 16 \
    --dpi 300 \
    --zMin 0 --zMax "${YMAX_TSS_HEAT}" \
    --yMin 0 --yMax "${YMAX_TSS_PROFILE}" \
    --refPointLabel "TSS" \
    --samplesLabel "${SAMPLES_LABEL[@]}" \
    --regionsLabel "${REGION_LABELS[@]}" \
    --plotTitle "H3K27me3 signal around TSS (Top${topn} markers per cell line)"
) &
wait

# =======================
# Profiles only (one PDF each)
# =======================
plotProfile \
  --matrixFile "$MAT_DIR/TSS_AC.mat.gz" \
  --outFileName "$FIG_DIR/Top${topn}_CL_TSS_AC_profile.pdf" \
  --numPlotsPerRow 4 \
  --perGroup \
  --legendLocation upper-right \
  --yAxisLabel "mean RPGC" \
  --yMin 0 --yMax "${YMAX_TSS_PROFILE}" \
  --samplesLabel "${SAMPLES_LABEL[@]}" \
  --regionsLabel "${REGION_LABELS[@]}" \
  --colors "${COLORS_PROFILE[@]}" \
  --plotTitle "H3K27ac: aggregate signal (Top${topn} markers per cell line)"

plotProfile \
  --matrixFile "$MAT_DIR/TSS_ME3.mat.gz" \
  --outFileName "$FIG_DIR/Top${topn}_CL_TSS_ME3_profile.pdf" \
  --numPlotsPerRow 4 \
  --perGroup \
  --legendLocation upper-right \
  --yAxisLabel "mean RPGC" \
  --yMin 0 --yMax "${YMAX_TSS_PROFILE}" \
  --samplesLabel "${SAMPLES_LABEL[@]}" \
  --regionsLabel "${REGION_LABELS[@]}" \
  --colors "${COLORS_PROFILE[@]}" \
  --plotTitle "H3K27me3: aggregate signal (Top${topn} markers per cell line)"
