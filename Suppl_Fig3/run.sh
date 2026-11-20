#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export DATA_ROOT="${DATA_ROOT:-$SCRIPT_DIR/inputs}"
export WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
PYTHON="${PYTHON:-python}"
BAM_DIR="${BAM_DIR:-$DATA_ROOT}"
RNA_DIR="${RNA_DIR:-$DATA_ROOT}"

usage() { echo "Usage: bash run.sh [--check-inputs]"; }
if (( $# > 1 )); then usage >&2; exit 2; fi
case "${1:-}" in
  ""|--check-inputs) ;;
  -h|--help) usage; exit 0 ;;
  *) usage >&2; exit 2 ;;
esac

missing=()
need_file() {
  if [[ ! -f "$1" || ! -r "$1" ]]; then missing+=("$1"); fi
}
for cell in GM12878 JJN3 Karpas422 WSU A20; do
  for mark in H3K27ac H3K27me3; do
    need_file "$BAM_DIR/Sc_CMD_${cell}_${mark}_NoDup.bam"
  done
  prefix="$RNA_DIR/CellLineMix_${cell}_RNA_"
  if [[ -f "${prefix}matrix.mtx.gz" && -f "${prefix}features.tsv.gz" && -f "${prefix}barcodes.tsv.gz" ]]; then
    for file in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz; do need_file "$prefix$file"; done
    continue
  fi
  matrix_dir="$RNA_DIR/${cell}.Solo.outGeneFull/filtered"
  if [[ -f "$matrix_dir/genes.tsv" ]]; then
    for file in matrix.mtx genes.tsv barcodes.tsv; do need_file "$matrix_dir/$file"; done
  else
    for file in matrix.mtx features.tsv barcodes.tsv; do
      if [[ -f "$matrix_dir/$file.gz" || ! -f "$matrix_dir/$file" ]]; then
        need_file "$matrix_dir/$file.gz"
      else
        need_file "$matrix_dir/$file"
      fi
    done
  fi
done
need_file "${CELL_CYCLE_GENES:-$DATA_ROOT/regev_lab_cell_cycle_genes.txt}"
need_file "${DNA_BLACKLIST:-$DATA_ROOT/hg38-blacklist.bed.gz}"
need_file "${BLACKLIST_BED_HS:-$DATA_ROOT/hg38-blacklist.bed}"
need_file "${BLACKLIST_BED_MM:-$DATA_ROOT/mm10_grcm39LiftOver-blacklist.v2.bed}"
need_file "${EPD_BED:-$DATA_ROOT/human_epdnew_mjdVA.bed}"
if (( ${#missing[@]} )); then
  printf 'Missing or unreadable inputs:\n' >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi
echo "Found all 30 required input files."
if [[ "${1:-}" == --check-inputs ]]; then exit 0; fi

"$PYTHON" "$SCRIPT_DIR/SupplFig3_0_ProcessingCellLinesMix.py"
"$PYTHON" "$SCRIPT_DIR/SupplFig3_1_HeatmapMarkerGenes.py" --top-n 50 --min-logfc 1
bash "$SCRIPT_DIR/SupplFig3_2_SubsetBAMs.sh"
bash "$SCRIPT_DIR/SupplFig3_3_CallPeaks.sh"
bash "$SCRIPT_DIR/SupplFig3_4_DeepToolsLaunch.sh"
