#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/results}"
DATA_ROOT="${DATA_ROOT:-$SCRIPT_DIR/inputs}"
BAM_DIR="${BAM_DIR:-$DATA_ROOT}"
THREADS="${THREADS:-64}"
SAMPLES=(GM12878 Karpas422 JJN3 WSU A20)
MARKS=(H3K27ac H3K27me3)

[[ "$BAM_DIR" = /* ]] || BAM_DIR="$PWD/$BAM_DIR"
cd -- "$WORK_DIR"
for sample in "${SAMPLES[@]}"; do
  rg_file="${sample}_CountFilteredShared_cell_names.lst"
  if [[ ! -f "$rg_file" ]]; then
    echo "Missing read-group list: $rg_file" >&2
    exit 1
  fi

  for mark in "${MARKS[@]}"; do
    bam="$BAM_DIR/Sc_CMD_${sample}_${mark}_NoDup.bam"
    out="${sample}_${mark}_FiltShared.bam"
    if [[ ! -f "$bam" ]]; then
      echo "Skipping missing BAM: $bam" >&2
      continue
    fi
    echo "Subsetting $bam to $out using $rg_file"
    samtools view --threads "$THREADS" --with-header --read-group-file "$rg_file" --output "$out" "$bam"
    samtools index --threads "$THREADS" "$out"
  done
done
