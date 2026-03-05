#!/bin/bash
set -euo pipefail

THREADS=64
BAM_DIR="../processed"

# Map "short" sample names (used in .lst + output) to BAM prefixes
declare -A PREFIX
PREFIX[A20]="Sc_CMD_A20"
PREFIX[GM12878]="Sc_CMD_GM12878"
PREFIX[JJN2]="Sc_CMD_JJN2"
PREFIX[Karpas422]="Sc_CMD_Karpas422"
PREFIX[WSU]="Sc_CMD_WSU"

# Samples & marks we want
SAMPLES=(GM12878 Karpas422 JJN2 WSU A20)
MARKS=(H3K27ac H3K27me3)

for sample in "${SAMPLES[@]}"; do
    rg_file="${sample}_CountFilteredShared_cell_names.lst"
    if [[ ! -f "${rg_file}" ]]; then
        echo "ERROR: missing read-group file ${rg_file}" >&2
        exit 1
    fi

    for mark in "${MARKS[@]}"; do
        in_bam="${BAM_DIR}/${PREFIX[$sample]}_${mark}_NoDup.bam"
        out_bam="${sample}_${mark}_FiltShared.bam"

        # Only run if input exists and output not already there
        if [[ ! -f "${in_bam}" ]]; then
            echo "WARNING: input BAM not found, skipping: ${in_bam}" >&2
            continue
        fi

        echo ">>> Subsetting ${in_bam}"
        echo "    using RG list: ${rg_file}"
        echo "    -> ${out_bam}"

        samtools view \
            --threads "${THREADS}" \
            --with-header \
            --read-group-file "${rg_file}" \
            --output "${out_bam}" \
            "${in_bam}"

        samtools index --threads "${THREADS}" "${out_bam}"
    done
done
