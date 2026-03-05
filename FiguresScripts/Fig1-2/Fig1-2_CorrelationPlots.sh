#DNA
#!/bin/bash
set -euo pipefail

# ---------- Params you can tweak ----------
BWS="*bw"   # bigWig glob
BLACKLIST="/mnt/dataFast/ahrmad/hg38-blacklist.v2.bed"
CORES_TOTAL=90
CORES_PER_JOB=16              # 5 jobs × 16 = 80 cores
BINSIZES=(10000)
# ------------------------------------------

NJOBS=${#BINSIZES[@]}
if (( CORES_PER_JOB * NJOBS > CORES_TOTAL )); then
  echo "Warning: requested ${CORES_PER_JOB}×${NJOBS}=${CORES_PER_JOB*NJOBS} cores > ${CORES_TOTAL} total."
  echo "Proceeding anyway."
fi

run_one() {
  local bin="$1"
  if [[ -z "${bin}" ]]; then
    echo "ERROR: empty bin size" >&2
    return 2
  fi

  echo ">>> Running bin size: ${bin}"
  local OUTPREFIX="bin${bin}_H2RvScKDMA_ptm"

  # multiBigwigSummary
  multiBigwigSummary bins \
    -b $BWS \
    -o "${OUTPREFIX}.npz" \
    --smartLabels \
    --binSize "${bin}" \
    --blackListFileName "${BLACKLIST}" \
    --numberOfProcessors "${CORES_PER_JOB}"

  # scatterplot
  plotCorrelation \
    -in "${OUTPREFIX}.npz" \
    -c spearman \
    -p scatterplot \
    -o "${OUTPREFIX}_scatter.pdf" \
    --log1p

  # heatmap
  plotCorrelation \
    -in "${OUTPREFIX}.npz" \
    -c spearman \
    -p heatmap \
    -o "${OUTPREFIX}_heatmap.pdf" \
    --plotNumbers

  echo ">>> Finished bin size: ${bin}"
}

export -f run_one
export BWS BLACKLIST CORES_PER_JOB

# One job per bin size
parallel -j "${NJOBS}" run_one ::: "${BINSIZES[@]}"




#RNA

#!/bin/bash
set -euo pipefail

##############################
# params
##############################

BWS="*.bw"              
MAX_CORES=80            
BINSIZES=(10000)

# Label for outputs 
#OUTLABEL="TrESCorrelation_RNA"

##############################
# Derived settings
##############################

NJOBS=${#BINSIZES[@]}
if (( NJOBS == 0 )); then
  echo "ERROR: no bin sizes defined" >&2
  exit 1
fi

# Choose cores per job so that all jobs can run at once within MAX_CORES
CORES_PER_JOB=$(( MAX_CORES / NJOBS ))
if (( CORES_PER_JOB < 1 )); then
  CORES_PER_JOB=1
fi

TOTAL_USED=$(( CORES_PER_JOB * NJOBS ))

echo "=== deepTools RNA correlation ==="
echo "Found ${NJOBS} bin sizes: ${BINSIZES[*]}"
echo "Using ${CORES_PER_JOB} threads per job"
echo "Total threads used: ${TOTAL_USED} (limit ${MAX_CORES})"
echo "BigWigs: $BWS"
echo "Out label: ${OUTLABEL}"
echo "================================="

run_one() {
  local bin="$1"
  if [[ -z "${bin}" ]]; then
    echo "ERROR: empty bin size" >&2
    return 2
  fi

  local OUTPREFIX="bin${bin}_${OUTLABEL}"

  echo ">>> [bin=${bin}] Output prefix: ${OUTPREFIX}"

  # multiBigwigSummary (skipped if npz already exists)
  if [[ -f "${OUTPREFIX}.npz" ]]; then
    echo ">>> [bin=${bin}] ${OUTPREFIX}.npz exists, skipping multiBigwigSummary"
  else
    echo ">>> [bin=${bin}] Running multiBigwigSummary"
    multiBigwigSummary bins \
      -b ${BWS} \
      -o "${OUTPREFIX}.npz" \
      --smartLabels \
      --binSize "${bin}" \
      --numberOfProcessors "${CORES_PER_JOB}"
  fi

  if [[ ! -f "${OUTPREFIX}.npz" ]]; then
    echo ">>> [bin=${bin}] WARNING: ${OUTPREFIX}.npz missing, skipping plots"
    return 0
  fi

  # scatterplot
  echo ">>> [bin=${bin}] Plotting scatter"
  plotCorrelation \
    -in "${OUTPREFIX}.npz" \
    -c spearman \
    -p scatterplot \
    -o "${OUTPREFIX}_scatter.pdf" \
    --log1p

  # heatmap
  echo ">>> [bin=${bin}] Plotting heatmap"
  plotCorrelation \
    -in "${OUTPREFIX}.npz" \
    -c spearman \
    -p heatmap \
    -o "${OUTPREFIX}_heatmap.pdf" \
    --plotNumbers

  echo ">>> [bin=${bin}] Done"
}

export -f run_one
export BWS OUTLABEL CORES_PER_JOB

parallel -j "${NJOBS}" run_one ::: "${BINSIZES[@]}"




















































#ac
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_ac.npz --smartLabels --binSize 10000 --blackListFileName /mnt/dataFast/ahrmad/hg38-blacklist.v2.bed --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_ac.npz --corMethod spearman --whatToPlot scatterplot --plotFile 10xTrES_ac_scatter.pdf --log1p
plotCorrelation --corData H2RvScKDMA_ac.npz --corMethod spearman --whatToPlot heatmap --plotFile 10xTrES_ac_heatmap.pdf --plotNumbers

#me3
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_me3.npz --smartLabels --binSize 10000 --blackListFileName /mnt/dataFast/ahrmad/hg38-blacklist.v2.bed --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_me3.npz --corMethod spearman --whatToPlot scatterplot --plotFile 10xTrES_me3_scatter.pdf --log1p
plotCorrelation --corData H2RvScKDMA_me3.npz --corMethod spearman --whatToPlot heatmap --plotFile 10xTrES_me3_heatmap.pdf --plotNumbers


#me3 and ac
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_ptm.npz --smartLabels --binSize 10000 --blackListFileName /mnt/dataFast/ahrmad/hg38-blacklist.v2.bed --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_ptm.npz --corMethod spearman --whatToPlot scatterplot --plotFile 10xTrES_h3k27_scatter.pdf --log1p
plotCorrelation --corData H2RvScKDMA_ptm.npz --corMethod spearman --whatToPlot heatmap --plotFile 10xTrES_h3k27_heatmap.pdf --plotNumbers

#h3k9me3
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_me3.npz --smartLabels --binSize 20000 --blackListFileName /mnt/dataFast/ahrmad/hg38-blacklist.v2.bed --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_me3.npz --corMethod spearman --whatToPlot scatterplot --plotFile CHiPTrES_h3k9me3_scatter.pdf --log1p
plotCorrelation --corData H2RvScKDMA_me3.npz --corMethod spearman --whatToPlot heatmap --plotFile CHiPTrES_h3k9me3_heatmap.pdf --plotNumbers


#me3 and ac and h3k9me3
multiBigwigSummary bins --bwfiles *bw --outFileName H2R.npz --smartLabels --binSize 20000 --blackListFileName /mnt/dataFast/ahrmad/hg38-blacklist.v2.bed --numberOfProcessors 64 
plotCorrelation --corData H2R.npz --corMethod spearman --whatToPlot scatterplot --plotFile TrES_DNA_scatter_20kb.pdf --log1p
plotCorrelation --corData H2R.npz --corMethod spearman --whatToPlot heatmap --plotFile TrES_DNA_heatmap_20kb.pdf --plotNumbers








bigwigCompare --bigwig1 H2RRNA_Human.stranded_Signal.Unique.str1.out.bw --bigwig2 H2RRNA_Human.stranded_Signal.Unique.str2.out.bw -o H2RRNA_Human.unstrandedCALC.bw --operation add --skipNAs --binSize 20 -p 64 -v 
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_RNA.npz --smartLabels --binSize 1000 --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_RNA.npz --corMethod spearman --whatToPlot heatmap --plotFile 10xTrES_RNA_heatmap.pdf --plotNumbers




#rna
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_RNA.npz --smartLabels --binSize 1000 --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_RNA.npz --corMethod spearman --whatToPlot scatterplot --plotFile 10xTrES_RNA_scatter.pdf --log1p
plotCorrelation --corData H2RvScKDMA_RNA.npz --corMethod spearman --whatToPlot heatmap --plotFile 10xTrES_RNA_heatmap.pdf --plotNumbers


#rna with DT
multiBigwigSummary bins --bwfiles *bw --outFileName H2RvScKDMA_RNA_WithDT.npz --smartLabels --binSize 1000 --numberOfProcessors 64 
plotCorrelation --corData H2RvScKDMA_RNA_WithDT.npz --corMethod spearman --whatToPlot scatterplot --plotFile 10xTrES_RNA_scatter_WithDT.pdf --log1p
plotCorrelation --corData H2RvScKDMA_RNA_WithDT.npz --corMethod spearman --whatToPlot heatmap --plotFile 10xTrES_RNA_heatmap_WithDT.pdf --plotNumbers






#CELL LINES

#!/usr/bin/env bash
set -euo pipefail

export THREADS_PER_TASK=12

do_one() {
  e="$1"
  samtools view -@ "$THREADS_PER_TASK" \
    --with-header \
    --read-group-file "Sc_CMD_${e}_rna_CountFiltered_cell_names.lst" \
    -o "rna_${e}_Filt.bam" \
    "Sc_CMr_${e}_UMI.bam"

  samtools index -@ "$THREADS_PER_TASK" "rna_${e}_Filt.bam"
}
export -f do_one

# Adjust -j to how many you want to run concurrently
parallel -j 64 do_one ::: GM12878 JJN2 Karpas422 WSU



#!/usr/bin/env bash
set -euo pipefail

export THREADS_PER_TASK=12

do_one() {
  e="$1"
  hg38_eff_size=2913022398

    bamCoverage -p 16 -bs 10 -b rna_${e}_Filt.bam -o rna_${e}_Filt_CPM_BothStrands.bw -of bigwig --normalizeUsing CPM --effectiveGenomeSize ${hg38_eff_size}

    STAR --runMode inputAlignmentsFromBAM \
    --inputBAMfile rna_${e}_Filt.bam \
    --outWigType bedGraph \
    --outWigNorm RPM \
    --outWigStrand Stranded \
    --outWigReferencesPrefix chr \
    --outFileNamePrefix rna_${e}_Filt

    rm rna_${e}_FiltSignal.Unique.str*.out.bg
    sortBed -i rna_${e}_FiltSignal.UniqueMultiple.str1.out.bg > rna_${e}_FiltSignal.Str1Sorted.bg
    sortBed -i rna_${e}_FiltSignal.UniqueMultiple.str2.out.bg > rna_${e}_FiltSignal.Str2Sorted.bg
    bedGraphToBigWig rna_${e}_FiltSignal.Str1Sorted.bg /mnt/dataFast/ahrmad/hg38.chrom.sizes rna_${e}_Filt_CPM_minus.bw
    bedGraphToBigWig rna_${e}_FiltSignal.Str2Sorted.bg /mnt/dataFast/ahrmad/hg38.chrom.sizes rna_${e}_Filt_CPM_plus.bw
}
export -f do_one

# Adjust -j to how many you want to run concurrently
parallel -j 70 do_one ::: GM12878 JJN2 Karpas422 WSU
