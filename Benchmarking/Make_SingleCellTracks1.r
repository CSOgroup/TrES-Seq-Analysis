#!/usr/bin/env Rscript

# Run with Rscript Fig1_MakeSingleCellTracks.r from any directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPT_DIR <- dirname(normalizePath(script_file, mustWork = TRUE))

# Editable paths. Use coordinate-sorted locus BAMs with both histone mates retained.
# The BAMs should contain the cells eligible for plotting.
DATA_ROOT <- file.path(SCRIPT_DIR, "inputs")
OUTPUT_ROOT <- file.path(SCRIPT_DIR, "single_cell_stacks")
INPUTS <- c(
  TrES_H3K27ac = file.path(DATA_ROOT, "TrES_H3K27ac.loci.bam"),
  TrES_H3K27me3 = file.path(DATA_ROOT, "TrES_H3K27me3.loci.bam"),
  TrES_RNA = file.path(DATA_ROOT, "TrES_RNA.cohort.loci.bam")
)
SOURCES <- c(TrES_H3K27ac = "", TrES_H3K27me3 = "", TrES_RNA = "")  # URLs or accessions
n_cells_target <- 200L

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "--smoke")) {
  stop("Usage: Rscript Fig1_MakeSingleCellTracks.r [--smoke]", call. = FALSE)
}
missing <- INPUTS[!file.exists(INPUTS)]
if (length(missing)) {
  stop("Missing BAM inputs:\n", paste(missing, collapse = "\n"), call. = FALSE)
}
if (identical(args, "--smoke")) {
  OUTPUT_ROOT <- tempfile("fig1_stacks_", tmpdir = dirname(tempdir()))
}

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(Rsamtools)
  library(ggplot2)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Locus coordinates are hg38, 1-based inclusive.
REGIONS <- list(
  MYC   = c("chr8", 127723000, 127753000),
  BCL6  = c("chr3", 187713000, 187755000),
  IRF4  = c("chr6",      336448,      436683),
  PRDM1 = c("chr6", 106062900, 106133200)
)
if (identical(args, "--smoke")) REGIONS <- REGIONS["MYC"]

sample_names <- c("TrES H3K27ac", "TrES H3K27me3", "TrES RNA")
sample_colors <- c(
  "TrES H3K27ac" = "#3D66AC",
  "TrES H3K27me3" = "#AF3134",
  "TrES RNA" = "#449952"
)

read_tag <- function(alignments, tag) {
  ids <- mcols(alignments)[[tag]]
  if (is.null(ids)) return(rep(NA_character_, length(alignments)))
  ids[!nzchar(ids)] <- NA_character_
  ids
}

# Prefer CB on either mate. RG is a fallback for BAMs that use it for cell IDs.
cell_ids <- function(reads, mates = NULL) {
  cells <- read_tag(reads, "CB")
  groups <- read_tag(reads, "RG")
  if (!is.null(mates)) {
    mate_cells <- read_tag(mates, "CB")
    mate_groups <- read_tag(mates, "RG")
    missing <- is.na(cells)
    cells[missing] <- mate_cells[missing]
    if (any(cells != mate_cells, na.rm = TRUE) ||
        any(is.na(cells) & groups != mate_groups, na.rm = TRUE)) {
      stop("Histone mates have different cell IDs", call. = FALSE)
    }
    missing <- is.na(groups)
    groups[missing] <- mate_groups[missing]
  }
  missing <- is.na(cells)
  cells[missing] <- groups[missing]
  cells
}

# Read each compact BAM once, keeping both mates before filtering by locus.
read_segments <- function(bam_file, is_histone) {
  message("  BAM: ", bam_file)
  bf <- BamFile(bam_file, index = character())
  param <- ScanBamParam(tag = c("CB", "RG"), flag = scanBamFlag(
    isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE
  ))

  if (is_histone) {
    galp <- readGAlignmentPairs(bf, param = param)
    if (!length(galp)) {
      stop("No complete histone read pairs in ", bam_file, call. = FALSE)
    }
    cell <- cell_ids(GenomicAlignments::first(galp), GenomicAlignments::last(galp))
    frags <- granges(galp)
    chrom <- as.character(seqnames(frags))
    xs <- start(frags)
    xe <- end(frags)
  } else {
    gal <- readGAlignments(bf, param = param)
    exon_list <- cigarRangesAlongReferenceSpace(
      cigar(gal),
      pos = start(gal),
      ops = c("M", "=", "X"),
      drop.empty.ranges = TRUE
    )
    exon_ir <- unlist(exon_list)
    chrom <- rep(as.character(seqnames(gal)), lengths(exon_list))
    xs <- start(exon_ir)
    xe <- end(exon_ir)
    cell <- rep(cell_ids(gal), lengths(exon_list))
  }
  keep <- !is.na(cell)
  if (!any(keep)) {
    stop("No segments with CB or RG cell IDs in ", bam_file, call. = FALSE)
  }
  data.frame(chrom = chrom[keep], start = xs[keep], end = xe[keep], cell = cell[keep])
}

get_read_df <- function(segments, sample_name, chrom, start, end) {
  overlaps <- segments$chrom == chrom & segments$start <= end & segments$end >= start
  segments <- segments[overlaps, , drop = FALSE]
  if (!nrow(segments)) return(NULL)

  counts <- table(segments$cell)
  ranked <- names(counts)[order(-as.integer(counts), names(counts), method = "radix")]
  kept <- head(ranked, n_cells_target)
  rank <- match(segments$cell, kept)
  segments <- segments[!is.na(rank), , drop = FALSE]
  message("  ", sample_name, ": ", length(kept), " of ", length(counts), " cells")
  # Count whole overlapping segments, then clip them for display.
  data.frame(
    sample = sample_name,
    xstart = pmax(segments$start, start),
    xend   = pmin(segments$end, end),
    y      = rank[!is.na(rank)]
  )
}

segments <- Map(read_segments, INPUTS, c(TRUE, TRUE, FALSE))
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUTPUT_ROOT)) stop("Cannot create output directory: ", OUTPUT_ROOT, call. = FALSE)

for (locus in names(REGIONS)) {
  message("Locus: ", locus)

  chrom <- REGIONS[[locus]][1]
  start <- as.numeric(REGIONS[[locus]][2])
  end   <- as.numeric(REGIONS[[locus]][3])

  df_all <- do.call(rbind, Map(function(rows, name) {
    get_read_df(rows, name, chrom, start, end)
  }, segments, sample_names))
  if (is.null(df_all)) {
    warning("No segments overlap ", locus, call. = FALSE)
    next
  }

  df_all$sample <- factor(df_all$sample, levels = sample_names)

  p <- ggplot(df_all) +
    geom_segment(
      aes(x = xstart, xend = xend, y = y, yend = y, color = sample),
      linewidth = 0.4,
      alpha     = 0.85,
      lineend   = "round"
    ) +
    facet_grid(
      rows   = vars(sample),
      scales = "free_y",
      switch = "y"
    ) +
    scale_color_manual(values = sample_colors) +
    scale_y_reverse() +
    scale_x_continuous(
      limits = c(start, end),
      breaks = c(start, end),
      labels = label_comma(accuracy = 1)
    ) +
    labs(
      x = paste0("Genomic position (", chrom, ")"),
      y = NULL,
      title = locus
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title        = element_text(hjust = 0, face = "bold", size = 12),
      panel.grid        = element_blank(),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.title.y      = element_blank(),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
      strip.background  = element_blank(),
      legend.position   = "none",
      panel.spacing.y   = unit(0.25, "lines"),
      plot.margin       = margin(5.5, 5.5, 5.5, 5.5, "pt")
    )

  out_file <- file.path(OUTPUT_ROOT, paste0("reads_fragments_", locus, "_stacked_top.pdf"))
  ggsave(out_file, p, width = 8, height = 6, useDingbats = FALSE)
  cat("  → Saved:", out_file, "\n")
}

cat("\nAll loci processed.\n")
