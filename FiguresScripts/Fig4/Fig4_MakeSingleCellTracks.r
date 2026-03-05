#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(Rsamtools)
  library(rtracklayer)
  library(ggplot2)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE, scipen = 999)

###############################################################
## USER TOGGLES
###############################################################

## 1) How to pick cells?
##    "top"    = top N cells by segment/fragment count in region
##    "random" = randomly sample N cells (uniform over cells)
cell_selection_mode <- "top"   # "top" or "random"
n_cells_target <- 200L

## If using random cell selection, optional seed for reproducibility
random_seed <- 1L
if (cell_selection_mode == "random" && !is.null(random_seed)) {
  set.seed(random_seed)
}

###############################################################
## REGION & SAMPLE DEFINITIONS
###############################################################

REGIONS <- list(
  Region1  = c("chr19", 36000000, 39000000),
  Region2  = c("chr17", 49400000, 51540000)
)

## Base sample definitions matching your actual filenames
sample_names <- c(
  "TrES H3K9me3",
  "TrES H3K27me3",
  "TrES RNA"
)

bam_files_fixed <- c(
  "Sc_K9D_H3K9me3_Filt.bam",
  "Sc_K9D_H3K27me3_Filt.bam",
  "Sc_K9D_RNA_Filt.bam"
)

sample_colors <- c(
  "TrES H3K9me3"     = "#E68619",  # Your new orange
  "TrES H3K27me3"    = "#AF3134",  # Brick Red
  "TrES RNA"         = "#449952"   # Forest Green
)

###############################################################
## CORE FUNCTION: per-cell segments for one BAM & locus
###############################################################

get_read_df <- function(bam_file,
                        sample_name,
                        chrom,
                        start,
                        end,
                        n_cells = 200L,
                        cell_mode = c("top", "random")) {

  cell_mode <- match.arg(cell_mode)

  message("  BAM: ", bam_file, " (", sample_name, ")")
  if (!file.exists(bam_file)) {
    message("    [WARN] File does not exist, skipping.")
    return(NULL)
  }

  region <- GRanges(seqnames = chrom,
                    ranges   = IRanges(start = start, end = end))

  bf <- BamFile(bam_file)
  # CHANGED: tag is now CB
  param <- ScanBamParam(which = region, tag = "CB")

  # CHANGED: Broadened to include H3K9 and H3K27
  is_histone <- grepl("H3K", sample_name, fixed = TRUE)

  if (is_histone) {
    ############################################################
    ## Histone: paired-end FRAGMENTS
    ############################################################
    message("    Using fragments (read pairs).")
    galp <- readGAlignmentPairs(bf, param = param)

    if (length(galp) == 0L) {
      message("    [INFO] No read pairs in region.")
      return(NULL)
    }

    first_gal <- GenomicAlignments::first(galp)
    cb <- mcols(first_gal)[["CB"]] # CHANGED to CB

    keep <- !is.na(cb)
    galp <- galp[keep]
    cb   <- cb[keep]

    if (length(galp) == 0L) {
      message("    [INFO] No read pairs with CB tag.")
      return(NULL)
    }

    frags <- granges(galp)
    xs <- start(frags)
    xe <- end(frags)

  } else {
    ############################################################
    ## RNA: plot each exon chunk separately
    ############################################################
    message("    Using per-read exon segments (no intron spans).")

    gal <- import(bf, param = param)

    if (length(gal) == 0L) {
      message("    [INFO] No reads in region.")
      return(NULL)
    }

    cb <- mcols(gal)$CB # CHANGED to CB
    keep <- !is.na(cb)
    gal  <- gal[keep]
    cb   <- cb[keep]

    if (length(gal) == 0L) {
      message("    [INFO] No reads with CB tag.")
      return(NULL)
    }

    exon_list <- cigarRangesAlongReferenceSpace(
      cigar(gal),
      pos = start(gal),
      ops = c("M", "=", "X"),
      drop.empty.ranges = TRUE
    )

    if (length(exon_list) == 0L) {
      message("    [INFO] No exon blocks found.")
      return(NULL)
    }

    exon_ir <- unlist(exon_list)
    xs <- start(exon_ir)
    xe <- end(exon_ir)

    cb <- rep(cb, lengths(exon_list))
  }

  if (length(xs) == 0L) {
    message("    [INFO] No segments to plot.")
    return(NULL)
  }

  ############################################################
  ## Select cells (top N or random N)
  ############################################################

  cb_counts      <- table(cb)
  all_cells      <- names(cb_counts)
  n_cells_total  <- length(all_cells)
  n_keep         <- min(n_cells, n_cells_total)

  if (cell_mode == "top") {
    cb_order <- names(sort(cb_counts, decreasing = TRUE))
    cb_keep  <- cb_order[seq_len(n_keep)]
    message("    Cells in region: ", n_cells_total,
            " | keeping TOP ", n_keep, " cells by count.")
  } else {
    cb_keep  <- sample(all_cells, n_keep)
    message("    Cells in region: ", n_cells_total,
            " | keeping RANDOM ", n_keep, " cells.")
  }

  keep_cb <- cb %in% cb_keep
  xs <- xs[keep_cb]
  xe <- xe[keep_cb]
  cb <- cb[keep_cb]

  if (length(xs) == 0L) {
    message("    [INFO] No segments left after CB selection.")
    return(NULL)
  }

  ## Assign y-position per cell
  if (cell_mode == "top") {
    ordering <- cb_keep
  } else {
    ordering <- sort(cb_keep)
  }
  cb_to_row <- setNames(seq_along(ordering), ordering)
  yvals <- cb_to_row[as.character(cb)]

  data.frame(
    sample = sample_name,
    xstart = as.numeric(xs),
    xend   = as.numeric(xe),
    y      = as.numeric(yvals)
  )
}

###############################################################
## MAIN LOOP: one plot per locus
###############################################################

for (locus in names(REGIONS)) {
  cat("\n==============================\n")
  cat("Locus:", locus, "\n")
  cat("==============================\n")

  chrom <- REGIONS[[locus]][1]
  start <- as.numeric(REGIONS[[locus]][2])
  end   <- as.numeric(REGIONS[[locus]][3])

  all_df_list <- vector("list", length(bam_files_fixed))

  for (i in seq_along(bam_files_fixed)) {
    all_df_list[[i]] <- get_read_df(
      bam_file    = bam_files_fixed[i],
      sample_name = sample_names[i],
      chrom       = chrom,
      start       = start,
      end         = end,
      n_cells     = n_cells_target,
      cell_mode   = cell_selection_mode
    )
  }

  non_null <- !vapply(all_df_list, is.null, logical(1))
  if (!any(non_null)) {
    message("[WARN] No data to plot for locus ", locus, ". Skipping.")
    next
  }

  df_all <- do.call(rbind, all_df_list[non_null])
  df_all$sample <- factor(df_all$sample, levels = sample_names)

  ############################################################
  ## Plot
  ############################################################

  p <- ggplot(df_all) +
    geom_segment(
      aes(x = xstart, xend = xend, y = y, yend = y, color = sample),
      linewidth = 1,
      alpha     = 0.85,
      lineend   = "round"
    ) +
    facet_grid(
      rows   = vars(sample),
      scales = "free_y",
      switch = "y"
    ) +
    scale_color_manual(values = sample_colors) +
    scale_x_continuous(
      limits = c(start, end),
      expand = c(0, 0),
      labels = label_comma()
    ) +
    labs(
      x = paste0("Genomic position (", chrom, ")"),
      y = NULL,
      title = paste0("Locus: ", locus)
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title        = element_text(hjust = 0, face = "bold"),
      panel.grid        = element_blank(),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      strip.background  = element_blank(),
      legend.position   = "none",
      panel.spacing.y   = unit(0.5, "lines")
    )

  mode_label <- if (cell_selection_mode == "top") "top" else "random"
  out_file <- paste0("reads_fragments_", locus, "_stacked_", mode_label, ".pdf")
  ggsave(out_file, p, width = 10, height = 7)
  cat("  → Saved:", out_file, "\n")
}

cat("\nAll loci processed.\n")