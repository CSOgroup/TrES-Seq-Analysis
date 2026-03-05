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

## 1) Include NanoCNT samples?
include_nanocnt <- FALSE  # set to TRUE to include NanoCNT, FALSE to drop them

## 2) How to pick cells?
##    "top"    = top N cells by segment/fragment count in region
##    "random" = randomly sample N cells (uniform over cells)
cell_selection_mode <- "random"   # "top" or "random"
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
  MYC   = c("chr8", 127723000, 127753000),
  BCL6  = c("chr3", 187713000, 187755000),
  IRF4  = c("chr6",      336448,      436683),
  PRDM1 = c("chr6", 106062900, 106133200)
)

## Base sample definitions (full set)
base_sample_names <- c(
  "TrES H3K27ac",
  "TrES H3K27me3",
  "NanoCNT H3K27ac",
  "NanoCNT H3K27me3",
  "TrES RNA"
)

base_bam_templates <- c(
  "TrES_H3K27ac.subset.FX.%s.bam",
  "TrES_H3K27me3.subset.FX.%s.bam",
  "NanoCNT_H3K27ac.subset.FX.%s.bam",
  "NanoCNT_H3K27me3.subset.FX.%s.bam",
  "TrES_RNA.subset.FX.%s.bam"
)

base_sample_colors <- c(
  "TrES H3K27ac"     = "#3D66AC",  # 61,102,172
  "TrES H3K27me3"    = "#AF3134",  # 175,49,52
  "NanoCNT H3K27ac"  = "#63BDC2",  # 99,189,194
  "NanoCNT H3K27me3" = "#E1706C",  # 225,112,108
  "TrES RNA"         = "#449952"   # 68,153,82
)

## Apply NanoCNT toggle
if (include_nanocnt) {
  sample_names   <- base_sample_names
  bam_templates  <- base_bam_templates
} else {
  keep <- !grepl("NanoCNT", base_sample_names, fixed = TRUE)
  sample_names   <- base_sample_names[keep]
  bam_templates  <- base_bam_templates[keep]
}
sample_colors <- base_sample_colors[sample_names]

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
  param <- ScanBamParam(which = region, tag = "FX")

  is_h3k27 <- grepl("H3K27", sample_name, fixed = TRUE)

  if (is_h3k27) {
    ############################################################
    ## H3K27ac / H3K27me3: paired-end FRAGMENTS
    ############################################################
    message("    Using fragments (read pairs).")
    galp <- readGAlignmentPairs(bf, param = param)

    if (length(galp) == 0L) {
      message("    [INFO] No read pairs in region.")
      return(NULL)
    }

    first_gal <- GenomicAlignments::first(galp)
    fx <- mcols(first_gal)[["FX"]]

    keep <- !is.na(fx)
    galp <- galp[keep]
    fx   <- fx[keep]

    if (length(galp) == 0L) {
      message("    [INFO] No read pairs with FX.")
      return(NULL)
    }

    frags <- granges(galp)
    xs <- start(frags)
    xe <- end(frags)

  } else {
    ############################################################
    ## RNA: plot each exon chunk separately (no intron-spanning)
    ############################################################
    message("    Using per-read exon segments (no intron spans).")

    gal <- import(bf, param = param)

    if (length(gal) == 0L) {
      message("    [INFO] No reads in region.")
      return(NULL)
    }

    fx <- mcols(gal)$FX
    keep <- !is.na(fx)
    gal  <- gal[keep]
    fx   <- fx[keep]

    if (length(gal) == 0L) {
      message("    [INFO] No reads with FX.")
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

    fx <- rep(fx, lengths(exon_list))
  }

  if (length(xs) == 0L) {
    message("    [INFO] No segments to plot.")
    return(NULL)
  }

  ############################################################
  ## Select cells (top N or random N)
  ############################################################

  fx_counts      <- table(fx)
  all_cells      <- names(fx_counts)
  n_cells_total  <- length(all_cells)
  n_keep         <- min(n_cells, n_cells_total)

  if (cell_mode == "top") {
    fx_order <- names(sort(fx_counts, decreasing = TRUE))
    fx_keep  <- fx_order[seq_len(n_keep)]
    message("    Cells in region: ", n_cells_total,
            " | keeping TOP ", n_keep, " cells by count.")
  } else {
    fx_keep  <- sample(all_cells, n_keep)
    message("    Cells in region: ", n_cells_total,
            " | keeping RANDOM ", n_keep, " cells.")
  }

  keep_fx <- fx %in% fx_keep
  xs <- xs[keep_fx]
  xe <- xe[keep_fx]
  fx <- fx[keep_fx]

  if (length(xs) == 0L) {
    message("    [INFO] No segments left after FX selection.")
    return(NULL)
  }

  ## Assign y-position per cell (compressed rows, in a stable order)
  ## For "top": preserve abundance order; for "random": sort by FX ID.
  if (cell_mode == "top") {
    ordering <- fx_keep
  } else {
    ordering <- sort(fx_keep)
  }
  fx_to_row <- setNames(seq_along(ordering), ordering)
  yvals <- fx_to_row[as.character(fx)]

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

  ## BAMs for this locus, in the current sample order
  bam_files <- sprintf(bam_templates, locus)

  all_df_list <- vector("list", length(bam_files))

  for (i in seq_along(bam_files)) {
    all_df_list[[i]] <- get_read_df(
      bam_file    = bam_files[i],
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

  if (nrow(df_all) == 0L) {
    message("[WARN] Empty data.frame for locus ", locus, ". Skipping.")
    next
  }

  df_all$sample <- factor(df_all$sample, levels = sample_names)

  ############################################################
  ## Plot
  ############################################################

  p <- ggplot(df_all) +
    geom_segment(
      aes(x = xstart, xend = xend, y = y, yend = y, color = sample),
      linewidth = 0.4,
      alpha     = 0.85,
      lineend   = "round"
    ) +
    facet_grid(
      rows   = vars(sample),
      cols   = vars(),       # single column
      scales = "free_y",
      switch = "y"
    ) +
    scale_color_manual(values = sample_colors) +
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

  mode_label <- if (cell_selection_mode == "top") "top" else "random"
  out_file <- paste0("reads_fragments_", locus, "_stacked_", mode_label, ".pdf")
  ggsave(out_file, p, width = 8, height = 6, useDingbats = FALSE)
  cat("  → Saved:", out_file, "\n")
}

cat("\nAll loci processed.\n")
