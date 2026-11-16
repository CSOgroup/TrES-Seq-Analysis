#!/usr/bin/env Rscript
# DESeq2 contrasts on donor pseudobulks or disjoint exploratory cell subsets.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(4, 6), args[1] %in% c("pair", "zone"),
          args[3] %in% c("paired", "random"))
pair <- args[1] == "pair"
expected_args <- if (pair) 4L else 6L
stopifnot(length(args) == expected_args)
folder <- args[2]
paired <- args[3] == "paired"
alpha <- as.numeric(args[4])
suppressPackageStartupMessages(library(DESeq2))

counts <- as.matrix(read.delim(gzfile(file.path(folder, "counts.tsv.gz")),
                              row.names = 1, check.names = FALSE, quote = "", comment.char = ""))
meta <- read.delim(file.path(folder, "units.tsv"), row.names = 1,
                   check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
if (!identical(colnames(counts), rownames(meta))) stop("Counts/metadata order mismatch")
if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)) ||
    max(counts) > .Machine$integer.max) stop("Raw integer counts in R integer range required")
storage.mode(counts) <- "integer"

variable <- if (pair) "condition" else "zone"
levels <- if (pair) c("rest", "cluster") else c("LZ", "DZ")
meta[[variable]] <- factor(meta[[variable]], levels = levels)
if (anyNA(meta[[variable]]) || any(table(meta[[variable]]) < 2L)) stop("Too few replicates per group")
if (paired) {
  meta$donor <- factor(meta$donor)
  if (nlevels(meta$donor) < 2L || any(table(meta$donor, meta[[variable]]) != 1L)) {
    stop("One complete pair per donor required")
  }
}
formula <- reformulate(c(if (paired) "donor", variable))
design <- model.matrix(formula, meta)
if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Non-estimable design")

keep <- if (pair) rep(TRUE, nrow(counts)) else rowSums(counts >= as.integer(args[5])) >= as.integer(args[6])
if (!pair && sum(keep) < 2L) stop("Too few genes pass the count filter")
if (any(colSums(counts[keep, , drop = FALSE]) == 0)) stop("Empty pseudobulk after count filtering")
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = meta, design = formula)
dds <- if (pair) DESeq(dds, quiet = TRUE) else DESeq(dds, sfType = "poscounts", quiet = TRUE)
res <- results(dds, contrast = c(variable, levels[2], levels[1]), alpha = alpha)
if (pair) {
  out <- data.frame(gene = rownames(res), logFC = res$log2FoldChange,
                    baseMean = res$baseMean, logCPM = log2(res$baseMean + 0.5),
                    PValue = res$pvalue, FDR = res$padj,
                    analysis_unit = if (paired) "biological donor (paired)" else "random cell subset; exploratory")
} else {
  out <- data.frame(gene = rownames(counts), logFC = NA_real_, baseMean = NA_real_,
                    PValue = NA_real_, FDR = NA_real_, count_filter_pass = keep,
                    lfcSE = NA_real_, stat = NA_real_, dispersion = NA_real_)
  pos <- match(rownames(res), out$gene)
  out$logFC[pos] <- res$log2FoldChange
  out$baseMean[pos] <- res$baseMean
  out$PValue[pos] <- res$pvalue
  out$FDR[pos] <- res$padj
  out$lfcSE[pos] <- res$lfcSE
  out$stat[pos] <- res$stat
  out$dispersion[match(rownames(dds), out$gene)] <- dispersions(dds)
}
write.table(out, gzfile(file.path(folder, "DESeq2_results.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
