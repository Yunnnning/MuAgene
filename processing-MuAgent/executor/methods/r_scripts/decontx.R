#!/usr/bin/env Rscript
# DecontX wrapper invoked by executor/methods/ambient.py.
#
# Usage:
#   Rscript decontx.R <input_dir> <output_dir> <max_contamination>
#
# Inputs (under <input_dir>):
#   counts.mtx        — cells x genes sparse matrix (Matrix Market). The R side
#                        transposes this internally to genes x cells.
#   barcodes.tsv      — cell barcodes (rownames of counts.mtx)
#   features.tsv      — gene names (colnames of counts.mtx)
#
# Outputs (under <output_dir>):
#   decontaminated.mtx       — cells x genes corrected counts (rounded integers)
#   contamination.tsv        — per-cell contamination fraction (no header)
#   summary.json             — diagnostic summary
suppressPackageStartupMessages({
    library(Matrix)
    library(celda)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript decontx.R <input_dir> <output_dir> [max_contamination]")
}
in_dir <- args[1]
out_dir <- args[2]
max_contam <- if (length(args) >= 3) as.numeric(args[3]) else 1.0

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readMM(file.path(in_dir, "counts.mtx"))   # cells x genes
barcodes <- scan(file.path(in_dir, "barcodes.tsv"), what = character(),
                                 quiet = TRUE, sep = "\n")
features <- scan(file.path(in_dir, "features.tsv"), what = character(),
                                 quiet = TRUE, sep = "\n")

# DecontX expects genes x cells.
counts_t <- t(counts)
rownames(counts_t) <- features
colnames(counts_t) <- barcodes

set.seed(0)
res <- celda::decontX(counts_t, verbose = FALSE)

decon_t <- res$decontXcounts          # genes x cells (sparse, possibly fractional)
contam <- as.numeric(res$contamination)
names(contam) <- barcodes

# Optional cap: refuse to scrub more than `max_contam` from any single cell.
if (!is.null(max_contam) && is.finite(max_contam) && max_contam < 1.0) {
    capped <- pmin(contam, max_contam)
    if (any(capped < contam)) {
        # Re-derive decontaminated counts using the capped fraction:
        #   decon_i = round(orig_i * (1 - capped_i))
        # Apply per cell.
        ratio <- (1 - capped) / pmax(1 - contam, 1e-9)
        ratio[!is.finite(ratio)] <- 1.0
        decon_t <- decon_t %*% Diagonal(x = ratio)
    }
    contam <- capped
}

# Cells x genes again, rounded to non-negative integer counts.
decon <- t(decon_t)
decon@x <- pmax(round(decon@x), 0)
decon <- drop0(decon)

writeMM(decon, file.path(out_dir, "decontaminated.mtx"))
write.table(contam, file.path(out_dir, "contamination.tsv"),
                        quote = FALSE, row.names = FALSE, col.names = FALSE)

summary_obj <- list(
    method = "DecontX",
    n_cells = ncol(counts_t),
    n_genes = nrow(counts_t),
    contamination = list(
        mean = mean(contam),
        median = median(contam),
        min = min(contam),
        max = max(contam),
        q10 = unname(quantile(contam, 0.10)),
        q90 = unname(quantile(contam, 0.90))
    ),
    n_high_contamination = sum(contam > 0.20),
    max_contam_cap = max_contam,
    total_counts_pre = sum(counts_t),
    total_counts_post = sum(decon_t)
)
write_json(summary_obj, file.path(out_dir, "summary.json"),
                     auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("DecontX done: %d cells, median contamination=%.3f\n",
                        length(contam), median(contam)))
