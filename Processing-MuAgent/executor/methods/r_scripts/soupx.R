#!/usr/bin/env Rscript
# SoupX wrapper invoked by executor/methods/ambient.py.
#
# Usage:
#   Rscript soupx.R <input_dir> <output_dir> [max_contamination]
#
# Inputs (under <input_dir>):
#   filtered_counts.mtx   — cells x genes sparse matrix
#   filtered_barcodes.tsv
#   raw_counts.mtx        — droplets x genes sparse matrix (raw, includes empties)
#   raw_barcodes.tsv
#   features.tsv          — shared gene order between filtered & raw
#
# Outputs (under <output_dir>):
#   decontaminated.mtx       — cells x genes corrected counts (integers)
#   contamination.tsv        — per-cell contamination fraction (no header)
#   summary.json             — diagnostic summary
#
# SoupX semantics:
#   - Estimate the ambient-RNA "soup" profile from raw drops outside the cell-
#     called set (low-count droplets).
#   - Quick-cluster the filtered cells (Seurat-style) so SoupX has cluster
#     labels for its profile-aware correction. We cluster with a deliberately
#     coarse Louvain at low resolution to keep runtime modest.
#   - autoEstCont() picks rho per cluster; adjustCounts() applies it.
suppressPackageStartupMessages({
    library(Matrix)
    library(SoupX)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript soupx.R <input_dir> <output_dir> [max_contamination]")
}
in_dir <- args[1]
out_dir <- args[2]
max_contam <- if (length(args) >= 3) as.numeric(args[3]) else 1.0

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

filtered <- readMM(file.path(in_dir, "filtered_counts.mtx"))    # cells x genes
filt_bc <- scan(file.path(in_dir, "filtered_barcodes.tsv"),
                                what = character(), quiet = TRUE, sep = "\n")
features <- scan(file.path(in_dir, "features.tsv"),
                                 what = character(), quiet = TRUE, sep = "\n")
raw <- readMM(file.path(in_dir, "raw_counts.mtx"))                  # droplets x genes
raw_bc <- scan(file.path(in_dir, "raw_barcodes.tsv"),
                             what = character(), quiet = TRUE, sep = "\n")

# SoupX expects genes x cells / genes x droplets.
tod <- t(raw); rownames(tod) <- features; colnames(tod) <- raw_bc
toc <- t(filtered); rownames(toc) <- features; colnames(toc) <- filt_bc

# A coarse cluster label is enough for autoEstCont().
# We use a quick log-CPM PCA + k-means as a SoupX-internal stand-in (pure-base R)
# to avoid pulling Seurat in as a hard dep.
quick_clusters <- function(counts_mat, k = 8) {
    libsize <- pmax(colSums(counts_mat), 1)
    norm <- log1p(t(t(counts_mat) / libsize) * 1e4)
    # Variance filter to keep PCA cheap.
    rv <- apply(norm, 1, function(x) var(x))
    top_idx <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
    norm <- norm[top_idx, , drop = FALSE]
    set.seed(0)
    pcs <- prcomp(t(as.matrix(norm)), center = TRUE, scale. = FALSE,
                                rank. = min(20, ncol(norm) - 1))
    set.seed(0)
    n_k <- min(k, max(2, ncol(counts_mat) %/% 30))
    km <- kmeans(pcs$x, centers = n_k, nstart = 10, iter.max = 50)
    setNames(as.integer(km$cluster), colnames(counts_mat))
}

clusters <- quick_clusters(toc)

sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)
sc <- setClusters(sc, clusters)
sc <- tryCatch(autoEstCont(sc, doPlot = FALSE),
                             error = function(e) {
                                 message(sprintf("SoupX autoEstCont failed (%s); falling back to setContaminationFraction(0.10)", e$message))
                                 setContaminationFraction(sc, 0.10, forceAccept = TRUE)
                             })

# Optional cap on per-cell rho.
if (!is.null(max_contam) && is.finite(max_contam) && max_contam < 1.0) {
    rho <- sc$metaData$rho
    rho_capped <- pmin(rho, max_contam)
    if (any(rho_capped < rho)) {
        sc$metaData$rho <- rho_capped
    }
}

decon_t <- adjustCounts(sc, roundToInt = TRUE)   # genes x cells
decon <- t(decon_t)                                                            # cells x genes
decon <- drop0(decon)

contam <- sc$metaData$rho
names(contam) <- filt_bc

writeMM(decon, file.path(out_dir, "decontaminated.mtx"))
write.table(contam, file.path(out_dir, "contamination.tsv"),
                        quote = FALSE, row.names = FALSE, col.names = FALSE)

summary_obj <- list(
    method = "SoupX",
    n_cells = ncol(toc),
    n_genes = nrow(toc),
    n_raw_droplets = ncol(tod),
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
    total_counts_pre = sum(toc),
    total_counts_post = sum(decon_t)
)
write_json(summary_obj, file.path(out_dir, "summary.json"),
                     auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("SoupX done: %d cells, median contamination=%.3f\n",
                        length(contam), median(contam)))
