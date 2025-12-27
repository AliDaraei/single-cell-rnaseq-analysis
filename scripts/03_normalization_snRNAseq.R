# 03_normalization_snRNAseq.R
# ------------------------------------------------------------
# Normalization strategies for single-nucleus RNA-seq
#
# This script performs expression normalization and
# variable feature selection prior to integration.
#
# Two approaches are supported:
#   1) LogNormalize — lightweight, exploratory
#   2) SCTransform  — variance-stabilized, integration-ready
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
})

# ---------------------------
# 1. Load post-doublet object
# ---------------------------

input_path <- "results/seurat_post_doublet.rds"

if (!file.exists(input_path)) {
  stop("Post-doublet Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

cat("Loaded Seurat object with",
    ncol(seu), "nuclei\n")

# ---------------------------
# 2. Choose normalization method
# ---------------------------

# Options:
#   "log" → LogNormalize
#   "sct" → SCTransform
#
# Change this value depending on analysis goals

normalization_method <- "log"

# ---------------------------
# 3A. LogNormalize workflow
# ---------------------------

if (normalization_method == "log") {

  cat("Running LogNormalize\n")

  seu <- NormalizeData(
    seu,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )

  output_path <- "results/seurat_log_normalized.rds"
}

# ---------------------------
# 3B. SCTransform workflow
# ---------------------------

if (normalization_method == "sct") {

  cat("Running SCTransform\n")

  # For snRNA-seq, mitochondrial regression is optional
  # and often omitted unless strong technical effects exist

  seu <- SCTransform(
    seu,
    vars.to.regress = NULL,
    verbose = FALSE
  )

  output_path <- "results/seurat_sct_normalized.rds"
}

# ---------------------------
# 4. Save normalized object
# ---------------------------

saveRDS(seu, file = output_path)

cat("Normalized Seurat object saved to:",
    output_path, "\n")
