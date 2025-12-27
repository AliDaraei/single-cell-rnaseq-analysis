# 04_integration_harmony_snRNAseq.R
# ------------------------------------------------------------
# Dataset integration and batch correction for snRNA-seq
#
# This script performs:
#   - scaling on variable genes
#   - PCA for dimensionality reduction
#   - batch correction using Harmony
#
# Harmony is used to align multiple samples or batches
# while preserving biological structure.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

# ---------------------------
# 1. Load normalized object
# ---------------------------

# Use log-normalized object by default
input_path <- "results/seurat_log_normalized.rds"

if (!file.exists(input_path)) {
  stop("Normalized Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

cat("Loaded normalized object with",
    ncol(seu), "nuclei\n")

# ---------------------------
# 2. Metadata requirements
# ---------------------------

# Batch variable required for Harmony
# Common examples:
#   - sample
#   - mouse_id
#   - sequencing_run

batch_var <- "sample"

if (!batch_var %in% colnames(seu@meta.data)) {
  stop("Required metadata column missing: ", batch_var)
}

# ---------------------------
# 3. Scaling and PCA
# ---------------------------

# Scale only variable features to reduce memory usage
seu <- ScaleData(
  seu,
  features = VariableFeatures(seu),
  verbose = FALSE
)

seu <- RunPCA(
  seu,
  features = VariableFeatures(seu),
  npcs = 50,
  verbose = FALSE
)

cat("PCA completed\n")

# ---------------------------
# 4. Harmony integration
# ---------------------------

seu <- RunHarmony(
  object = seu,
  group.by.vars = batch_var,
  reduction = "pca",
  assay.use = DefaultAssay(seu),
  verbose = TRUE
)

cat("Harmony integration completed\n")

# ---------------------------
# 5. Save integrated object
# ---------------------------

output_path <- "results/seurat_harmony_integrated.rds"
saveRDS(seu, file = output_path)

cat("Harmony-integrated object saved to:",
    output_path, "\n")
