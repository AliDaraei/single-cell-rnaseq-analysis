# 02_doublet_filtering_snRNAseq.R
# ------------------------------------------------------------
# Doublet-aware preprocessing for single-nucleus RNA-seq
#
# This script documents and applies doublet mitigation
# strategies depending on data provenance:
#
#   Strategy A: CellBender
#     - Used when starting from raw Cell Ranger outputs
#     - Models ambient RNA and background contamination
#     - Does NOT explicitly classify doublets
#
#   Strategy B: scds
#     - Applied to Seurat objects post-QC
#     - Scores potential doublets using gene co-expression
#
# The choice of strategy is dataset-dependent.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load QC-filtered object
# ---------------------------

input_path <- "results/seurat_qc_filtered.rds"

if (!file.exists(input_path)) {
  stop("QC-filtered Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

cat("Loaded QC-filtered object with",
    ncol(seu), "nuclei\n")

# ---------------------------
# 2. Strategy A: CellBender
# ---------------------------

# IMPORTANT:
# CellBender should be run BEFORE this script when raw
# Cell Ranger outputs are available.
#
# CellBender removes ambient RNA and technical background
# but does not explicitly label or score doublets.
#
# If CellBender was applied upstream:
#   - No additional action is required here
#   - Proceed directly to normalization

cellbender_applied <- FALSE  # set TRUE if applicable

if (cellbender_applied) {
  cat("CellBender preprocessing assumed; skipping scds\n")
}

# ---------------------------
# 3. Strategy B: scds (optional)
# ---------------------------

# scds is applied ONLY when:
# - Starting from an existing Seurat object
# - CellBender was not used
# - Additional doublet control is desired

if (!cellbender_applied) {

  cat("Running scds doublet scoring\n")

  suppressPackageStartupMessages({
    library(scds)
    library(SingleCellExperiment)
  })

  sce <- as.SingleCellExperiment(seu)

  # Co-expression based doublet score
  sce <- cxds(sce)

  # Binary classification based doublet score
  sce <- bcds(sce)

  # Combined score
  sce <- cxds_bcds_hybrid(sce)

  # Add scores to Seurat metadata
  seu$doublet_score <- colData(sce)$hybrid_score

  # Conservative thresholding
  cutoff <- quantile(seu$doublet_score, 0.95, na.rm = TRUE)

  seu <- subset(
    seu,
    subset = doublet_score < cutoff
  )

  cat("scds filtering applied\n")
  cat("Doublet score cutoff:", round(cutoff, 3), "\n")
  cat("Remaining nuclei:", ncol(seu), "\n")
}

# ---------------------------
# 4. Save post-doublet object
# ---------------------------

output_path <- "results/seurat_post_doublet.rds"
saveRDS(seu, file = output_path)

cat("Post-doublet Seurat object saved to:",
    output_path, "\n")
