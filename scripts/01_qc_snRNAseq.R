# 01_qc_snRNAseq.R
# ------------------------------------------------------------
# Quality control for single-nucleus RNA-seq (snRNA-seq)
#
# This script performs initial quality assessment and
# conservative filtering of nuclei prior to downstream
# integration and clustering.
#
# QC metrics are intentionally lightweight to accommodate
# nuclear RNA profiles, which differ from whole-cell RNA-seq.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load input Seurat object
# ---------------------------

# Expected input:
# - Raw or minimally processed Seurat object
# - Contains RNA assay with raw counts
# - Metadata may already include sample information

input_path <- "data/seurat_raw_example.rds"

if (!file.exists(input_path)) {
  stop("Input Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

cat("Loaded Seurat object with",
    ncol(seu), "nuclei and",
    nrow(seu), "features\n")

# ---------------------------
# 2. Compute QC metrics
# ---------------------------

# Mitochondrial percentage
# For snRNA-seq, this metric is informative but not
# always used as a strict filtering criterion

seu[["percent.mt"]] <- PercentageFeatureSet(
  seu,
  pattern = "^mt-"
)

qc_summary <- seu@meta.data %>%
  select(nFeature_RNA, nCount_RNA, percent.mt)

print(summary(qc_summary))

# ---------------------------
# 3. Visual QC inspection
# ---------------------------

# These plots are useful when running interactively
# They are not saved automatically

if (interactive()) {
  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )
}

# ---------------------------
# 4. Filter low-quality nuclei
# ---------------------------

# Thresholds should be adjusted per dataset.
# These values represent conservative starting points
# commonly used in snRNA-seq studies.

min_features <- 200
max_features <- 6000
max_mito     <- 5

seu_qc <- subset(
  seu,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_mito
)

cat("QC filtering complete\n")
cat("Nuclei before filtering:", ncol(seu), "\n")
cat("Nuclei after filtering :", ncol(seu_qc), "\n")

# ---------------------------
# 5. Save QC-filtered object
# ---------------------------

output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

output_path <- file.path(
  output_dir,
  "seurat_qc_filtered.rds"
)

saveRDS(seu_qc, file = output_path)

cat("QC-filtered Seurat object saved to:", output_path, "\n")
