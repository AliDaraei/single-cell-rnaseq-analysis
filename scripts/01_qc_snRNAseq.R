# 01_qc_snRNAseq.R
# ------------------------------------------------------------
# Quality control for single-nucleus RNA-seq
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

cat("QC script initialized\n")
