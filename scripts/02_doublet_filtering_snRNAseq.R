# 02_doublet_filtering_snRNAseq.R
# ------------------------------------------------------------
# Doublet mitigation strategies for snRNA-seq data
#
# This script documents and optionally applies
# doublet-aware preprocessing approaches depending
# on dataset provenance.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
})

cat("Doublet filtering module initialized\n")
