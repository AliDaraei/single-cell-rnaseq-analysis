# 08_hepatocyte_subcluster_annotation_snRNAseq.R
# ------------------------------------------------------------
# Annotation of hepatocyte subclusters in snRNA-seq
#
# This script assigns biologically meaningful labels to
# hepatocyte subclusters based on known zonation,
# metabolic, and stress-associated marker genes.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ---------------------------
# 1. Load hepatocyte subclustered object
# ---------------------------

input_path <- "results/hepatocyte_subclustered.rds"

if (!file.exists(input_path)) {
  stop("Hepatocyte subclustered object not found: ", input_path)
}

hep <- readRDS(input_path)

cat("Loaded hepatocyte object with",
    length(unique(hep$seurat_clusters)),
    "subclusters\n")

# ---------------------------
# 2. Reference marker sets
# ---------------------------

# These markers reflect well-established liver zonation
# and stress-associated hepatocyte programs.

marker_reference <- list(
  Periportal = c("Cps1", "Ass1", "Hal"),
  Midzonal   = c("Hsd17b13", "Gstm1"),
  Pericentral = c("Glul", "Cyp2e1", "Cyp1a2"),
  Stress     = c("Jun", "Fos", "Atf3")
)

# ---------------------------
# 3. Visual validation
# ---------------------------

if (interactive()) {
  DotPlot(
    hep,
    features = marker_reference
  ) +
    RotatedAxis() +
    ggtitle("Hepatocyte subcluster marker expression")
}

# ---------------------------
# 4. Assign subcluster identities
# ---------------------------

# IMPORTANT:
# This mapping must be updated after inspecting markers.
# The example below is a biologically reasonable placeholder.

hep$hep_subtype <- case_when(
  hep$seurat_clusters %in% c("0", "2") ~ "Periportal-like",
  hep$seurat_clusters %in% c("1")      ~ "Midzonal-like",
  hep$seurat_clusters %in% c("3")      ~ "Pericentral-like",
  hep$seurat_clusters %in% c("4")      ~ "Stress-associated",
  TRUE                                 ~ "Unassigned"
)

hep$hep_subtype <- factor(
  hep$hep_subtype,
  levels = c(
    "Periportal-like",
    "Midzonal-like",
    "Pericentral-like",
    "Stress-associated",
    "Unassigned"
  )
)

cat("Hepatocyte subcluster labels assigned\n")

# ---------------------------
# 5. Save annotated object
# ---------------------------

output_path <- "results/hepatocyte_subcluster_annotated.rds"
saveRDS(hep, file = output_path)

cat("Annotated hepatocyte subclusters saved to:",
    output_path, "\n")
