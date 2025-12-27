# 06_celltype_annotation_snRNAseq.R
# ------------------------------------------------------------
# Marker-based cell type annotation for snRNA-seq
#
# This script assigns biological labels to clusters using
# canonical marker genes and manual inspection.
#
# Annotations are intentionally conservative and should be
# revisited as additional data or validation becomes available.
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
# 1. Load clustered object
# ---------------------------

input_path <- "results/seurat_harmony_clustered.rds"

if (!file.exists(input_path)) {
  stop("Clustered Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

cat("Loaded clustered object with",
    length(unique(seu$seurat_clusters)), "clusters\n")

# ---------------------------
# 2. Identify cluster markers
# ---------------------------

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

dir.create("results", showWarnings = FALSE)

write.csv(
  markers,
  file = "results/cluster_markers_all.csv",
  row.names = FALSE
)

cat("Cluster marker table written to results/\n")

# ---------------------------
# 3. Define canonical markers
# ---------------------------

# NOTE:
# Marker sets should be adapted to tissue and species.
# These are example markers commonly used in liver snRNA-seq.

marker_reference <- list(
  Hepatocyte     = c("Alb", "Ttr", "Apoa1"),
  Kupffer        = c("Lyz2", "Clec4f", "Adgre1"),
  LSEC           = c("Kdr", "Pecam1", "Rspo3"),
  Stellate       = c("Col1a1", "Des", "Lrat"),
  Cholangiocyte  = c("Krt19", "Krt7"),
  T_cell         = c("Cd3d", "Cd3e"),
  B_cell         = c("Cd79a", "Ms4a1")
)

# ---------------------------
# 4. Marker visualization
# ---------------------------

if (interactive()) {
  DotPlot(
    seu,
    features = marker_reference
  ) +
    RotatedAxis() +
    ggtitle("Canonical marker expression by cluster")
}

# ---------------------------
# 5. Assign cell-type labels
# ---------------------------

# IMPORTANT:
# Update this mapping after reviewing marker expression.
# This example mapping is a placeholder.

cluster_to_celltype <- c(
  "0" = "Hepatocyte",
  "1" = "Kupffer",
  "2" = "LSEC",
  "3" = "Stellate",
  "4" = "T_cell",
  "5" = "B_cell"
)

seu$cell_type <- dplyr::recode(
  seu$seurat_clusters,
  !!!cluster_to_celltype,
  .default = "Unassigned"
)

seu$cell_type <- factor(seu$cell_type)

cat("Cell-type annotation added to metadata\n")

# ---------------------------
# 6. Save annotated object
# ---------------------------

output_path <- "results/seurat_annotated.rds"
saveRDS(seu, file = output_path)

cat("Annotated Seurat object saved to:",
    output_path, "\n")
