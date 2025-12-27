# 13_umap_plots_snRNAseq.R
# ------------------------------------------------------------
# UMAP visualization utilities for snRNA-seq analyses
#
# This script generates publication-ready UMAP figures
# from clustered and annotated Seurat objects.
#
# Outputs include:
#   - UMAP by cluster
#   - UMAP by cell type
#   - UMAP by condition
#   - Optional hepatocyte subcluster UMAP
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ---------------------------
# 1. Load clustered / annotated object
# ---------------------------

input_path <- "results/seurat_annotated.rds"

if (!file.exists(input_path)) {
  stop("Annotated Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

dir.create("figures", showWarnings = FALSE)

# ---------------------------
# 2. Global UMAP: clusters
# ---------------------------

p_cluster <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("snRNA-seq UMAP — Clusters") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "figures/UMAP_clusters.png",
  plot = p_cluster,
  width = 7,
  height = 6,
  dpi = 300
)

# ---------------------------
# 3. Global UMAP: cell types
# ---------------------------

if ("cell_type" %in% colnames(seu@meta.data)) {

  p_celltype <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("snRNA-seq UMAP — Cell types") +
    theme_minimal(base_size = 14)

  ggsave(
    filename = "figures/UMAP_cell_types.png",
    plot = p_celltype,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# ---------------------------
# 4. Global UMAP: condition
# ---------------------------

if ("condition" %in% colnames(seu@meta.data)) {

  p_condition <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "condition"
  ) +
    ggtitle("snRNA-seq UMAP — Condition") +
    theme_minimal(base_size = 14)

  ggsave(
    filename = "figures/UMAP_condition.png",
    plot = p_condition,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# ---------------------------
# 5. Optional: hepatocyte subcluster UMAP
# ---------------------------

hep_path <- "results/hepatocyte_subcluster_annotated.rds"

if (file.exists(hep_path)) {

  hep <- readRDS(hep_path)

  if ("hep_subtype" %in% colnames(hep@meta.data)) {

    p_hep <- DimPlot(
      hep,
      reduction = "umap_hepatocyte",
      group.by = "hep_subtype",
      label = TRUE,
      repel = TRUE
    ) +
      ggtitle("Hepatocyte subclusters") +
      theme_minimal(base_size = 14)

    ggsave(
      filename = "figures/UMAP_hepatocyte_subclusters.png",
      plot = p_hep,
      width = 7,
      height = 6,
      dpi = 300
    )
  }
}

cat("UMAP figures generated in figures/\n")
