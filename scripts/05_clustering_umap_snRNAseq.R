# 05_clustering_umap_snRNAseq.R
# ------------------------------------------------------------
# Graph-based clustering and UMAP embedding for snRNA-seq
#
# This script operates on a Harmony-integrated Seurat object
# and performs:
#   - neighbor graph construction
#   - clustering at a defined resolution
#   - UMAP visualization
#
# Parameter choices should be tuned per dataset.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
})

# ---------------------------
# 1. Load integrated object
# ---------------------------

input_path <- "results/seurat_harmony_integrated.rds"

if (!file.exists(input_path)) {
  stop("Integrated Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

cat("Loaded integrated object with",
    ncol(seu), "nuclei\n")

# ---------------------------
# 2. Dimensionality selection
# ---------------------------

# Typical range for snRNA-seq:
#   20–40 Harmony dimensions
#
# Adjust based on variance explained and dataset complexity

dims_use <- 1:30

# ---------------------------
# 3. Construct neighbor graph
# ---------------------------

seu <- FindNeighbors(
  seu,
  reduction = "harmony",
  dims = dims_use,
  verbose = FALSE
)

cat("Neighbor graph constructed\n")

# ---------------------------
# 4. Clustering
# ---------------------------

# Resolution controls cluster granularity
# Conservative starting value recommended

cluster_resolution <- 0.6

seu <- FindClusters(
  seu,
  resolution = cluster_resolution,
  verbose = FALSE
)

cat("Clustering completed at resolution",
    cluster_resolution, "\n")

# ---------------------------
# 5. UMAP embedding
# ---------------------------

seu <- RunUMAP(
  seu,
  reduction = "harmony",
  dims = dims_use,
  verbose = FALSE
)

cat("UMAP embedding computed\n")

# ---------------------------
# 6. Save clustered object
# ---------------------------

output_path <- "results/seurat_harmony_clustered.rds"
saveRDS(seu, file = output_path)

cat("Clustered Seurat object saved to:",
    output_path, "\n")

