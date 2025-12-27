# 07_hepatocyte_subclustering_snRNAseq.R
# ------------------------------------------------------------
# Cell-type–specific subclustering of hepatocytes in snRNA-seq
#
# This script refines hepatocyte heterogeneity that is often
# masked during global clustering by:
#   - subsetting hepatocytes
#   - re-normalizing expression
#   - conservative reclustering
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load annotated object
# ---------------------------

input_path <- "results/seurat_annotated.rds"

if (!file.exists(input_path)) {
  stop("Annotated Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

if (!"cell_type" %in% colnames(seu@meta.data)) {
  stop("Required metadata column 'cell_type' is missing")
}

# ---------------------------
# 2. Subset hepatocytes
# ---------------------------

hep <- subset(seu, subset = cell_type == "Hepatocyte")

cat("Hepatocyte nuclei retained:",
    ncol(hep), "\n")

if (ncol(hep) < 100) {
  warning("Low number of hepatocyte nuclei — subclustering may be unstable")
}

# ---------------------------
# 3. Re-normalization
# ---------------------------

# Re-normalization is REQUIRED for valid subclustering
# to remove global variance structure

hep <- NormalizeData(
  hep,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

hep <- FindVariableFeatures(
  hep,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

hep <- ScaleData(
  hep,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)

# ---------------------------
# 4. PCA
# ---------------------------

hep <- RunPCA(
  hep,
  npcs = 50,
  verbose = FALSE
)

if (interactive()) {
  ElbowPlot(hep, ndims = 50)
}

# Conservative dimensionality for subclustering
pcs_use <- 1:25

# ---------------------------
# 5. Subclustering
# ---------------------------

hep <- FindNeighbors(
  hep,
  dims = pcs_use,
  verbose = FALSE
)

hep <- FindClusters(
  hep,
  resolution = 0.3,
  verbose = FALSE
)

hep <- RunUMAP(
  hep,
  dims = pcs_use,
  reduction.name = "umap_hepatocyte",
  verbose = FALSE
)

cat("Hepatocyte subclustering completed\n")

# ---------------------------
# 6. Marker discovery
# ---------------------------

hep_markers <- FindAllMarkers(
  hep,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  hep_markers,
  file = "results/hepatocyte_subcluster_markers.csv",
  row.names = FALSE
)

cat("Hepatocyte subcluster markers saved\n")

# ---------------------------
# 7. Save hepatocyte object
# ---------------------------

output_path <- "results/hepatocyte_subclustered.rds"
saveRDS(hep, file = output_path)

cat("Hepatocyte subclustered object saved to:",
    output_path, "\n")
