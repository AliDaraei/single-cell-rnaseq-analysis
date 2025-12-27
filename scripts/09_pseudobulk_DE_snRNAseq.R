# 09_pseudobulk_DE_snRNAseq.R
# ------------------------------------------------------------
# Pseudobulk differential expression analysis for snRNA-seq
#
# This script aggregates gene counts at the sample × cell_type
# level and performs differential expression using edgeR.
#
# Pseudobulk DE improves statistical robustness by respecting
# biological replication and reducing single-nucleus noise.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(edgeR)
})

# ---------------------------
# 1. Load annotated object
# ---------------------------

input_path <- "results/seurat_annotated.rds"

if (!file.exists(input_path)) {
  stop("Annotated Seurat object not found: ", input_path)
}

seu <- readRDS(input_path)

required_meta <- c("sample", "condition", "cell_type")
missing_meta <- setdiff(required_meta, colnames(seu@meta.data))

if (length(missing_meta) > 0) {
  stop("Missing required metadata: ",
       paste(missing_meta, collapse = ", "))
}

# ---------------------------
# 2. Pseudobulk aggregation
# ---------------------------

cat("Aggregating counts by sample and cell type\n")

pb_counts <- AggregateExpression(
  seu,
  group.by = c("sample", "cell_type"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)$RNA

# Build pseudobulk metadata
pb_meta <- expand.grid(
  sample    = unique(seu$sample),
  cell_type = unique(seu$cell_type)
) %>%
  filter(paste(sample, cell_type, sep = "_") %in% colnames(pb_counts))

pb_meta <- pb_meta %>%
  left_join(
    seu@meta.data %>% distinct(sample, condition),
    by = "sample"
  )

# ---------------------------
# 3. edgeR setup
# ---------------------------

dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge)

design <- model.matrix(~ condition, data = pb_meta)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# ---------------------------
# 4. Differential expression
# ---------------------------

qlf <- glmQLFTest(fit, coef = 2)

de_table <- topTags(qlf, n = Inf)$table
de_table$gene <- rownames(de_table)

# ---------------------------
# 5. Save results
# ---------------------------

dir.create("results", showWarnings = FALSE)

output_path <- "results/pseudobulk_DE_results.csv"

write.csv(
  de_table,
  file = output_path,
  row.names = FALSE
)

cat("Pseudobulk DE results saved to:",
    output_path, "\n")
