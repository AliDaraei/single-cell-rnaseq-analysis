# 12_celltype_specific_pseudobulk_DE_snRNAseq.R
# ------------------------------------------------------------
# Cell-type–specific pseudobulk differential expression (snRNA-seq)
#
# This script performs pseudobulk DE separately for each
# annotated cell type using edgeR, ensuring biological
# replication is respected.
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

cell_types <- sort(unique(seu$cell_type))

dir.create("results/celltype_DE", recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 2. Loop over cell types
# ---------------------------

for (ct in cell_types) {

  message("Processing cell type: ", ct)

  seu_ct <- subset(seu, subset = cell_type == ct)

  # Ensure replication
  meta_ct <- seu_ct@meta.data %>%
    distinct(sample, condition)

  if (nrow(meta_ct) < 4 ||
      min(table(meta_ct$condition)) < 2) {
    message("  Skipping ", ct, " (insufficient replication)")
    next
  }

  # ---------------------------
  # 3. Pseudobulk aggregation
  # ---------------------------

  pb_counts <- AggregateExpression(
    seu_ct,
    group.by = "sample",
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA

  pb_meta <- meta_ct %>%
    distinct(sample, condition) %>%
    arrange(match(sample, colnames(pb_counts)))

  # ---------------------------
  # 4. edgeR DE
  # ---------------------------

  dge <- DGEList(counts = pb_counts)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ condition, data = pb_meta)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)

  qlf <- glmQLFTest(fit, coef = 2)

  de_res <- topTags(qlf, n = Inf)$table
  de_res$gene <- rownames(de_res)
  de_res$cell_type <- ct

  # ---------------------------
  # 5. Save results
  # ---------------------------

  out_file <- file.path(
    "results/celltype_DE",
    paste0("DE_", gsub(" ", "_", ct), ".csv")
  )

  write.csv(de_res, out_file, row.names = FALSE)

  message("  Results written to: ", out_file)
}

cat("Cell-type–specific pseudobulk DE completed\n")
