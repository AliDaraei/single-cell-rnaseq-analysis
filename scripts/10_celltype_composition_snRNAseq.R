# 10_celltype_composition_snRNAseq.R
# ------------------------------------------------------------
# Cell-type composition analysis for single-nucleus RNA-seq
#
# This script computes and visualizes cell-type proportions
# across samples and experimental conditions.
#
# Composition analysis helps distinguish changes in
# cell abundance from transcriptional regulation.
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
# 2. Compute composition table
# ---------------------------

comp_table <- seu@meta.data %>%
  count(sample, condition, cell_type) %>%
  group_by(sample) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()

dir.create("results", showWarnings = FALSE)

write.csv(
  comp_table,
  file = "results/cell_type_composition_table.csv",
  row.names = FALSE
)

cat("Cell-type composition table saved\n")

# ---------------------------
# 3. Visualization
# ---------------------------

p <- ggplot(
  comp_table,
  aes(
    x = sample,
    y = proportion,
    fill = cell_type
  )
) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~ condition, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Sample",
    y = "Cell-type proportion",
    fill = "Cell type",
    title = "Cell-type composition across samples"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

dir.create("figures", showWarnings = FALSE)

ggsave(
  filename = "figures/cell_type_composition_barplot.pdf",
  plot = p,
  width = 10,
  height = 6
)

cat("Cell-type composition plot saved to figures/\n")
