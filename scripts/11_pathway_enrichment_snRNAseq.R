# 11_pathway_enrichment_snRNAseq.R
# ------------------------------------------------------------
# Pathway-level analysis for snRNA-seq pseudobulk results
#
# This script performs:
#   1) Gene Set Enrichment Analysis (GSEA)
#   2) Over-representation analysis (ORA)
#
# GSEA is applied to ranked genes to avoid arbitrary cutoffs,
# while ORA summarizes pathways enriched among significant genes.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# ---------------------------
# 1. Load pseudobulk DE results
# ---------------------------

de_path <- "results/pseudobulk_DE_results.csv"

if (!file.exists(de_path)) {
  stop("Pseudobulk DE results not found: ", de_path)
}

de <- read.csv(de_path, stringsAsFactors = FALSE)

if (!all(c("gene", "logFC") %in% colnames(de))) {
  stop("Required columns missing in DE table (gene, logFC)")
}

# ---------------------------
# 2. Prepare ranked gene list
# ---------------------------

ranked_genes <- de %>%
  filter(!is.na(logFC)) %>%
  arrange(desc(logFC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(gene, logFC)

gene_ranks <- ranked_genes$logFC
names(gene_ranks) <- ranked_genes$gene

# ---------------------------
# 3. GSEA using Hallmark gene sets
# ---------------------------

hallmark <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

hallmark_list <- split(
  hallmark$gene_symbol,
  hallmark$gs_name
)

gsea_res <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark_list,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

write.csv(
  as.data.frame(gsea_res),
  file = "results/GSEA_hallmark_results.csv",
  row.names = FALSE
)

pdf(
  file = "figures/GSEA_hallmark_dotplot.pdf",
  width = 8,
  height = 6
)
dotplot(gsea_res, showCategory = 20) +
  ggtitle("Hallmark GSEA (pseudobulk)")
dev.off()

cat("GSEA analysis completed\n")

# ---------------------------
# 4. ORA (GO Biological Process)
# ---------------------------

sig_genes <- de %>%
  filter(FDR < 0.05, abs(logFC) > 0.5) %>%
  pull(gene) %>%
  unique()

if (length(sig_genes) < 10) {
  warning("Few significant genes detected — ORA may be underpowered")
}

entrez <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
  gene          = entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(
  as.data.frame(ego),
  file = "results/ORA_GO_results.csv",
  row.names = FALSE
)

p <- dotplot(ego, showCategory = 20) +
  ggtitle("GO Biological Process ORA")

ggsave(
  filename = "figures/ORA_GO_dotplot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("ORA analysis completed\n")
