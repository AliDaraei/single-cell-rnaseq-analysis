# single-cell-rnaseq-analysis

This repository contains a modular, research-grade workflow for analyzing
single-nucleus RNA-sequencing (snRNA-seq) data.

The pipeline reflects real analytical practices used in biological research
and is designed to prioritize:

- biological interpretability  
- statistical rigor  
- reproducibility  
- clean, readable code  

The workflows are written as sequential scripts rather than a single
monolithic package, making each analytical decision transparent and easy
to adapt to new datasets.

---

## 🔬 Scope of the repository

This repository covers the full snRNA-seq analysis lifecycle:

- quality control and filtering  
- doublet handling and ambient RNA considerations  
- normalization and batch correction  
- clustering and cell type annotation  
- cell-type–specific subclustering  
- pseudobulk differential expression  
- pathway enrichment  
- cell-type composition analysis  
- visualization utilities  
- optional latent-space modeling with scVI  

Raw sequencing data are **not included**.

---

## 📁 Repository structure

single-cell-rnaseq-analysis/
├── scripts/
│ ├── 01_qc_snRNAseq.R
│ ├── 02_doublet_filtering_snRNAseq.R
│ ├── 03_normalization_snRNAseq.R
│ ├── 04_integration_harmony_snRNAseq.R
│ ├── 05_clustering_umap_snRNAseq.R
│ ├── 06_annotation_snRNAseq.R
│ ├── 07_hepatocyte_subclustering_snRNAseq.R
│ ├── 08_hepatocyte_subcluster_annotation_snRNAseq.R
│ ├── 09_pseudobulk_DE_snRNAseq.R
│ ├── 10_celltype_composition_snRNAseq.R
│ ├── 11_pathway_enrichment_snRNAseq.R
│ ├── 12_celltype_specific_pseudobulk_DE_snRNAseq.R
│ ├── 13_umap_plots_snRNAseq.R
│ └── 14_scVI_latent_space_snRNAseq.py
├── results/ # Generated result objects and tables (not tracked)
├── figures/ # Generated figures (PNG / PDF, not tracked)
├── envs/ # Environment documentation
├── README.md
└── METHODS.md

---

## 🔁 Analysis workflow (high level)

Quality control  
→ Doublet handling  
→ Normalization  
→ Batch correction (Harmony)  
→ Clustering and UMAP  
→ Cell type annotation  
→ Cell-type–specific subclustering  
→ Pseudobulk differential expression  
→ Pathway enrichment  
→ Visualization and interpretation  

---

## 🧬 Analysis philosophy

### Why pseudobulk?
Differential expression is performed using pseudobulk aggregation to respect
biological replication and avoid statistical pitfalls associated with
single-cell–level testing.

### Why cell-type–specific subclustering?
Global clustering can mask biologically meaningful heterogeneity, especially
within abundant populations such as hepatocytes. Subclustering enables
resolution of zonation, stress responses, and disease-associated programs.

### Why scVI?
scVI is included as an **optional exploratory module** to validate
transcriptional structure in a probabilistic latent space. It is not used
for annotation or differential expression.

---

## ⚙️ Requirements

### R (≥ 4.2)
- Seurat  
- harmony  
- edgeR  
- dplyr  
- ggplot2  
- clusterProfiler  
- msigdbr  
- ComplexHeatmap  

### Python (≥ 3.9)
- scanpy  
- scvi-tools  
- anndata  
- PyTorch backend  

Exact versions may vary by system and dataset.

---

## ♻️ Reproducibility notes

- Scripts are designed to be run sequentially  
- Parameters (QC thresholds, clustering resolution) are dataset-specific  
- Manual biological interpretation is required for annotation steps  
- Intermediate objects are saved to support inspection and reuse  

---

## 🎯 Intended audience

This repository is intended for:

- computational biologists  
- bioinformaticians  
- researchers working with scRNA-seq or snRNA-seq  
- scientists interested in reproducible single-cell workflows  

---

## 👤 Author

**Ali Daraei**  
Computational biology · Single-cell transcriptomics
