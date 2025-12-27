# 14_scVI_latent_space_snRNAseq.py
# ------------------------------------------------------------
# scVI latent space modeling for single-nucleus RNA-seq
#
# This script trains an scVI model to learn a low-dimensional,
# denoised latent representation of snRNA-seq data.
#
# scVI is used here for:
#   - exploratory structure validation
#   - robustness checks
#   - downstream ML or visualization
#
# It does NOT replace Seurat/Harmony integration.
#
# Author: Ali Daraei
# Repository: single-cell-rnaseq-analysis
# ------------------------------------------------------------

import scanpy as sc
import scvi
import anndata as ad
import os

# ---------------------------
# 1. Input / output paths
# ---------------------------

INPUT_FILE = "results/seurat_for_scvi.h5ad"
OUTPUT_FILE = "results/scvi_latent_space.h5ad"

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(
        f"Input AnnData file not found: {INPUT_FILE}\n"
        "Export from Seurat using SeuratDisk before running scVI."
    )

# ---------------------------
# 2. Load data
# ---------------------------

adata = sc.read_h5ad(INPUT_FILE)

print(f"Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes")

# ---------------------------
# 3. Setup scVI
# ---------------------------

# Required metadata:
#   - batch_key (e.g., sample)
if "sample" not in adata.obs.columns:
    raise ValueError("Metadata column 'sample' is required for scVI")

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="sample"
)

# ---------------------------
# 4. Train model
# ---------------------------

model = scvi.model.SCVI(
    adata,
    n_latent=30
)

model.train(
    max_epochs=100,
    early_stopping=True
)

# ---------------------------
# 5. Extract latent space
# ---------------------------

adata.obsm["X_scVI"] = model.get_latent_representation()

print("Latent space learned:", adata.obsm["X_scVI"].shape)

# ---------------------------
# 6. Save output
# ---------------------------

adata.write(OUTPUT_FILE)

print(f"scVI latent space saved to: {OUTPUT_FILE}")
