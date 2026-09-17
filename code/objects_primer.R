#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# objects_primer.R — understand & convert the object (the universal currency).
#                                                       (Tutorial Part 6)
# Not part of the pipeline flow — a reference you run to learn the data model.
# Inputs  : results/objects/05_annotated.rds
# Run     : Rscript code/objects_primer.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
cfg <- load_config()
obj <- readRDS(obj_path(cfg, "05_annotated"))

# ---- Seurat anatomy ----
print(obj)                                   # assays, cells, features
cat("\nLayers:\n"); print(Layers(obj))        # counts / data / scale.data
cat("\nmeta.data columns (where YOUR biology lives):\n")
print(colnames(obj[[]]))
cat("\nReductions:\n"); print(Reductions(obj)) # pca / harmony / umap

# accessors you'll use constantly:
counts <- GetAssayData(obj, layer = "counts") # raw matrix (genes x cells)
norm   <- GetAssayData(obj, layer = "data")   # log-normalized
emb    <- Embeddings(obj, "umap")             # cell coordinates

# ---- Seurat <-> SingleCellExperiment (R) ----
suppressPackageStartupMessages(library(SingleCellExperiment))
sce <- as.SingleCellExperiment(obj)           # for Bioconductor spokes (SingleR, slingshot)
obj2 <- as.Seurat(sce)                         # and back

# ---- Seurat -> AnnData (.h5ad, for the Python spokes) ----
# done in 05_annotate.R via export_h5ad(); shown here for reference:
# sceasy::convertFormat(obj, from="seurat", to="anndata",
#                       outFile=out_path(cfg,"objects","05_annotated.h5ad"))

cat("\nKey idea: Seurat (R), SingleCellExperiment (R/Bioc), AnnData (Python)\n",
    "hold the SAME data. Every spoke just needs the format its tool expects.\n")
