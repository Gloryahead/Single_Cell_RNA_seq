#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 04_cluster_integrate.R — filtered -> integrated, clustered object.
#                                                       (Tutorial Part 3)
# Inputs  : results/objects/03_filtered.rds
# Outputs : results/objects/04_clustered.rds + UMAPs in results/plots/
# Run     : Rscript code/04_cluster_integrate.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
cfg <- load_config(); set.seed(cfg$seed)
cc  <- cfg$cluster

obj <- readRDS(obj_path(cfg, "03_filtered"))

# ---- standard normalization -> HVG -> scale -> PCA ----
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = cc$n_hvg)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = cc$n_pcs)

# ---- batch integration (Part 3): only correct nuisance batch, never condition ----
reduction <- "pca"
if (cc$integration == "harmony") {
  suppressPackageStartupMessages(library(harmony))
  obj <- RunHarmony(obj, group.by.vars = cc$batch_key)
  reduction <- "harmony"
}

# ---- neighbors -> clusters -> UMAP, all on the integrated reduction ----
dims <- 1:cc$n_pcs
obj <- FindNeighbors(obj, reduction = reduction, dims = dims)
obj <- FindClusters(obj, resolution = cc$resolution)   # resolution is a knob, not truth
obj <- RunUMAP(obj, reduction = reduction, dims = dims)

saveRDS(obj, obj_path(cfg, "04_clustered"))
# sanity check: cells should mix by SAMPLE but separate by biology
ggsave(out_path(cfg, "plots", "04_umap_clusters_samples.png"),
       DimPlot(obj, label = TRUE) + DimPlot(obj, group.by = "sample"),
       width = 13, height = 5)
log_session(cfg, "04_cluster_integrate")
