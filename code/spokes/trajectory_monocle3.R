#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# trajectory_monocle3.R — trajectory & pseudotime (Monocle 3).   (Part 7)
# Needs GitHub install: cole-trapnell-lab/monocle3 (+ leidenbase).
# In: results/objects/05_annotated.rds   Out: results/objects/spoke_monocle3.rds
# Run: Rscript code/spokes/trajectory_monocle3.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(monocle3); library(SeuratWrappers) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

cds <- as.cell_data_set(obj)               # Seurat -> Monocle3 (via SeuratWrappers)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
# order_cells needs a root; set programmatically or via root_cells = <barcodes>:
# cds <- order_cells(cds, root_cells = colnames(cds)[obj$celltype == "<root type>"])
saveRDS(cds, out_path(cfg, "objects", "spoke_monocle3.rds"))
log_session(cfg, "spoke_monocle3")
# plot: plot_cells(cds, color_cells_by="pseudotime")
