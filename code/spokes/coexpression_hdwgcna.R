#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# coexpression_hdwgcna.R — gene co-expression networks (hdWGCNA).  (Part 12)
# Needs GitHub install: smorabit/hdWGCNA (ref="dev").
# In: results/objects/05_annotated.rds   Out: results/objects/spoke_hdwgcna.rds
# Run: Rscript code/spokes/coexpression_hdwgcna.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(hdWGCNA); library(WGCNA) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

obj <- SetupForWGCNA(obj, gene_select = "fraction", fraction = 0.05, wgcna_name = "scrna")
# metacells reduce sparsity; group by the cell types you want networks within:
obj <- MetacellsByGroups(obj, group.by = c("celltype", "condition"),
                         k = 25, ident.group = "celltype")
obj <- NormalizeMetacells(obj)
obj <- SetDatExpr(obj, group_name = "CD14 Mono", group.by = "celltype")  # TODO: target cell type
obj <- TestSoftPowers(obj)
obj <- ConstructNetwork(obj, setDatExpr = FALSE)
obj <- ModuleEigengenes(obj)
obj <- ModuleConnectivity(obj)
saveRDS(obj, out_path(cfg, "objects", "spoke_hdwgcna.rds"))
log_session(cfg, "spoke_hdwgcna")
# downstream: ModuleTraitCorrelation() vs condition to find disease-linked modules.
