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

# Auto-detect most abundant cell type (HumanPrimaryCellAtlasData labels vary).
ct_counts <- sort(table(obj$celltype), decreasing = TRUE)
target_ct <- names(ct_counts)[1]
message("hdWGCNA target cell type: ", target_ct, " (", ct_counts[1], " cells)")

# Subset to target cell type before the WGCNA pipeline: ModuleEigengenes
# densifies a target_ct × genes matrix on all single cells while the full
# Seurat object is resident, OOMing even at 90 GB. Subsetting first drops
# ~35 GB of non-target cells; MetacellsByGroups only uses target_ct metacells
# downstream anyway, so the result is identical.
obj <- subset(obj, subset = celltype == target_ct)

obj <- SetupForWGCNA(obj, gene_select = "fraction", fraction = 0.05, wgcna_name = "scrna")
obj <- MetacellsByGroups(obj, group.by = c("celltype", "condition"),
                         k = 25, ident.group = "celltype")
obj <- NormalizeMetacells(obj)
obj <- SetDatExpr(obj, group_name = target_ct, group.by = "celltype")
obj <- TestSoftPowers(obj)
# overwrite_tom = TRUE: re-runs don't block on a stale TOM file from a prior
# failed attempt (ConstructNetwork writes TOM to disk before ModuleEigengenes).
obj <- ConstructNetwork(obj, setDatExpr = FALSE, overwrite_tom = TRUE)
obj <- ModuleEigengenes(obj)
obj <- ModuleConnectivity(obj)
saveRDS(obj, out_path(cfg, "objects", "spoke_hdwgcna.rds"))
log_session(cfg, "spoke_hdwgcna")
# downstream: ModuleTraitCorrelation() vs condition to find disease-linked modules.
