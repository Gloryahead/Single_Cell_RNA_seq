#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# ccc_cellchat.R — cell-cell communication (CellChat).            (Part 9)
# Needs GitHub install: jinworks/CellChat.
# In: results/objects/05_annotated.rds   Out: results/objects/spoke_cellchat.rds
# Run: Rscript code/spokes/ccc_cellchat.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(CellChat))
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# TIP: for disease-vs-control comparison, build one CellChat object PER condition
#      (subset obj by condition) and compare with mergeCellChat().
cc <- createCellChat(obj, group.by = "celltype", assay = "RNA")
cc@DB <- if (cfg$project$organism == "mouse") CellChatDB.mouse else CellChatDB.human
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc)
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
saveRDS(cc, out_path(cfg, "objects", "spoke_cellchat.rds"))
log_session(cfg, "spoke_cellchat")
# viz: netVisual_circle(cc@net$weight); netVisual_bubble(cc, ...)
