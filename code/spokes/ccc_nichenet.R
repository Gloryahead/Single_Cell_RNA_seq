#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# ccc_nichenet.R — ligand-target prediction (NicheNet).          (Part 10)
# Needs GitHub install: saeyslab/nichenetr  + prior-model files from Zenodo.
# In: results/objects/05_annotated.rds   Out: results/tables/spoke_nichenet_ligands.csv
# Run: Rscript code/spokes/ccc_nichenet.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(nichenetr); library(tidyverse) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# Download once from Zenodo (https://zenodo.org/record/7074291) and point here:
lr_network    <- readRDS("references/nichenet/lr_network_human_21122021.rds")
ligand_target <- readRDS("references/nichenet/ligand_target_matrix_nsga2r_final.rds")

# Define sender/receiver cell types and the condition contrast of interest:
res <- nichenet_seuratobj_aggregate(
  seurat_obj      = obj,
  receiver        = "CD8 T",                        # TODO: your receiver cell type
  sender          = c("CD14 Mono", "FCGR3A Mono"),  # TODO: your sender cell type(s)
  condition_colname = "condition", condition_oi = "pre", condition_reference = "ctrl",
  ligand_target_matrix = ligand_target, lr_network = lr_network,
  weighted_networks = NULL)
write.csv(res$ligand_activities,
          out_path(cfg, "tables", "spoke_nichenet_ligands.csv"), row.names = FALSE)
log_session(cfg, "spoke_nichenet")
