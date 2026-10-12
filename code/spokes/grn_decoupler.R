#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# grn_decoupler.R — transcription-factor activity (decoupleR + CollecTRI). (Part 16)
# Packages on conda: bioconductor-decoupler, bioconductor-omnipathr.
# In: results/objects/05_annotated.rds   Out: results/tables/spoke_tf_activity.csv
# Run: Rscript code/spokes/grn_decoupler.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(decoupleR) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# CollecTRI via get_collectri() crashes (OmnipathR unnest_evidences bug on static table).
# DoRothEA uses a different OmnipathR code path — well-tested, no unnest_evidences call.
net <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"))
mat <- GetAssayData(obj, layer = "data")   # keep sparse — as.matrix() would materialize ~28 GB
rm(obj); gc()   # free Seurat object (~20 GB) before run_ulm
acts <- decoupleR::run_ulm(mat = mat, net = net,
                           .source = "tf", .target = "target", .mor = "mor",
                           minsize = 5)
write.csv(acts, out_path(cfg, "tables", "spoke_tf_activity.csv"), row.names = FALSE)
log_session(cfg, "spoke_decoupler_tf")
# downstream: average TF activity per cell type / condition to find active regulators.
