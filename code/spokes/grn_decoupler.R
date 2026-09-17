#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# grn_decoupler.R — transcription-factor activity (decoupleR + CollecTRI). (Part 16)
# Packages on conda: bioconductor-decoupler, bioconductor-omnipathr.
# In: results/objects/05_annotated.rds   Out: results/tables/spoke_tf_activity.csv
# Run: Rscript code/spokes/grn_decoupler.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(decoupleR); library(OmnipathR) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)  # TF->target network
mat <- as.matrix(GetAssayData(obj, layer = "data"))                            # log-norm expression
acts <- decoupleR::run_ulm(mat = mat, net = net,
                           .source = "source", .target = "target", .mor = "mor",
                           minsize = 5)
write.csv(acts, out_path(cfg, "tables", "spoke_tf_activity.csv"), row.names = FALSE)
log_session(cfg, "spoke_decoupler_tf")
# downstream: average TF activity per cell type / condition to find active regulators.
