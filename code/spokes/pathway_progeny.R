#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# pathway_progeny.R — signaling-pathway activity (decoupleR + PROGENy). (Part 17)
# Packages on conda: bioconductor-decoupler, bioconductor-progeny.
# In: results/objects/05_annotated.rds   Out: results/tables/spoke_pathway_activity.csv
# Run: Rscript code/spokes/pathway_progeny.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(decoupleR))
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

net <- decoupleR::get_progeny(organism = "human", top = 500)   # 14 pathway signatures
mat <- as.matrix(GetAssayData(obj, layer = "data"))
acts <- decoupleR::run_mlm(mat = mat, net = net,
                           .source = "source", .target = "target", .mor = "weight",
                           minsize = 5)
write.csv(acts, out_path(cfg, "tables", "spoke_pathway_activity.csv"), row.names = FALSE)
log_session(cfg, "spoke_progeny")
# downstream: compare pathway activity (e.g. TNFa/NFkB) across ctrl/pre/post.
