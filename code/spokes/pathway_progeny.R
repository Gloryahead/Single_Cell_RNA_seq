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

# Pseudobulk: mean log-norm expression per celltype×sample (~100 columns total).
# Per-cell run_mlm (105K cols) OOMs even at 60 GB — internal dense conversion.
# Pseudobulk is also the scientifically preferred PROGENy usage (Schubert 2018).
obj$pb_id <- paste(obj$celltype, obj$sample, sep = "..")
groups    <- split(seq_len(ncol(obj)), obj$pb_id)
pb_mat    <- do.call(cbind, lapply(groups, function(idx)
  Matrix::rowMeans(GetAssayData(obj, layer = "data")[, idx, drop = FALSE])))
colnames(pb_mat) <- names(groups)
rm(obj); gc()

acts <- decoupleR::run_mlm(mat = pb_mat, net = net,
                           .source = "source", .target = "target", .mor = "weight",
                           minsize = 5)
write.csv(acts, out_path(cfg, "tables", "spoke_pathway_activity.csv"), row.names = FALSE)
log_session(cfg, "spoke_progeny")
# downstream: compare pathway activity (e.g. TNFa/NFkB) across ctrl/pre/post.
