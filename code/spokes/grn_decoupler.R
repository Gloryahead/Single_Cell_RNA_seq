#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# grn_decoupler.R — transcription-factor activity (decoupleR + CollecTRI). (Part 16)
# Packages on conda: bioconductor-decoupler, bioconductor-omnipathr.
# In: results/objects/05_annotated.rds   Out: results/tables/spoke_tf_activity.csv
# Run: Rscript code/spokes/grn_decoupler.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(decoupleR); library(msigdbr) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# OmnipathR::unnest_evidences crashes on both collectri and dorothea static tables
# (version bug in the container). Use GTRD TF target gene sets from MSigDB C3 instead —
# no OmniPath dependency, msigdbr is confirmed installed.
sp  <- if (cfg$project$organism == "mouse") "Mus musculus" else "Homo sapiens"
tf_df <- msigdbr(species = sp, category = "C3", subcategory = "TFT:GTRD")
net   <- data.frame(source = tf_df$gs_name, target = tf_df$gene_symbol,
                    mor = 1L, stringsAsFactors = FALSE)

# Pseudobulk: mean log-norm expression per celltype×sample (~100 columns total).
# Per-cell run_ulm (105K cols) would OOM; pseudobulk is also more interpretable.
obj$pb_id <- paste(obj$celltype, obj$sample, sep = "..")
groups    <- split(seq_len(ncol(obj)), obj$pb_id)
pb_mat    <- do.call(cbind, lapply(groups, function(idx)
  Matrix::rowMeans(GetAssayData(obj, layer = "data")[, idx, drop = FALSE])))
colnames(pb_mat) <- names(groups)
rm(obj); gc()

acts <- decoupleR::run_ulm(mat = pb_mat, net = net,
                           .source = "source", .target = "target", .mor = "mor",
                           minsize = 5)
write.csv(acts, out_path(cfg, "tables", "spoke_tf_activity.csv"), row.names = FALSE)
log_session(cfg, "spoke_decoupler_tf")
# downstream: average TF activity per cell type / condition to find active regulators.
