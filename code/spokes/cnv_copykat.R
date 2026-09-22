#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# cnv_copykat.R — CNV inference to flag aneuploid cells (CopyKAT).  (Part 11)
# Needs GitHub install: navinlabcode/copykat.
# NOTE: built for tumor data. On PBMCs (all normal) it's a methods demo —
#       expect mostly "diploid". Most meaningful on a tumor-containing sample.
# In: results/objects/05_annotated.rds   Out: results/objects/spoke_copykat.rds
# Run: Rscript code/spokes/cnv_copykat.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(copykat))
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

counts <- as.matrix(GetAssayData(obj, layer = "counts"))
res <- copykat(rawmat = counts, id.type = "S", sam.name = cfg$project$name, n.cores = 8)
# norm.cell.names = <known normal barcodes> improves calls if available.
saveRDS(res, out_path(cfg, "objects", "spoke_copykat.rds"))
log_session(cfg, "spoke_copykat")
# res$prediction: aneuploid/diploid per cell -> map back to obj$cnv_status
