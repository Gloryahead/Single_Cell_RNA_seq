#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# trajectory_slingshot.R — trajectory & pseudotime (Slingshot). (Part 7-2)
# Pattern: load hub -> run tool -> save. (Every spoke is this shape.)
# In: results/objects/05_annotated.rds   Out: results/objects/spoke_slingshot.rds
# Run: Rscript code/spokes/trajectory_slingshot.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({ library(slingshot); library(SingleCellExperiment) })
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# OPTIONAL: subset to the compartment with a real trajectory (e.g. a lineage).
# obj <- subset(obj, subset = celltype %in% c("Progenitors","Intermediate","Mature"))

# Drop NA and singleton clusters — slingshot's MST fails when a centroid has NA
# coordinates (happens with clusters of < 2 cells or unexpected NA labels).
valid_ct <- !is.na(obj$celltype)
ct_counts <- table(obj$celltype[valid_ct])
keep_ct   <- names(ct_counts)[ct_counts >= 10]
if (length(keep_ct) < ncol(obj)) {
  message(sprintf("Dropping %d cells (NA or tiny clusters) before slingshot",
                  ncol(obj) - sum(obj$celltype %in% keep_ct)))
  obj <- obj[, obj$celltype %in% keep_ct]
}
sce <- as.SingleCellExperiment(obj)
sce <- slingshot(sce, clusterLabels = "celltype", reducedDim = "UMAP",
                 start.clus = NULL)        # TODO: set biologically-known root cluster
saveRDS(sce, out_path(cfg, "objects", "spoke_slingshot.rds"))
log_session(cfg, "spoke_slingshot")
# next: tradeSeq::fitGAM for trajectory-dependent DE along lineages.
