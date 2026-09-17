#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 05_annotate.R  ★ THE HUB — clustered -> annotated object.
#                                                       (Tutorial Part 4)
# SingleR labels are a DRAFT. The marker dotplot + diagnose_labels.R are
# what you trust. Curation happens in the LABEL CURATION block below.
# Inputs  : results/objects/04_clustered.rds
# Outputs : results/objects/05_annotated.rds  (+ .h5ad for Python spokes)
# Run     : Rscript code/05_annotate.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
cfg <- load_config(); set.seed(cfg$seed)

obj <- readRDS(obj_path(cfg, "04_clustered"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)   # v5: collapse layers

# ---- automated annotation: SingleR + an immune-tuned reference (Part 4) ----
if (cfg$annotate$method == "singler") {
  suppressPackageStartupMessages({ library(SingleR); library(celldex); library(SingleCellExperiment) })
  ref  <- get(cfg$annotate$singler_ref, asNamespace("celldex"))()
  sce  <- as.SingleCellExperiment(obj)
  pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
  obj$celltype_singler <- pred$labels          # keep the raw draft for the record
  obj$celltype         <- pred$labels
} else {
  manual_map <- c("0" = "CD4 T", "1" = "CD14 Mono", "2" = "B")   # TODO
  obj$celltype <- manual_map[as.character(Idents(obj))]
}

# =====================================================================
# LABEL CURATION — driven by config.yaml (annotate: rename: / drop:).
# Run diagnose_labels.R FIRST, then encode your decisions in config.
# Evidence for the suggested defaults:
#   * "Progenitors" expressing PPBP/PF4 are PLATELETS (unambiguous).
#   * Neutrophils don't survive PBMC density-gradient prep -> artifact.
#   * Basophils here are a handful of cells -> likely artifact.
# The CD8 <-> "T cells" swap is NOT auto-applied — confirm it with
# diagnose_labels.R, then add it to `rename:` yourself.
# =====================================================================
ren <- cfg$annotate$rename                     # named list: old -> new
if (!is.null(ren)) for (old in names(ren)) {
  n <- sum(obj$celltype == old)
  obj$celltype[obj$celltype == old] <- ren[[old]]
  message(sprintf("renamed %-18s -> %-18s (%d cells)", old, ren[[old]], n))
}
drop <- unlist(cfg$annotate$drop)              # vector of labels to remove
if (!is.null(drop)) {
  keep <- !(obj$celltype %in% drop)
  message(sprintf("dropping %d cells labelled: %s",
                  sum(!keep), paste(drop, collapse = ", ")))
  obj <- obj[, keep]
}
Idents(obj) <- "celltype"
cat("\nfinal cell types:\n"); print(sort(table(obj$celltype)))

# ---- canonical PBMC markers: ALWAYS eyeball these to validate labels ----
pbmc_markers <- c("CD3D","IL7R","CD4","CD8A","CD8B","GZMK","GNLY","NKG7",
                  "MS4A1","CD79A","CD14","LYZ","FCGR3A","MS4A7",
                  "FCER1A","CST3","PPBP","PF4")
pbmc_markers <- pbmc_markers[pbmc_markers %in% rownames(obj)]
ggsave(out_path(cfg, "plots", "05_marker_dotplot.png"),
       DotPlot(obj, features = pbmc_markers) + RotatedAxis(), width = 12, height = 6)

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, out_path(cfg, "tables", "05_cluster_markers.csv"), row.names = FALSE)

# ---- save the canonical hub object (+ AnnData copy for velocity/fate) ----
saveRDS(obj, obj_path(cfg, "05_annotated"))
export_h5ad(obj, out_path(cfg, "objects", "05_annotated.h5ad"))
ggsave(out_path(cfg, "plots", "05_umap_celltypes.png"),
       DimPlot(obj, label = TRUE), width = 7, height = 6)
log_session(cfg, "05_annotate")
