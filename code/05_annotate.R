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
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(patchwork))
cfg <- load_config(); set.seed(cfg$seed)

obj <- readRDS(obj_path(cfg, "04_clustered"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)

# ---- automated annotation: three SingleR references (Part 4) ----
# HPCA = broad immune atlas (primary label driving config curation).
# Monaco + DICE = immune-specific refs for cross-validation only;
# their labels are stored as supplementary columns, not curated.
if (cfg$annotate$method == "singler") {
  suppressPackageStartupMessages({
    library(SingleR); library(celldex); library(SingleCellExperiment)
  })
  sce <- as.SingleCellExperiment(obj)

  ref  <- get(cfg$annotate$singler_ref, asNamespace("celldex"))()
  pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
  obj$celltype_singler <- pred$labels
  obj$celltype         <- pred$labels
  rm(ref, pred); invisible(gc())

  ref  <- MonacoImmuneData()
  pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine)
  obj$celltype_monaco <- pred$labels
  rm(ref, pred); invisible(gc())

  ref  <- DatabaseImmuneCellExpressionData()
  pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
  obj$celltype_dice <- pred$labels
  rm(ref, pred, sce); invisible(gc())
} else {
  manual_map <- c("0" = "CD4 T", "1" = "CD14 Mono", "2" = "B")   # TODO
  obj$celltype <- manual_map[as.character(Idents(obj))]
}

# =====================================================================
# LABEL CURATION — driven by config.yaml (annotate: rename: / drop:).
# Applies ONLY to obj$celltype (the primary HPCA-based label).
# Monaco / DICE / scType columns are left uncurated for comparison.
# =====================================================================
ren <- cfg$annotate$rename
if (!is.null(ren)) for (old in names(ren)) {
  n <- sum(obj$celltype == old)
  obj$celltype[obj$celltype == old] <- ren[[old]]
  message(sprintf("renamed %-18s -> %-18s (%d cells)", old, ren[[old]], n))
}
drop <- unlist(cfg$annotate$drop)
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
dp <- DotPlot(obj, features = pbmc_markers) +
      RotatedAxis() +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y  = element_text(size = 9),
            legend.title = element_text(size = 9),
            legend.text  = element_text(size = 8)) +
      labs(x = "Features", y = "Identity")
ggsave(out_path(cfg, "plots", "05_marker_dotplot.png"), dp, width = 12, height = 6)
ggsave(out_path(cfg, "plots", "05_marker_dotplot.pdf"), dp, width = 12, height = 6)
rm(dp); invisible(gc())

# ---- cluster markers: filter ribo/MT, save top 5 per cluster ----
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                          max.cells.per.ident = 500)
markers <- markers[!grepl("^RP[SL]|^MT-", markers$gene), ]
write.csv(markers, out_path(cfg, "tables", "05_cluster_markers.csv"), row.names = FALSE)

top5 <- do.call(rbind, lapply(split(markers, markers$cluster), function(x)
  head(x[order(-x$avg_log2FC), ], 5)))
write.csv(top5, out_path(cfg, "tables", "05_top5_markers_per_cluster.csv"), row.names = FALSE)
rm(markers); invisible(gc())

# ---- top-5 marker heatmap (200 cells/cluster subsample for readability) ----
set.seed(cfg$seed)
cells_hm <- unlist(lapply(levels(Idents(obj)), function(cl)
  sample(WhichCells(obj, idents = cl), min(200L, sum(Idents(obj) == cl)))))
p_hm <- DoHeatmap(obj[, cells_hm], features = unique(top5$gene)) + NoLegend()
ggsave(out_path(cfg, "plots", "05_top5_marker_heatmap.png"), p_hm, width = 14, height = 10)
ggsave(out_path(cfg, "plots", "05_top5_marker_heatmap.pdf"), p_hm, width = 14, height = 10)
rm(p_hm, cells_hm, top5); invisible(gc())

# ---- scType: marker-scoring annotation (Immune system gene sets) ----
# Helper scripts and DB are cached to xdisk on the first run, then reused.
tryCatch({
  message("Running scType...")
  sctype_dir <- file.path(cfg$project$outdir, "references", "sctype")
  dir.create(sctype_dir, showWarnings = FALSE, recursive = TRUE)
  base_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master"
  for (fn in c("R/gene_sets_prepare.R", "R/sctype_score_.R")) {
    lp <- file.path(sctype_dir, basename(fn))
    if (!file.exists(lp)) download.file(paste0(base_url, "/", fn), lp, quiet = TRUE)
    source(lp)
  }
  db_path <- file.path(sctype_dir, "ScTypeDB_full.xlsx")
  if (!file.exists(db_path))
    download.file(paste0(base_url, "/ScTypeDB_full.xlsx"), db_path, mode = "wb", quiet = TRUE)

  gs <- gene_sets_prepare(db_path, "Immune system")

  # Scale only the genes scType needs — avoids full 30K-gene ScaleData
  st_genes <- intersect(unique(c(unlist(gs$gs_positive), unlist(gs$gs_negative))),
                        rownames(obj))
  obj_st   <- ScaleData(obj, features = st_genes, verbose = FALSE)
  smat     <- tryCatch(GetAssayData(obj_st, layer = "scale.data"),
                       error = function(e) GetAssayData(obj_st, slot = "scale.data"))
  es <- sctype_score(scRNAseqData = smat, scaled = TRUE,
                     gs = gs$gs_positive, gs2 = gs$gs_negative)
  rm(obj_st, smat); invisible(gc())

  clust_res <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
    cells <- WhichCells(obj, idents = cl)
    sc    <- sort(rowSums(es[, cells, drop = FALSE]), decreasing = TRUE)
    data.frame(cluster = cl, type = names(sc)[1], score = sc[1],
               ncells = length(cells), stringsAsFactors = FALSE)
  }))
  clust_res$type_final <- ifelse(clust_res$score < clust_res$ncells / 4,
                                 "Unknown", clust_res$type)
  obj$celltype_sctype <- clust_res$type_final[
    match(as.character(Idents(obj)), clust_res$cluster)]
  write.csv(clust_res, out_path(cfg, "tables", "05_sctype_cluster_scores.csv"),
            row.names = FALSE)
  cat("\n=== scType cluster assignments ===\n")
  print(clust_res[, c("cluster", "type_final", "score", "ncells")])
  rm(es, gs, clust_res); invisible(gc())
}, error = function(e) {
  message("scType skipped: ", conditionMessage(e))
})

# ---- save the canonical hub object + lightweight CSV export for velocity ----
saveRDS(obj, obj_path(cfg, "05_annotated"))
export_meta_for_velocity(obj, cfg)

# ---- UMAPs ----
# Primary celltype UMAP
p_ct <- DimPlot(obj, label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.3) +
        theme(plot.title = element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "05_umap_celltypes.png"), p_ct, width = 10, height = 7)
ggsave(out_path(cfg, "plots", "05_umap_celltypes.pdf"), p_ct, width = 10, height = 7)
rm(p_ct); invisible(gc())

# Condition-split UMAP: check biological plausibility per condition
n_cond  <- length(unique(obj$condition))
p_split <- DimPlot(obj, group.by = "celltype", split.by = "condition",
                   label = TRUE, repel = TRUE, label.size = 2.5, pt.size = 0.2) +
           NoLegend() +
           theme(strip.text = element_text(size = 11))
ggsave(out_path(cfg, "plots", "05_umap_celltypes_by_condition.png"),
       p_split, width = 5 * n_cond, height = 5)
ggsave(out_path(cfg, "plots", "05_umap_celltypes_by_condition.pdf"),
       p_split, width = 5 * n_cond, height = 5)
rm(p_split); invisible(gc())

# Multi-reference comparison: HPCA / Monaco / DICE / scType side by side
ref_map <- c(HPCA   = "celltype_singler",
             Monaco = "celltype_monaco",
             DICE   = "celltype_dice",
             scType = "celltype_sctype")
ref_map <- ref_map[ref_map %in% colnames(obj@meta.data)]
if (length(ref_map) > 0) {
  panels <- lapply(names(ref_map), function(nm)
    DimPlot(obj, group.by = ref_map[nm], label = TRUE, repel = TRUE,
            label.size = 2.5, pt.size = 0.2) +
    ggtitle(nm) + NoLegend() +
    theme(plot.title = element_text(hjust = 0.5)))
  ncol_p  <- min(2L, length(panels))
  nrow_p  <- ceiling(length(panels) / ncol_p)
  p_multi <- wrap_plots(panels, ncol = ncol_p)
  ggsave(out_path(cfg, "plots", "05_umap_all_methods.png"),
         p_multi, width = 10 * ncol_p, height = 7 * nrow_p)
  ggsave(out_path(cfg, "plots", "05_umap_all_methods.pdf"),
         p_multi, width = 10 * ncol_p, height = 7 * nrow_p)
  rm(panels, p_multi); invisible(gc())

  # All-methods comparison table (cell-level)
  write.csv(obj@meta.data[, c("sample", "condition", "celltype", unname(ref_map))],
            out_path(cfg, "tables", "05_all_methods_comparison.csv"))
}

log_session(cfg, "05_annotate")
