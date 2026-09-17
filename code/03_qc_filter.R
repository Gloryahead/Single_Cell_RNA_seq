#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 03_qc_filter.R — Cell Ranger output -> clean, filtered object.
#                                                  (Tutorial Parts 2, 2-2)
# Inputs  : results/counts/<id>/outs/filtered_feature_bc_matrix/
# Outputs : results/objects/03_filtered.rds
#           results/qc/  (violins, doublet UMAP, doublet_summary.csv)
# Run     : Rscript code/03_qc_filter.R
#
# NOTE: Cell Ranger's filtered_feature_bc_matrix has ALREADY removed empty
# droplets, so "pre-filter" below means "before OUR thresholds", not raw.
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
cfg <- load_config(); set.seed(cfg$seed)

# ---- load each sample, tag metadata, compute QC metrics ----
objs <- lapply(cfg$samples, function(s) {
  mtx <- file.path(cfg$project$outdir, "counts", s$id,
                   "outs", "filtered_feature_bc_matrix")
  obj <- CreateSeuratObject(Read10X(mtx), project = s$id, min.cells = cfg$qc$min_cells)
  obj$sample <- s$id; obj$condition <- s$condition
  add_qc_metrics(obj, cfg)          # percent.mt + percent.ribo
})
names(objs) <- vapply(cfg$samples, `[[`, "", "id")

# ---- QC plots BEFORE filtering: thresholds are decisions made from these ----
pre <- merge(objs[[1]], objs[-1])
if (inherits(pre[["RNA"]], "Assay5")) pre <- JoinLayers(pre)
ggsave(out_path(cfg, "qc", "qc_violin_prefilter.png"),
       VlnPlot(pre, c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
               group.by = "sample", pt.size = 0, ncol = 4), width = 16, height = 4)
# scatter reveals the doublet corner (high counts AND high genes)
ggsave(out_path(cfg, "qc", "qc_scatter_prefilter.png"),
       FeatureScatter(pre, "nCount_RNA", "nFeature_RNA", group.by = "sample"),
       width = 7, height = 5)
rm(pre); invisible(gc())

# ---- filter + flag doublets, per sample (Parts 2 / 2-2) ----
dbl_stats <- list()
objs <- lapply(names(objs), function(id) {
  obj <- apply_qc_filter(objs[[id]], cfg)          # nFeature + percent.mt cutoffs
  n_before <- ncol(obj)
  if (isTRUE(cfg$qc$remove_doublets)) {
    obj <- run_scdblfinder(obj)                    # adds $doublet: singlet/doublet
    n_dbl <- sum(obj$doublet == "doublet")
    dbl_stats[[id]] <<- data.frame(sample = id, cells = n_before,
                                   doublets = n_dbl,
                                   pct = round(100 * n_dbl / n_before, 2))
  }
  obj
})
names(objs) <- vapply(cfg$samples, `[[`, "", "id")

merged <- if (length(objs) == 1) objs[[1]] else merge(objs[[1]], objs[-1], add.cell.ids = names(objs))
if (inherits(merged[["RNA"]], "Assay5")) merged <- JoinLayers(merged)

# ---- VALIDATE the doublet calls before trusting them ----------------
# Real doublets sit at cluster boundaries / form intermediate blobs. If they're
# scattered randomly through cluster interiors, detection failed.
if (isTRUE(cfg$qc$remove_doublets)) {
  stats <- do.call(rbind, dbl_stats)
  write.csv(stats, out_path(cfg, "qc", "doublet_summary.csv"), row.names = FALSE)
  cat("\n=== doublets flagged per sample ===\n"); print(stats)
  cat(sprintf("\nTOTAL: %d doublets / %d cells (%.2f%%)\n",
              sum(stats$doublets), sum(stats$cells),
              100 * sum(stats$doublets) / sum(stats$cells)))
  cat("Expect roughly ~0.8%% per 1000 cells loaded (10x). Far less => under-calling.\n")

  # quick embedding purely to LOOK at where the doublets landed
  tmp <- NormalizeData(merged, verbose = FALSE) |>
         FindVariableFeatures(verbose = FALSE) |> ScaleData(verbose = FALSE) |>
         RunPCA(npcs = 20, verbose = FALSE) |> RunUMAP(dims = 1:20, verbose = FALSE)
  ggsave(out_path(cfg, "qc", "doublet_umap.png"),
         DimPlot(tmp, group.by = "doublet", order = "doublet",
                 cols = c(singlet = "grey85", doublet = "red")) +
           ggtitle("scDblFinder calls (red = doublet)"), width = 7, height = 6)
  # doublets should have HIGHER counts/genes than singlets — sanity check
  ggsave(out_path(cfg, "qc", "doublet_violin.png"),
         VlnPlot(tmp, c("nCount_RNA","nFeature_RNA"), group.by = "doublet", pt.size = 0),
         width = 9, height = 4)
  rm(tmp); invisible(gc())

  merged <- subset(merged, subset = doublet == "singlet")   # NOW remove them
}

cat(sprintf("\nCells after QC: %d | genes: %d\n", ncol(merged), nrow(merged)))
saveRDS(merged, obj_path(cfg, "03_filtered"))
log_session(cfg, "03_qc_filter")
