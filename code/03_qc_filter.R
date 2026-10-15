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
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(patchwork))
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
# capture per-sample raw counts before any QC filtering
cells_raw <- vapply(objs, ncol, integer(1L))

# ---- QC plots BEFORE filtering: thresholds are decisions made from these ----
pre <- merge(objs[[1]], objs[-1])
if (inherits(pre[["RNA"]], "Assay5")) pre <- JoinLayers(pre)
# build panels individually so each gets its own threshold lines (red dashed)
p_nfeat  <- VlnPlot(pre, "nFeature_RNA",  group.by = "sample", pt.size = 0) + NoLegend() +
             geom_hline(yintercept = c(cfg$qc$min_genes, cfg$qc$max_genes),
                        linetype = "dashed", colour = "red", linewidth = 0.5)
p_ncount <- VlnPlot(pre, "nCount_RNA",    group.by = "sample", pt.size = 0) + NoLegend()
p_mt     <- VlnPlot(pre, "percent.mt",    group.by = "sample", pt.size = 0) + NoLegend() +
             geom_hline(yintercept = cfg$qc$max_mito_pct,
                        linetype = "dashed", colour = "red", linewidth = 0.5)
p_ribo   <- VlnPlot(pre, "percent.ribo",  group.by = "sample", pt.size = 0) + NoLegend()
p_pre    <- wrap_plots(p_nfeat, p_ncount, p_mt, p_ribo, nrow = 1)
ggsave(out_path(cfg, "qc", "qc_violin_prefilter.png"), p_pre, width = 16, height = 4)
ggsave(out_path(cfg, "qc", "qc_violin_prefilter.pdf"), p_pre, width = 16, height = 4)
# scatter reveals the doublet corner (high counts AND high genes)
p_scat <- FeatureScatter(pre, "nCount_RNA", "nFeature_RNA", group.by = "sample")
ggsave(out_path(cfg, "qc", "qc_scatter_prefilter.png"), p_scat, width = 7, height = 5)
ggsave(out_path(cfg, "qc", "qc_scatter_prefilter.pdf"), p_scat, width = 7, height = 5)
rm(pre, p_nfeat, p_ncount, p_mt, p_ribo, p_pre, p_scat); invisible(gc())

# ---- filter + flag doublets, per sample (Parts 2 / 2-2) ----
dbl_stats    <- list()
filter_stats <- list()
objs <- lapply(names(objs), function(id) {
  obj    <- apply_qc_filter(objs[[id]], cfg)       # nFeature + percent.mt cutoffs
  n_post <- ncol(obj)
  filter_stats[[id]] <<- data.frame(
    sample         = id,
    cells_raw      = cells_raw[[id]],
    cells_after_qc = n_post,
    pct_kept       = round(100 * n_post / cells_raw[[id]], 1)
  )
  if (isTRUE(cfg$qc$remove_doublets)) {
    obj   <- run_scdblfinder(obj)                  # adds $doublet: singlet/doublet
    n_dbl <- sum(obj$doublet == "doublet")
    dbl_stats[[id]] <<- data.frame(sample = id, cells = n_post,
                                   doublets = n_dbl,
                                   pct = round(100 * n_dbl / n_post, 2))
  }
  obj
})
names(objs) <- vapply(cfg$samples, `[[`, "", "id")
filter_summary <- do.call(rbind, filter_stats)
write.csv(filter_summary, out_path(cfg, "qc", "filter_summary.csv"), row.names = FALSE)
cat("\n=== cells per sample after QC filtering ===\n"); print(filter_summary)

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
  p_dumap <- DimPlot(tmp, group.by = "doublet", order = "doublet",
                     cols = c(singlet = "grey85", doublet = "red")) +
             ggtitle("scDblFinder calls (red = doublet)")
  ggsave(out_path(cfg, "qc", "doublet_umap.png"), p_dumap, width = 7, height = 6)
  ggsave(out_path(cfg, "qc", "doublet_umap.pdf"), p_dumap, width = 7, height = 6)
  # doublets should have HIGHER counts/genes than singlets — sanity check
  p_dvln <- VlnPlot(tmp, c("nCount_RNA","nFeature_RNA"), group.by = "doublet", pt.size = 0)
  ggsave(out_path(cfg, "qc", "doublet_violin.png"), p_dvln, width = 9, height = 4)
  ggsave(out_path(cfg, "qc", "doublet_violin.pdf"), p_dvln, width = 9, height = 4)
  rm(tmp, p_dumap, p_dvln); invisible(gc())

  merged <- subset(merged, subset = doublet == "singlet")   # NOW remove them
}

# ---- QC plots AFTER filtering: confirm distributions are clean ----
p_post <- wrap_plots(
  VlnPlot(merged, "nFeature_RNA",  group.by = "sample", pt.size = 0) + NoLegend(),
  VlnPlot(merged, "nCount_RNA",    group.by = "sample", pt.size = 0) + NoLegend(),
  VlnPlot(merged, "percent.mt",    group.by = "sample", pt.size = 0) + NoLegend(),
  VlnPlot(merged, "percent.ribo",  group.by = "sample", pt.size = 0) + NoLegend(),
  nrow = 1
)
ggsave(out_path(cfg, "qc", "qc_violin_postfilter.png"), p_post, width = 16, height = 4)
ggsave(out_path(cfg, "qc", "qc_violin_postfilter.pdf"), p_post, width = 16, height = 4)
rm(p_post); invisible(gc())

cat(sprintf("\nCells after QC: %d | genes: %d\n", ncol(merged), nrow(merged)))
saveRDS(merged, obj_path(cfg, "03_filtered"))
log_session(cfg, "03_qc_filter")
