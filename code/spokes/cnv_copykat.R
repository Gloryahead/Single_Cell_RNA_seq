#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# cnv_copykat.R — CNV inference to flag aneuploid cells (CopyKAT). (Part 11)
#
# NOTE: CopyKAT is NOT in the current scrna_r.sif. Rebuild with:
#   devtools::install_github("navinlabcode/copykat")
# before enabling this rule in Snakefile rule all.
#
# PBMC note: CopyKAT is designed for tumor data. On healthy PBMCs all
# cells should be diploid; elevated aneuploid fractions in a cell type
# indicate either a quality issue or an unexpected population. This
# run is a methods demo and a QC cross-check.
#
# Adaptation: use known-diploid cell types (T_cells, B_cells, NK_cells)
# as norm.cell.names — more reliable than the tutorial's tumor/normal
# sample split for a non-cancer dataset.
#
# In : results/objects/05_annotated.rds
# Out: results/objects/spoke_copykat.rds  (Seurat obj + cnv_status column)
#      results/copykat/<sample>/          (CopyKAT native files)
#      results/tables/spoke_copykat_summary.csv
#      results/plots/spoke_copykat_*.{png,pdf}
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(copykat); library(patchwork)
})
cfg      <- load_config(); set.seed(cfg$seed)
n_cores  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
obj      <- readRDS(obj_path(cfg, "05_annotated"))

# ---- Define diploid reference cell types --------------------------------
# Known diploid populations in PBMCs: T cells, B cells, NK cells.
# CopyKAT uses these as the diploid baseline for CNV calling.
diploid_pattern <- "T_cell|NK_cell|B_cell|T cell|NK cell|B cell"
diploid_types   <- grep(diploid_pattern, unique(obj$celltype),
                        value = TRUE, ignore.case = TRUE)
cat("Diploid reference cell types: ", paste(diploid_types, collapse = ", "), "\n")

# ---- Run CopyKAT per sample ---------------------------------------------
# Running on all cells at once is prohibitively slow (hours per 100K cells).
# Per-sample execution (8-10K cells each) is ~10-20 min/sample at 12 cores.
samples     <- unique(obj$sample)
all_pred    <- list()    # barcode -> prediction per sample
all_cna     <- list()    # CNAmat per sample (for subclone analysis)

ck_base_dir <- file.path(cfg$project$outdir, "copykat")

for (samp in samples) {
  message("\n=== CopyKAT: ", samp, " ===")
  pred_file <- file.path(ck_base_dir, samp,
                         paste0(samp, "_copykat_prediction.txt"))

  if (file.exists(pred_file)) {
    message("  Prediction file exists — loading cached result.")
    pred <- read.table(pred_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
    all_pred[[samp]] <- pred
    cna_file <- file.path(ck_base_dir, samp,
                          paste0(samp, "_copykat_CNA_results.txt"))
    if (file.exists(cna_file))
      all_cna[[samp]] <- read.table(cna_file, header = TRUE,
                                    sep = "\t", check.names = FALSE)
    next
  }

  samp_cells  <- colnames(obj)[obj$sample == samp]
  samp_obj    <- obj[, samp_cells]
  raw_counts  <- as.matrix(
    tryCatch(GetAssayData(samp_obj, layer = "counts"),
             error = function(e) GetAssayData(samp_obj, slot  = "counts")))

  # Diploid reference barcodes from this sample
  norm_bcs <- samp_cells[samp_obj$celltype %in% diploid_types]
  cat(sprintf("  %s: %d cells total, %d diploid references\n",
              samp, ncol(samp_obj), length(norm_bcs)))

  out_dir <- file.path(ck_base_dir, samp)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sam_pref <- file.path(out_dir, samp)

  ck_res <- tryCatch(
    copykat(rawmat          = raw_counts,
            id.type         = "S",
            genome          = "hg20",
            ngene.chr       = 5,
            win.size        = 25,
            KS.cut          = 0.1,
            LOW.DR          = 0.05,
            UP.DR           = 0.1,
            distance        = "euclidean",
            norm.cell.names = norm_bcs,
            cell.line       = "no",
            n.cores         = n_cores,
            sam.name        = sam_pref),
    error = function(e) {
      message("  ERROR running CopyKAT on ", samp, ": ", conditionMessage(e))
      NULL
    }
  )

  rm(samp_obj, raw_counts); invisible(gc())
  if (is.null(ck_res)) next

  all_pred[[samp]] <- ck_res$prediction
  if (!is.null(ck_res$CNAmat)) all_cna[[samp]] <- ck_res$CNAmat
  rm(ck_res); invisible(gc())
}

if (length(all_pred) == 0)
  stop("No CopyKAT predictions were produced for any sample.")

# ---- Merge predictions back to Seurat object ----------------------------
# Combine all per-sample prediction tables
pred_all <- do.call(rbind, all_pred)
cat(sprintf("\nCopyKAT predictions: %d cells\n", nrow(pred_all)))
cat("Label distribution:\n"); print(table(pred_all$copykat.pred))

# Map to obj by barcode
obj$cnv_status <- "not.defined"
m <- match(pred_all$cell.names, colnames(obj))
valid <- !is.na(m)
obj$cnv_status[m[valid]] <- pred_all$copykat.pred[valid]
# Fill cells without a prediction (rare, can happen if sample failed)
obj$cnv_status[is.na(obj$cnv_status)] <- "not.defined"
cat("\nAneuploidy by condition:\n")
print(table(obj$condition, obj$cnv_status))

# ---- Aneuploid fraction per cell type -----------------------------------
aneu_frac <- tapply(obj$cnv_status, obj$celltype, function(x)
  round(100 * mean(x == "aneuploid"), 2))
aneu_frac <- sort(aneu_frac, decreasing = TRUE)
cat("\nAneuploid fraction per cell type (%):\n"); print(aneu_frac)

# Per-sample aneuploid fraction summary
aneu_summary <- do.call(rbind, lapply(samples, function(s) {
  idx <- obj$sample == s
  tab <- table(obj$celltype[idx], obj$cnv_status[idx])
  tot <- rowSums(tab)
  aneu_n <- if ("aneuploid" %in% colnames(tab)) tab[, "aneuploid"] else
            rep(0L, nrow(tab))
  data.frame(sample    = s,
             condition = unique(obj$condition[idx]),
             celltype  = rownames(tab),
             n_cells   = tot,
             n_aneuploid = aneu_n,
             pct_aneuploid = round(100 * aneu_n / pmax(1L, tot), 2),
             stringsAsFactors = FALSE)
}))
write.csv(aneu_summary,
          out_path(cfg, "tables", "spoke_copykat_summary.csv"),
          row.names = FALSE)

# ---- Plots --------------------------------------------------------------
ct_pal  <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                     "#A65628","#F781BF","#999999"))(length(unique(obj$celltype))),
  sort(unique(obj$celltype)))
cnv_pal <- c(aneuploid = "#E41A1C", diploid = "#377EB8",
             not.defined = "grey80")

# UMAP: cell type vs CopyKAT label, split by condition
p_ct  <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE,
                 label.size = 3, pt.size = 0.3) +
          scale_color_manual(values = ct_pal) +
          ggtitle("Cell types") + NoLegend()
p_cnv <- DimPlot(obj, group.by = "cnv_status", pt.size = 0.3) +
          scale_color_manual(values = cnv_pal) +
          ggtitle("CopyKAT status")
ggsave(out_path(cfg, "plots", "spoke_copykat_umap.png"),
       p_ct + p_cnv, width = 14, height = 6)
ggsave(out_path(cfg, "plots", "spoke_copykat_umap.pdf"),
       p_ct + p_cnv, width = 14, height = 6)
rm(p_ct, p_cnv); invisible(gc())

# Aneuploid fraction barplot per cell type
aneu_df <- data.frame(celltype = names(aneu_frac),
                      pct      = as.numeric(aneu_frac))
aneu_df$celltype <- factor(aneu_df$celltype, levels = rev(names(aneu_frac)))
p_bar <- ggplot2::ggplot(aneu_df, ggplot2::aes(x = celltype, y = pct,
                                                fill = celltype)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::scale_fill_manual(values = ct_pal) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Aneuploid cells (%)",
                title = "CopyKAT aneuploid fraction by cell type") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "spoke_copykat_aneuploidy_barplot.png"),
       p_bar, width = 8, height = 5)
ggsave(out_path(cfg, "plots", "spoke_copykat_aneuploidy_barplot.pdf"),
       p_bar, width = 8, height = 5)
rm(p_bar); invisible(gc())

# Per-sample aneuploid fraction heatmap
tryCatch({
  hm_mat <- reshape(aneu_summary[, c("sample","celltype","pct_aneuploid")],
                    idvar = "celltype", timevar = "sample", direction = "wide")
  rownames(hm_mat) <- hm_mat$celltype; hm_mat$celltype <- NULL
  colnames(hm_mat) <- gsub("pct_aneuploid\\.", "", colnames(hm_mat))
  hm_mat <- as.matrix(hm_mat)

  col_fun <- circlize::colorRamp2(
    c(0, 25, 50, 100),
    c("white", "#FDAE6B", "#E41A1C", "#67000D"))
  png(out_path(cfg, "plots", "spoke_copykat_heatmap.png"),
      width = 1100, height = 700, res = 120)
  print(ComplexHeatmap::Heatmap(
    hm_mat, name = "% aneuploid",
    col  = col_fun,
    row_names_gp  = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    column_title  = "CopyKAT aneuploid fraction (%) per sample × cell type",
    cluster_rows  = FALSE, cluster_columns = FALSE))
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message("[skip] per-sample heatmap: ", conditionMessage(e))
})

# ---- Subclone detection (only if ≥100 aneuploid cells) -----------------
n_aneu <- sum(obj$cnv_status == "aneuploid")
if (n_aneu >= 100 && length(all_cna) > 0) {
  tryCatch({
    message(sprintf("\nSubclone detection: %d aneuploid cells", n_aneu))
    # Pool CNAmat across samples; columns 1:6 are annotation, 7:ncol are cells
    cna_pool <- do.call(cbind, lapply(all_cna, function(m) {
      cell_cols <- seq(7, ncol(m))
      m[, cell_cols, drop = FALSE]
    }))
    aneu_cells <- colnames(obj)[obj$cnv_status == "aneuploid"]
    aneu_cells <- intersect(aneu_cells, colnames(cna_pool))
    if (length(aneu_cells) >= 50) {
      cna_aneu <- t(cna_pool[, aneu_cells, drop = FALSE])
      hc  <- hclust(dist(cna_aneu), method = "ward.D2")
      sub <- cutree(hc, k = 2)
      obj$subclone <- "non-aneuploid"
      obj$subclone[match(names(sub), colnames(obj))] <- paste0("subclone_", sub)
      cat("\nSubclone distribution:\n"); print(table(obj$subclone))

      p_sub <- DimPlot(obj, group.by = "subclone", pt.size = 0.3) +
               ggtitle("CopyKAT subclones")
      ggsave(out_path(cfg, "plots", "spoke_copykat_subclones.png"),
             p_sub, width = 8, height = 6)
      ggsave(out_path(cfg, "plots", "spoke_copykat_subclones.pdf"),
             p_sub, width = 8, height = 6)
      rm(p_sub, cna_aneu, hc, sub); invisible(gc())

      # DE between subclones
      sub_obj <- subset(obj, subclone != "non-aneuploid")
      Idents(sub_obj) <- "subclone"
      sub_markers <- tryCatch(
        FindAllMarkers(sub_obj, assay = "RNA", only.pos = TRUE,
                       min.pct = 0.25, logfc.threshold = 0.25,
                       verbose = FALSE),
        error = function(e) { message("[skip] subclone DE: ", conditionMessage(e)); NULL })
      if (!is.null(sub_markers) && nrow(sub_markers) > 0) {
        write.csv(sub_markers,
                  out_path(cfg, "tables", "spoke_copykat_subclone_markers.csv"),
                  row.names = FALSE)
        cat("\nTop subclone markers:\n")
        print(do.call(rbind, lapply(split(sub_markers, sub_markers$cluster),
                                   function(x) head(x[order(-x$avg_log2FC), ], 5))))
      }
      rm(sub_obj); invisible(gc())
    }
    rm(cna_pool); invisible(gc())
  }, error = function(e)
    message("[skip] subclone detection: ", conditionMessage(e)))
} else {
  message(sprintf("Subclone detection skipped: %d aneuploid cells (need ≥100)", n_aneu))
}

# ---- Save ---------------------------------------------------------------
saveRDS(obj, out_path(cfg, "objects", "spoke_copykat.rds"))
log_session(cfg, "spoke_copykat")
