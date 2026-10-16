#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# viz_scplotter.R — publication-ready figures (scplotter).        (Part 15)
#
# NOTE: scplotter is NOT in the current scrna_r.sif. Rebuild with:
#   remotes::install_github("pwwang/scplotter")
#   install.packages(c("clustree","gglogger","metR","proxyC"))
#
# This spoke re-draws the key figures from across the pipeline using
# scplotter's one-function-per-category grammar — no patchwork gymnastics,
# consistent defaults, publication-ready output. It does NOT change any
# analysis result — only how the results are visualised.
#
# In : results/objects/05_annotated.rds
# Out: results/plots/spoke_scplotter_*.{png,pdf}
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(scplotter); library(Seurat); library(ggplot2); library(dplyr)
})
cfg <- load_config(); set.seed(cfg$seed)
obj  <- readRDS(obj_path(cfg, "05_annotated"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)

# Free the dense scaled layer (it is not needed for any plot here; keeping it
# risks OOM on the dot/heatmap steps that build dense cells-x-genes matrices)
tryCatch({
  if ("scale.data" %in% Layers(obj))
    obj[["RNA"]]$scale.data <- NULL
}, error = function(e) invisible(NULL))
invisible(gc())

# Detect UMAP reduction (Harmony → "umap"; fall back to first UMAP found)
avail_red   <- Reductions(obj)
umap_red    <- dplyr::case_when(
  "umap"         %in% avail_red ~ "umap",
  "umap.harmony" %in% avail_red ~ "umap.harmony",
  "umap.cca"     %in% avail_red ~ "umap.cca",
  TRUE                          ~ avail_red[grep("umap", tolower(avail_red))[1]]
)
message("Using UMAP reduction: ", umap_red)

# Detect cluster resolution columns (ClustreePlot)
res_cols <- grep("^RNA_snn_res\\.|^wsnn_res\\.|seurat_clusters",
                 colnames(obj@meta.data), value = TRUE)
res_prefix <- if (any(grepl("RNA_snn_res\\.", res_cols))) "RNA_snn_res."
              else if (any(grepl("wsnn_res\\.",    res_cols))) "wsnn_res."
              else NULL

# Primary contrast from config (Post vs Healthy)
primary_grp <- cfg$de_contrasts[[1]]$group  # e.g. "Post"
primary_ref <- cfg$de_contrasts[[1]]$ref    # e.g. "Healthy"

# Canonical PBMC marker genes (adjust if HGNC symbols differ in your genome)
pbmc_markers <- c(
  "CD3D",   # T cells (pan)
  "CD4",    # CD4+ T cells
  "CD8A",   # CD8+ T cells
  "NKG7",   # NK cells
  "GNLY",   # NK cells
  "CD79A",  # B cells
  "MS4A1",  # B cells (CD20)
  "CD14",   # Classical Monocytes
  "FCGR3A", # Non-classical Monocytes
  "S100A8"  # Myeloid
)
# Keep only markers actually present in the object
pbmc_markers <- intersect(pbmc_markers, rownames(obj))
if (length(pbmc_markers) < 3)
  message("WARNING: fewer than 3 canonical PBMC markers found. Check genome symbol format.")

# gg_sp: save ggplot objects
gg_sp <- function(p, name, w = 8, h = 7) {
  tryCatch(
    ggsave(out_path(cfg, "plots", paste0("spoke_scplotter_", name)), p,
           width = w, height = h, dpi = 150),
    error = function(e) message("[skip] ", name, ": ", conditionMessage(e))
  )
}

# png_sp: save ComplexHeatmap / grid objects (returned by dot/heatmap plot_type)
png_sp <- function(p, name, w = 900, h = 600) {
  path <- out_path(cfg, "plots", paste0("spoke_scplotter_", name))
  tryCatch({
    png(path, width = w, height = h, res = 120)
    print(p)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("[skip] ", name, ": ", conditionMessage(e))
  })
}

# ---- Step 5: QC metric violins (Part 2 recreation) -----------------------
tryCatch({
  p <- FeatureStatPlot(obj,
         features     = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
         ident        = "sample",
         plot_type    = "violin",
         add_box      = TRUE,
         x_text_angle = 45,
         facet_scales = "free_y")
  gg_sp(p, "01_qc_violins.png", w = 14, h = 5)
  message("Plot: QC violins")
}, error = function(e) message("[skip] QC violins: ", conditionMessage(e)))

# ---- Step 6: UMAP by cell type and by cluster ----------------------------
tryCatch({
  p_ct <- CellDimPlot(obj, group_by = "celltype", reduction = umap_red,
                      label = TRUE, label_repel = TRUE, label_insitu = TRUE)
  gg_sp(p_ct, "02_umap_celltypes.png")
  message("Plot: UMAP cell types")
}, error = function(e) message("[skip] celltype UMAP: ", conditionMessage(e)))

# Cluster UMAP — use seurat_clusters or first resolution column
clust_col <- if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters"
             else if (length(res_cols) > 0) res_cols[1]
             else NULL
if (!is.null(clust_col)) {
  tryCatch({
    p_cl <- CellDimPlot(obj, group_by = clust_col, reduction = umap_red, label = TRUE)
    gg_sp(p_cl, "03_umap_clusters.png")
    message("Plot: UMAP clusters (", clust_col, ")")
  }, error = function(e) message("[skip] cluster UMAP: ", conditionMessage(e)))
}

# ---- Step 7: UMAP by sample and by condition (integration check) ---------
tryCatch({
  p_s <- CellDimPlot(obj, group_by = "sample",    reduction = umap_red)
  p_c <- CellDimPlot(obj, group_by = "condition", reduction = umap_red, palette = "Set1")
  gg_sp(p_s, "04_umap_by_sample.png")
  gg_sp(p_c, "05_umap_by_condition.png")
  message("Plot: UMAP by sample + condition")
}, error = function(e) message("[skip] sample/condition UMAP: ", conditionMessage(e)))

# ---- Step 8: Split UMAP by condition (facet_by keeps consistent colors) --
tryCatch({
  p <- CellDimPlot(obj, group_by = "celltype", reduction = umap_red,
                   facet_by = "condition")
  gg_sp(p, "06_umap_split_condition.png", w = 18, h = 7)
  message("Plot: split UMAP by condition")
}, error = function(e) message("[skip] split UMAP: ", conditionMessage(e)))

# ---- Step 9: Highlight one cell type (auto-pick largest non-T population) -
tryCatch({
  ct_counts   <- sort(table(obj$celltype), decreasing = TRUE)
  t_like      <- grep("^T_|^T |CD4|CD8|Naive T|Memory T",
                       names(ct_counts), ignore.case = TRUE)
  highlight_ct <- if (length(t_like) < length(ct_counts))
                    names(ct_counts[setdiff(seq_along(ct_counts), t_like)[1]])
                  else names(ct_counts)[1]
  p <- CellDimPlot(obj, group_by = "celltype", reduction = umap_red,
                   highlight = paste0('celltype == "', highlight_ct, '"'))
  gg_sp(p, "07_umap_highlight.png")
  message("Plot: UMAP highlight (", highlight_ct, ")")
}, error = function(e) message("[skip] highlight UMAP: ", conditionMessage(e)))

# ---- Steps 11-14: Marker gene expression ---------------------------------
if (length(pbmc_markers) >= 3) {
  # Feature plots on UMAP
  tryCatch({
    p <- FeatureStatPlot(obj, features = pbmc_markers, plot_type = "dim",
                         reduction = umap_red, ncol = 5)
    gg_sp(p, "08_markers_featureplot.png", w = 20, h = ceiling(length(pbmc_markers)/5)*4)
    message("Plot: marker feature plots")
  }, error = function(e) message("[skip] feature plots: ", conditionMessage(e)))

  # Stacked violin
  tryCatch({
    p <- FeatureStatPlot(obj, features = pbmc_markers,
                         ident = "celltype", plot_type = "violin",
                         stack = TRUE, add_bg = TRUE, flip = TRUE)
    gg_sp(p, "09_markers_violin_stacked.png", w = 12, h = 8)
    message("Plot: stacked violin")
  }, error = function(e) message("[skip] stacked violin: ", conditionMessage(e)))

  # For dot/heatmap: subset to marker genes only (avoids dense-matrix OOM)
  tryCatch({
    obj_m <- subset(obj, features = pbmc_markers)

    # Dot plot (ComplexHeatmap object)
    tryCatch({
      p <- FeatureStatPlot(obj_m, features = pbmc_markers,
                           ident = "celltype", plot_type = "dot")
      png_sp(p, "10_markers_dotplot.png", w = 900, h = 600)
      message("Plot: marker dot plot")
    }, error = function(e) message("[skip] dot plot: ", conditionMessage(e)))

    # Heatmap (ComplexHeatmap object)
    tryCatch({
      p <- FeatureStatPlot(obj_m, features = pbmc_markers,
                           ident = "celltype", plot_type = "heatmap",
                           name = "Avg. Expression")
      png_sp(p, "11_markers_heatmap.png", w = 900, h = 600)
      message("Plot: marker heatmap")
    }, error = function(e) message("[skip] heatmap: ", conditionMessage(e)))

    rm(obj_m); invisible(gc())
  }, error = function(e) message("[skip] marker subset: ", conditionMessage(e)))
}

# ---- Steps 15-17: Cell composition plots ---------------------------------
tryCatch({
  p <- CellStatPlot(obj, ident = "celltype", plot_type = "bar", x_text_angle = 90)
  gg_sp(p, "12_celltype_counts.png", w = 9, h = 6)
  message("Plot: cell count bar chart")
}, error = function(e) message("[skip] count bar: ", conditionMessage(e)))

tryCatch({
  p <- CellStatPlot(obj, ident = "celltype", group_by = "condition",
                    frac = "group", swap = TRUE, position = "stack",
                    plot_type = "bar")
  gg_sp(p, "13_composition_by_condition.png", w = 8, h = 6)
  message("Plot: composition by condition")
}, error = function(e) message("[skip] composition bar: ", conditionMessage(e)))

tryCatch({
  p <- CellStatPlot(obj, group_by = "condition", plot_type = "ring",
                    palette = "Spectral")
  gg_sp(p, "14_composition_ring.png", w = 8, h = 6)
  message("Plot: ring chart")
}, error = function(e) message("[skip] ring chart: ", conditionMessage(e)))

tryCatch({
  p <- CellStatPlot(obj, group_by = c("condition", "celltype"),
                    plot_type = "sankey")
  gg_sp(p, "15_composition_sankey.png", w = 10, h = 6)
  message("Plot: Sankey diagram")
}, error = function(e) message("[skip] Sankey: ", conditionMessage(e)))

# ---- Step 19: Marker jitter plot (FindAllMarkers) ------------------------
celltype_markers <- NULL
tryCatch({
  message("\nRunning FindAllMarkers (downsampled 500 per cell type)...")
  Idents(obj) <- "celltype"
  celltype_markers <- FindAllMarkers(obj, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.25,
                                     max.cells.per.ident = 500, verbose = FALSE)
  celltype_markers <- dplyr::rename(celltype_markers, cell_type = cluster)
  message("FindAllMarkers: ", nrow(celltype_markers), " marker rows")
}, error = function(e) message("[skip] FindAllMarkers: ", conditionMessage(e)))

if (!is.null(celltype_markers) && nrow(celltype_markers) > 0) {
  tryCatch({
    p <- MarkersPlot(celltype_markers, subset_by = "cell_type",
                     plot_type = "jitter", select = 5)
    gg_sp(p, "16_markers_jitter.png", w = 14, h = 7)
    message("Plot: marker jitter plot")
  }, error = function(e) message("[skip] marker jitter: ", conditionMessage(e)))
}

# ---- Step 20: Condition DE volcano per cell type -------------------------
tryCatch({
  message(sprintf("\nRunning per-cell-type DE: %s vs %s ...", primary_grp, primary_ref))
  cell_types <- as.character(unique(obj$celltype))

  condition_de <- do.call(rbind, lapply(cell_types, function(ct) {
    de <- tryCatch(
      FindMarkers(obj,
                  ident.1             = primary_grp,
                  ident.2             = primary_ref,
                  group.by            = "condition",
                  subset.ident        = ct,
                  min.pct             = 0.1,
                  logfc.threshold     = 0,
                  max.cells.per.ident = 500,
                  verbose             = FALSE),
      error = function(e) NULL
    )
    if (is.null(de) || nrow(de) == 0) return(NULL)
    de$gene      <- rownames(de)
    de$cell_type <- ct
    de
  }))

  if (!is.null(condition_de) && nrow(condition_de) > 0) {
    message("Condition DE: ", nrow(condition_de), " rows across ", length(cell_types), " cell types")
    tryCatch({
      p <- MarkersPlot(condition_de, subset_by = "cell_type",
                       subset_as_facet = TRUE, plot_type = "volcano",
                       select = 5, ncol = 4)
      gg_sp(p, "17_condition_de_volcano.png", w = 18, h = 12)
      message("Plot: condition DE volcano")
    }, error = function(e) message("[skip] volcano: ", conditionMessage(e)))
  }
}, error = function(e) message("[skip] condition DE: ", conditionMessage(e)))

# ---- Step 22: Cluster tree (ClustreePlot) --------------------------------
tryCatch({
  # Generate clusters at multiple resolutions in-memory (does not overwrite obj)
  if (is.null(res_prefix) || length(res_cols) < 3) {
    message("\nGenerating clusters at multiple resolutions for ClustreePlot...")
    if (!is.null(obj@graphs[["RNA_snn"]])) {
      obj <- FindClusters(obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0),
                          verbose = FALSE)
      res_prefix <- "RNA_snn_res."
    } else {
      message("[skip] No SNN graph found — ClustreePlot skipped")
      res_prefix <- NULL
    }
  }

  if (!is.null(res_prefix)) {
    if (!requireNamespace("clustree", quietly = TRUE))
      stop("clustree not installed")
    p <- ClustreePlot(obj, prefix = res_prefix)
    gg_sp(p, "18_clustree.png", w = 10, h = 9)
    message("Plot: cluster tree (prefix: ", res_prefix, ")")
  }
}, error = function(e) message("[skip] ClustreePlot: ", conditionMessage(e)))

# Restore active identity
Idents(obj) <- "celltype"
log_session(cfg, "spoke_scplotter")
message("\nscplotter done. Plots in: ", dirname(out_path(cfg, "plots", "x")))
