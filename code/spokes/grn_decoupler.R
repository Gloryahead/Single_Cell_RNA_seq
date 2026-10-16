#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# grn_decoupler.R — TF activity (decoupleR + CollecTRI). (Part 16)
#
# Packages: bioconductor-decoupler, bioconductor-omnipathr, msigdbr (fallback)
#
# Network priority:
#   1. get_collectri() live download → cached to results/references/collectri_human.csv
#   2. Read from cache (offline / OmniPath down)
#   3. GTRD via msigdbr (fallback; no mor=-1 repression sign)
# Prior script used GTRD because OmnipathR::unnest_evidences crashed (container bug).
# New script tries CollecTRI first; GTRD is retained as final fallback.
#
# Memory strategy:
#   • Full object → pseudobulk (~100 columns) for the main spoke_tf_activity.csv
#   • 500 cells/celltype downsample for per-cell ULM + all visualizations
#   • rm(obj) before per-cell ULM; pseudobulk built from full obj first
#
# In : results/objects/05_annotated.rds
# Out: results/tables/spoke_tf_activity.csv           (pseudobulk ULM — main output)
#      results/tables/spoke_tf_activity_by_celltype.csv (mean activity per celltype)
#      results/plots/spoke_decoupler_*.{png,pdf}
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(decoupleR); library(Seurat); library(ggplot2); library(patchwork)
  library(dplyr); library(tibble); library(tidyr); library(pheatmap)
  library(RColorBrewer)
})
cfg <- load_config(); set.seed(cfg$seed)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# ---- Step 1: Load object + setup -----------------------------------------
obj <- readRDS(obj_path(cfg, "05_annotated"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
Idents(obj) <- "celltype"

# Detect UMAP reduction
umap_red <- dplyr::case_when(
  "umap"         %in% Reductions(obj) ~ "umap",
  "umap.harmony" %in% Reductions(obj) ~ "umap.harmony",
  "umap.cca"     %in% Reductions(obj) ~ "umap.cca",
  TRUE                                ~ Reductions(obj)[1]
)
message("Using UMAP reduction: ", umap_red)

# Primary contrast from config
primary_grp <- cfg$de_contrasts[[1]]$group  # e.g. "Post"
primary_ref <- cfg$de_contrasts[[1]]$ref    # e.g. "Healthy"

# ---- Step 2: Downsample for per-cell ULM (500/celltype) ------------------
# Stratified downsample to keep every celltype and condition represented.
# per-cell ULM on 106K cells → ~17 GB dense matrix → OOM; ~4K cells → ~0.6 GB.
obj_sub <- subset(obj, downsample = 500)
cat(sprintf("Downsampled: %d cells (%d celltypes)\n",
            ncol(obj_sub), length(unique(obj_sub$celltype))))

# ---- Step 3: CollecTRI network -------------------------------------------
ref_dir       <- file.path(cfg$project$outdir, "references")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
collectri_csv <- file.path(ref_dir, "collectri_human.csv")

net <- tryCatch({
  message("Step 3: Downloading CollecTRI from OmniPath...")
  n <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
  write.csv(n, collectri_csv, row.names = FALSE)
  message(sprintf("CollecTRI downloaded: %d interactions, %d TFs",
                  nrow(n), length(unique(n$source))))
  n
}, error = function(e) {
  message("OmniPath download failed: ", conditionMessage(e))
  if (file.exists(collectri_csv)) {
    message("Loading cached CollecTRI from: ", collectri_csv)
    read.csv(collectri_csv, stringsAsFactors = FALSE)
  } else {
    message("No cache found — falling back to GTRD (MSigDB C3, no repression sign).")
    NULL
  }
})

if (is.null(net)) {
  # GTRD fallback: no -1 MoR, but no OmniPath dependency
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("msigdbr not installed and CollecTRI download failed. Cannot build network.")
  sp    <- if (isTRUE(cfg$project$organism == "mouse")) "Mus musculus" else "Homo sapiens"
  tf_df <- msigdbr::msigdbr(species = sp, collection = "C3", subcollection = "TFT:GTRD")
  net   <- data.frame(source = tf_df$gs_name, target = tf_df$gene_symbol,
                      mor = 1L, stringsAsFactors = FALSE)
  net   <- dplyr::distinct(net, source, target, .keep_all = TRUE)
  message(sprintf("GTRD fallback: %d interactions, %d TFs", nrow(net), length(unique(net$source))))
}

# ---- Step 4a: Pseudobulk ULM on full object (main spoke output) ----------
# Mean log-norm expression per celltype×sample (~100 pseudobulk columns).
# This is the comprehensive, statistically robust output for downstream analyses.
message("\nStep 4a: Building pseudobulk matrix from full object...")
obj$pb_id <- paste(obj$celltype, obj$sample, sep = "..")
groups    <- split(seq_len(ncol(obj)), obj$pb_id)
pb_mat    <- do.call(cbind, lapply(groups, function(idx)
  Matrix::rowMeans(
    GetAssayData(obj, assay = "RNA", layer = "data")[, idx, drop = FALSE])))
colnames(pb_mat) <- names(groups)
cat(sprintf("Pseudobulk matrix: %d genes × %d columns\n", nrow(pb_mat), ncol(pb_mat)))

# Free full object — pseudobulk mat is small; obj_sub still in memory for viz
rm(obj); invisible(gc())

acts_pb <- decoupleR::run_ulm(mat = pb_mat, net = net,
                              .source = "source", .target = "target", .mor = "mor",
                              minsize = 5)
write.csv(acts_pb, out_path(cfg, "tables", "spoke_tf_activity.csv"), row.names = FALSE)
message(sprintf("Pseudobulk ULM: %d rows → spoke_tf_activity.csv", nrow(acts_pb)))
rm(pb_mat, acts_pb); invisible(gc())

# ---- Step 4b: Per-cell ULM on downsampled object -------------------------
# Dense matrix is feasible at ~4K cells × 20K genes (~0.6 GB).
message("\nStep 4b: Per-cell ULM on downsampled object...")
mat_sub <- as.matrix(GetAssayData(obj_sub, assay = "RNA", layer = "data"))
cat(sprintf("Dense matrix: %d genes × %d cells\n", nrow(mat_sub), ncol(mat_sub)))

acts <- decoupleR::run_ulm(mat = mat_sub, net = net,
                           .source = "source", .target = "target", .mor = "mor",
                           minsize = 5)
rm(mat_sub); invisible(gc())
message(sprintf("Per-cell ULM: %d rows (TF × cell)", nrow(acts)))

# ---- Step 5: Store activities as a new Seurat assay ---------------------
tf_mat <- acts %>%
  filter(statistic == "ulm") %>%
  tidyr::pivot_wider(id_cols = "source",
                     names_from  = "condition",  # "condition" = cell barcode in decoupleR
                     values_from = "score") %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Align columns to obj_sub barcodes (decoupleR may reorder or drop some)
shared_bcs <- intersect(colnames(tf_mat), colnames(obj_sub))
if (length(shared_bcs) < 10)
  stop("Too few shared barcodes between ULM output and Seurat object: ", length(shared_bcs))
tf_mat  <- tf_mat[, shared_bcs]
obj_sub <- obj_sub[, shared_bcs]

obj_sub[["tfsulm"]] <- CreateAssayObject(data = tf_mat)
DefaultAssay(obj_sub) <- "tfsulm"
# Scale so high- and low-magnitude TFs are visually comparable in heatmaps
obj_sub <- ScaleData(obj_sub, verbose = FALSE)
# Replace data layer with scaled values for plotting
LayerData(obj_sub, assay = "tfsulm", layer = "data") <-
  LayerData(obj_sub, assay = "tfsulm", layer = "scale.data")
message("tfsulm assay stored in downsampled object")

# Diverging colour helper (blue = inactive, white = 0, red = active)
div_pal <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))

# ---- Step 6: Activity vs expression comparison UMAP ----------------------
# Showcase TF: prefer PAX5 (B cells), fallback through known PBMC regulators
showcase_tf <- intersect(c("PAX5", "SPI1", "CEBPB", "RELA", "NFKB1", "EBF1"),
                         rownames(obj_sub[["tfsulm"]]))[1]

if (!is.null(showcase_tf) && !is.na(showcase_tf)) {
  tryCatch({
    p_ct <- DimPlot(obj_sub, reduction = umap_red, group.by = "celltype",
                    label = TRUE, repel = TRUE, pt.size = 0.3) +
            NoLegend() + ggtitle("Cell types")

    DefaultAssay(obj_sub) <- "tfsulm"
    p_act <- FeaturePlot(obj_sub, features = showcase_tf, reduction = umap_red,
                         pt.size = 0.3) +
             scale_colour_gradient2(low = div_pal[1], mid = "white",
                                    high = div_pal[11]) +
             ggtitle(paste0(showcase_tf, " activity"))

    DefaultAssay(obj_sub) <- "RNA"
    p_exp <- FeaturePlot(obj_sub, features = showcase_tf, reduction = umap_red,
                         pt.size = 0.3) +
             ggtitle(paste0(showcase_tf, " expression"))
    DefaultAssay(obj_sub) <- "tfsulm"

    p_cmp <- p_ct | p_act | p_exp
    ggsave(out_path(cfg, "plots", "spoke_decoupler_activity_vs_expression.png"),
           p_cmp, width = 18, height = 6, dpi = 150)
    message("Plot: activity vs expression (", showcase_tf, ")")
  }, error = function(e) message("[skip] activity vs expression: ", conditionMessage(e)))
}

# ---- Step 7: Mean activity heatmap of top variable TFs -------------------
tryCatch({
  n_tfs <- 25
  act_data <- as.matrix(LayerData(obj_sub, assay = "tfsulm", layer = "data"))

  activity_long <- t(act_data) %>%
    as.data.frame() %>%
    dplyr::mutate(cell_type = obj_sub$celltype) %>%
    tidyr::pivot_longer(cols = -cell_type,
                        names_to  = "source",
                        values_to = "score") %>%
    dplyr::group_by(cell_type, source) %>%
    dplyr::summarise(mean_activity = mean(score, na.rm = TRUE), .groups = "drop")

  top_tfs <- activity_long %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(variability = sd(mean_activity, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(variability)) %>%
    head(n_tfs) %>%
    dplyr::pull(source)

  heatmap_mat <- activity_long %>%
    dplyr::filter(source %in% top_tfs) %>%
    tidyr::pivot_wider(id_cols = "cell_type", names_from = "source",
                       values_from = "mean_activity") %>%
    tibble::column_to_rownames("cell_type") %>%
    as.matrix()

  # Export cell-type activity summary
  write.csv(activity_long,
            out_path(cfg, "tables", "spoke_tf_activity_by_celltype.csv"),
            row.names = FALSE)

  pal100 <- grDevices::colorRampPalette(div_pal)(100)
  pheatmap::pheatmap(
    mat          = heatmap_mat,
    color        = pal100,
    border_color = "white",
    cellwidth    = 14, cellheight = 14,
    treeheight_row = 20, treeheight_col = 20,
    main         = "Top variable TF activity by cell type",
    filename     = out_path(cfg, "plots", "spoke_decoupler_top_tf_heatmap.png")
  )
  message(sprintf("Plot: top-TF heatmap (%d TFs × %d cell types)", n_tfs, nrow(heatmap_mat)))
}, error = function(e) message("[skip] top-TF heatmap: ", conditionMessage(e)))

# ---- Step 8: Inflammatory TF activity by celltype × condition heatmap ---
tryCatch({
  # Canonical inflammatory regulators; filter to those present in the assay
  candidate_tfs <- c("NFKB1", "RELA", "STAT1", "STAT3", "IRF1", "JUN", "FOS",
                     "NFKB2", "RELB", "IRF4", "IRF8", "CEBPB")
  candidate_tfs <- intersect(candidate_tfs, rownames(obj_sub[["tfsulm"]]))

  if (length(candidate_tfs) >= 2) {
    act_data <- as.matrix(LayerData(obj_sub, assay = "tfsulm", layer = "data"))
    cond_long <- t(act_data) %>%
      as.data.frame() %>%
      dplyr::mutate(cell_type = obj_sub$celltype,
                    condition = obj_sub$condition) %>%
      tidyr::pivot_longer(cols = c(-cell_type, -condition),
                          names_to  = "source",
                          values_to = "score") %>%
      dplyr::filter(source %in% candidate_tfs) %>%
      dplyr::group_by(cell_type, condition, source) %>%
      dplyr::summarise(mean_activity = mean(score, na.rm = TRUE), .groups = "drop")

    cond_mat <- cond_long %>%
      dplyr::mutate(group = paste(cell_type, condition, sep = " :: ")) %>%
      tidyr::pivot_wider(id_cols = "group", names_from = "source",
                         values_from = "mean_activity") %>%
      tibble::column_to_rownames("group") %>%
      as.matrix()

    pal100 <- grDevices::colorRampPalette(div_pal)(100)
    pheatmap::pheatmap(
      mat          = cond_mat,
      color        = pal100,
      border_color = "white",
      cluster_rows = FALSE,     # keep celltype::condition groups readable
      cellwidth    = 16, cellheight = 14,
      main         = paste0("Inflammatory TF activity: cell type × condition\n",
                            "(", primary_grp, " vs ", primary_ref, ")"),
      filename     = out_path(cfg, "plots", "spoke_decoupler_condition_heatmap.png")
    )
    message(sprintf("Plot: condition heatmap (%d TFs)", length(candidate_tfs)))
  } else {
    message("[skip] condition heatmap: fewer than 2 candidate TFs found in assay")
  }
}, error = function(e) message("[skip] condition heatmap: ", conditionMessage(e)))

# ---- Step 9: Single-TF UMAP split by condition ---------------------------
tryCatch({
  tf_split <- if (!is.null(showcase_tf) && !is.na(showcase_tf)) showcase_tf else "NFKB1"
  tf_split  <- intersect(tf_split, rownames(obj_sub[["tfsulm"]]))[1]
  if (!is.na(tf_split)) {
    DefaultAssay(obj_sub) <- "tfsulm"
    p_split <- FeaturePlot(obj_sub, features = tf_split,
                           reduction = umap_red, split.by = "condition",
                           pt.size = 0.3) &
               scale_colour_gradient2(low = div_pal[1], mid = "white",
                                      high = div_pal[11])
    ggsave(out_path(cfg, "plots",
                    paste0("spoke_decoupler_", tf_split, "_by_condition.png")),
           p_split, width = 5 * length(unique(obj_sub$condition)), height = 6, dpi = 150)
    message("Plot: single-TF UMAP by condition (", tf_split, ")")
  }
}, error = function(e) message("[skip] single-TF UMAP: ", conditionMessage(e)))

# ---- Step 10: Export + save session --------------------------------------
DefaultAssay(obj_sub) <- "RNA"
log_session(cfg, "spoke_decoupler_tf")
message("\ndecoupler done.")
message("  → spoke_tf_activity.csv        (pseudobulk ULM, full dataset)")
message("  → spoke_tf_activity_by_celltype.csv (per-celltype mean, downsampled)")
