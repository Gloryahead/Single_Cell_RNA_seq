#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# pathway_progeny.R — signaling-pathway activity (decoupleR + PROGENy). (Part 17)
#
# Packages: bioconductor-decoupler, bioconductor-omnipathr
#
# Network: get_progeny(top=500) with CSV cache fallback.
#   PROGENy has 14 pathways with continuous responsiveness weights → run_mlm
#   (not run_ulm; pathways share responsive genes so MLM > ULM here).
#
# Memory strategy (mirrors Part 16):
#   • Full object → pseudobulk (~100 columns) → spoke_pathway_activity.csv
#   • 500 cells/celltype downsample → per-cell MLM → pathwaysmlm assay + plots
#   • rm(obj) before per-cell MLM; pseudobulk built from full obj first
#
# In : results/objects/05_annotated.rds
# Out: results/tables/spoke_pathway_activity.csv           (pseudobulk MLM)
#      results/tables/spoke_pathway_activity_by_celltype.csv
#      results/plots/spoke_progeny_*.png
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(decoupleR); library(Seurat); library(ggplot2); library(patchwork)
  library(dplyr); library(tibble); library(tidyr); library(pheatmap)
  library(RColorBrewer)
})
cfg <- load_config(); set.seed(cfg$seed)

# ---- Step 1: Load object + setup -----------------------------------------
obj <- readRDS(obj_path(cfg, "05_annotated"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
Idents(obj) <- "celltype"

umap_red <- dplyr::case_when(
  "umap"         %in% Reductions(obj) ~ "umap",
  "umap.harmony" %in% Reductions(obj) ~ "umap.harmony",
  "umap.cca"     %in% Reductions(obj) ~ "umap.cca",
  TRUE                                ~ Reductions(obj)[1]
)
message("Using UMAP reduction: ", umap_red)

primary_grp <- cfg$de_contrasts[[1]]$group  # e.g. "Post"
primary_ref <- cfg$de_contrasts[[1]]$ref    # e.g. "Healthy"

# ---- Step 2: Downsample for per-cell MLM (500/celltype) ------------------
obj_sub <- subset(obj, downsample = 500)
cat(sprintf("Downsampled: %d cells (%d celltypes)\n",
            ncol(obj_sub), length(unique(obj_sub$celltype))))

# ---- Step 3: PROGENy network (14 pathways, top=500 genes) ----------------
ref_dir      <- file.path(cfg$project$outdir, "references")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
progeny_csv  <- file.path(ref_dir, "progeny_human.csv")

net <- tryCatch({
  message("Step 3: Downloading PROGENy from OmniPath (top=500)...")
  n <- decoupleR::get_progeny(organism = "human", top = 500)
  write.csv(n, progeny_csv, row.names = FALSE)
  message(sprintf("PROGENy: %d interactions, %d pathways",
                  nrow(n), length(unique(n$source))))
  n
}, error = function(e) {
  message("OmniPath download failed: ", conditionMessage(e))
  if (file.exists(progeny_csv)) {
    message("Loading cached PROGENy from: ", progeny_csv)
    read.csv(progeny_csv, stringsAsFactors = FALSE)
  } else {
    stop("PROGENy download failed and no cache found at: ", progeny_csv)
  }
})

# ---- Step 4a: Pseudobulk MLM on full object (main spoke output) ----------
# Pseudobulk is the scientifically preferred PROGENy usage (Schubert 2018).
# Mean log-norm expression per celltype×sample keeps per-cell MLM as supplement.
message("\nStep 4a: Pseudobulk MLM on full object...")
obj$pb_id <- paste(obj$celltype, obj$sample, sep = "..")
groups    <- split(seq_len(ncol(obj)), obj$pb_id)
pb_mat    <- do.call(cbind, lapply(groups, function(idx)
  Matrix::rowMeans(
    GetAssayData(obj, assay = "RNA", layer = "data")[, idx, drop = FALSE])))
colnames(pb_mat) <- names(groups)
cat(sprintf("Pseudobulk: %d genes × %d columns\n", nrow(pb_mat), ncol(pb_mat)))

rm(obj); invisible(gc())

acts_pb <- decoupleR::run_mlm(mat = pb_mat, net = net,
                              .source = "source", .target = "target", .mor = "weight",
                              minsize = 5)
write.csv(acts_pb, out_path(cfg, "tables", "spoke_pathway_activity.csv"), row.names = FALSE)
message(sprintf("Pseudobulk MLM: %d rows → spoke_pathway_activity.csv", nrow(acts_pb)))
rm(pb_mat, acts_pb); invisible(gc())

# ---- Step 4b: Per-cell MLM on downsampled object -------------------------
# 14 pathways × top-500 genes; dense matrix for ~4K cells is ~0.6 GB
message("\nStep 4b: Per-cell MLM on downsampled object...")
mat_sub <- as.matrix(GetAssayData(obj_sub, assay = "RNA", layer = "data"))
cat(sprintf("Dense matrix: %d genes × %d cells\n", nrow(mat_sub), ncol(mat_sub)))

acts <- decoupleR::run_mlm(mat = mat_sub, net = net,
                           .source = "source", .target = "target", .mor = "weight",
                           minsize = 5)
rm(mat_sub); invisible(gc())
message(sprintf("Per-cell MLM: %d rows (pathway × cell)", nrow(acts)))

# ---- Step 5: Store activities as pathwaysmlm assay -----------------------
pw_mat <- acts %>%
  filter(statistic == "mlm") %>%
  tidyr::pivot_wider(id_cols = "source",
                     names_from  = "condition",  # "condition" = cell barcode
                     values_from = "score") %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

shared_bcs <- intersect(colnames(pw_mat), colnames(obj_sub))
if (length(shared_bcs) < 10)
  stop("Too few shared barcodes between MLM output and Seurat object: ", length(shared_bcs))
pw_mat  <- pw_mat[, shared_bcs]
obj_sub <- obj_sub[, shared_bcs]

obj_sub[["pathwaysmlm"]] <- CreateAssayObject(data = pw_mat)
DefaultAssay(obj_sub) <- "pathwaysmlm"
obj_sub <- ScaleData(obj_sub, verbose = FALSE)
LayerData(obj_sub, assay = "pathwaysmlm", layer = "data") <-
  LayerData(obj_sub, assay = "pathwaysmlm", layer = "scale.data")
message("pathwaysmlm assay stored; pathways: ", paste(rownames(obj_sub), collapse = ", "))

# Diverging palette (blue = low/suppressed, red = high/active)
div_pal <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
pal100  <- grDevices::colorRampPalette(div_pal)(100)

# ---- Step 6: Mean pathway activity heatmap by cell type ------------------
tryCatch({
  act_data <- as.matrix(LayerData(obj_sub, assay = "pathwaysmlm", layer = "data"))

  pathway_long <- t(act_data) %>%
    as.data.frame() %>%
    dplyr::mutate(cell_type = obj_sub$celltype) %>%
    tidyr::pivot_longer(cols = -cell_type, names_to = "pathway", values_to = "score") %>%
    dplyr::group_by(cell_type, pathway) %>%
    dplyr::summarise(mean_activity = mean(score, na.rm = TRUE), .groups = "drop")

  # Export cell-type summary
  write.csv(pathway_long,
            out_path(cfg, "tables", "spoke_pathway_activity_by_celltype.csv"),
            row.names = FALSE)

  heatmap_mat <- pathway_long %>%
    tidyr::pivot_wider(id_cols = "cell_type", names_from = "pathway",
                       values_from = "mean_activity") %>%
    tibble::column_to_rownames("cell_type") %>%
    as.matrix()

  pheatmap::pheatmap(
    mat          = heatmap_mat,
    color        = pal100,
    border_color = "white",
    cellwidth    = 18, cellheight = 14,
    treeheight_row = 20, treeheight_col = 20,
    main         = "PROGENy pathway activity by cell type (scaled MLM scores)",
    filename     = out_path(cfg, "plots", "spoke_progeny_01_pathway_by_celltype.png")
  )
  message(sprintf("Plot: pathway-by-celltype heatmap (%d pathways × %d cell types)",
                  ncol(heatmap_mat), nrow(heatmap_mat)))
}, error = function(e) message("[skip] celltype heatmap: ", conditionMessage(e)))

# ---- Step 7: UMAP — cell types + showcase pathway activity ---------------
# Auto-pick a relevant showcase pathway; prefer inflammation-related ones
showcase_pw <- intersect(c("TNFa", "NFkB", "JAK-STAT", "MAPK", "Hypoxia"),
                         rownames(obj_sub[["pathwaysmlm"]]))[1]

tryCatch({
  if (!is.null(showcase_pw) && !is.na(showcase_pw)) {
    p_ct <- DimPlot(obj_sub, reduction = umap_red, group.by = "celltype",
                    label = TRUE, repel = TRUE, pt.size = 0.3) +
            NoLegend() + ggtitle("Cell types")

    DefaultAssay(obj_sub) <- "pathwaysmlm"
    p_pw <- FeaturePlot(obj_sub, features = showcase_pw, reduction = umap_red,
                        pt.size = 0.3) +
            scale_colour_gradient2(low = div_pal[1], mid = "white",
                                   high = div_pal[11]) +
            ggtitle(paste0(showcase_pw, " activity"))

    ggsave(out_path(cfg, "plots",
                    paste0("spoke_progeny_02_", gsub("-", "", showcase_pw), "_umap.png")),
           p_ct | p_pw, width = 14, height = 6, dpi = 150)
    message("Plot: UMAP + ", showcase_pw, " activity")
  }
}, error = function(e) message("[skip] UMAP plot: ", conditionMessage(e)))

# ---- Step 8: Pathway activity by cell type × condition -------------------
tryCatch({
  # Inflammatory / disease-relevant PROGENy pathways
  condition_pws <- c("TNFa", "NFkB", "JAK-STAT", "MAPK", "Trail", "p53")
  condition_pws <- intersect(condition_pws, rownames(obj_sub[["pathwaysmlm"]]))

  if (length(condition_pws) >= 2) {
    DefaultAssay(obj_sub) <- "pathwaysmlm"
    act_data <- as.matrix(LayerData(obj_sub, assay = "pathwaysmlm", layer = "data"))

    cond_long <- t(act_data) %>%
      as.data.frame() %>%
      dplyr::mutate(cell_type = obj_sub$celltype,
                    condition = obj_sub$condition) %>%
      tidyr::pivot_longer(cols = c(-cell_type, -condition),
                          names_to = "pathway", values_to = "score") %>%
      dplyr::filter(pathway %in% condition_pws) %>%
      dplyr::group_by(cell_type, condition, pathway) %>%
      dplyr::summarise(mean_activity = mean(score, na.rm = TRUE), .groups = "drop")

    cond_mat <- cond_long %>%
      dplyr::mutate(group = paste(cell_type, condition, sep = " :: ")) %>%
      tidyr::pivot_wider(id_cols = "group", names_from = "pathway",
                         values_from = "mean_activity") %>%
      tibble::column_to_rownames("group") %>%
      as.matrix()

    pheatmap::pheatmap(
      mat          = cond_mat,
      color        = pal100,
      border_color = "white",
      cluster_rows = FALSE,
      cellwidth    = 18, cellheight = 14,
      main         = paste0("Inflammatory pathway activity: cell type × condition\n",
                            "(", primary_grp, " vs ", primary_ref, ")"),
      filename     = out_path(cfg, "plots", "spoke_progeny_03_pathway_by_condition.png")
    )
    message(sprintf("Plot: condition heatmap (%d pathways)", length(condition_pws)))
  } else {
    message("[skip] condition heatmap: fewer than 2 condition pathways in assay")
  }
}, error = function(e) message("[skip] condition heatmap: ", conditionMessage(e)))

# ---- Step 9: Export + session log ----------------------------------------
DefaultAssay(obj_sub) <- "RNA"
log_session(cfg, "spoke_progeny")
message("\nprogeny done.")
message("  → spoke_pathway_activity.csv           (pseudobulk MLM, full dataset)")
message("  → spoke_pathway_activity_by_celltype.csv (per-celltype mean, downsampled)")
