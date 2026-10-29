#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# trajectory_slingshot.R — T/NK trajectory & pseudotime (Slingshot). (Part 7-2)
#
# Strategy: analyse the T/NK lineage only — running Slingshot on all PBMC
# cell types creates hub nodes and spurious lineages.
#
# In : results/objects/05_annotated.rds
# Out: results/objects/spoke_slingshot.rds   (SCE with slingshot metadata)
#      results/tables/spoke_slingshot_pseudotime_genes.csv
#      results/plots/spoke_slingshot_*.{png,pdf}
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment); library(patchwork)
})
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# ---- T/NK subset --------------------------------------------------------
# Detect T/NK cell types by label pattern (works with HPCA, Monaco, manual)
tnk_types <- grep("T_cell|NK_cell|T cell|NK cell|CD4|CD8",
                  unique(obj$celltype), value = TRUE, ignore.case = TRUE)
if (length(tnk_types) == 0) {
  message("WARNING: no T/NK cell types found — running on all cell types")
  tnk_types <- unique(obj$celltype)
}
message("T/NK trajectory cell types: ", paste(tnk_types, collapse = ", "))
sub <- obj[, obj$celltype %in% tnk_types]
cat(sprintf("T/NK subset: %d cells across %d cell types\n",
            ncol(sub), length(tnk_types)))

# ---- Seurat → SCE, ensure UMAP slot is available -----------------------
sce <- as.SingleCellExperiment(sub)
if (!"UMAP" %in% reducedDimNames(sce)) {
  umap_hit <- grep("UMAP", reducedDimNames(sce), ignore.case = TRUE)[1]
  if (!is.na(umap_hit))
    reducedDims(sce)[["UMAP"]] <- reducedDim(sce, reducedDimNames(sce)[umap_hit])
  else
    stop("No UMAP reduction found — run RunUMAP() before this spoke.")
}

# ---- Slingshot trajectory inference ------------------------------------
# start.clus = NULL: auto-detect start cluster from topology.
# allow.breaks = FALSE prevents disconnected curves.
# extend = "n" stops curves from extrapolating beyond observed cells.
sce <- slingshot(sce,
                 clusterLabels = "celltype",
                 reducedDim    = "UMAP",
                 start.clus    = NULL,
                 allow.breaks  = FALSE,
                 extend        = "n",
                 reweight      = TRUE,
                 reassign      = TRUE)

lineages <- slingLineages(sce)
n_lin    <- length(lineages)
cat(sprintf("\nSlingshot detected %d lineage(s):\n", n_lin))
for (i in seq_along(lineages))
  cat(sprintf("  Lineage %d: %s\n", i, paste(lineages[[i]], collapse = " -> ")))

# ---- Assign primary lineage (highest weight) + pseudotime --------------
pt_mat  <- slingPseudotime(sce, na = FALSE)
wt_mat  <- slingCurveWeights(sce)
pri_idx <- apply(wt_mat, 1, which.max)
pri_pt  <- pt_mat[cbind(seq_len(nrow(pt_mat)), pri_idx)]

sce$slingshot_pseudotime <- pri_pt
sce$slingshot_lineage    <- paste0("Lineage_", pri_idx)

cat("\nMedian pseudotime by cell type:\n")
print(round(tapply(sce$slingshot_pseudotime, sce$celltype, median), 2))

# ---- Marker validation: naive markers should decrease, effector increase
val_markers  <- c("CCR7", "LEF1", "TCF7", "CD69", "GZMB", "PRF1", "NKG7")
val_present  <- val_markers[val_markers %in% rownames(sce)]
if (length(val_present) > 0) {
  marker_mat <- as.matrix(logcounts(sce)[val_present, , drop = FALSE])
  cor_val    <- sapply(val_present, function(g)
    cor(marker_mat[g, ], sce$slingshot_pseudotime, method = "spearman"))
  cat("\nMarker-pseudotime Spearman r (CCR7/LEF1 should be <0, GZMB/NKG7 >0):\n")
  print(round(cor_val, 3))
}

# ---- Condition comparison -----------------------------------------------
conds <- split(sce$slingshot_pseudotime, sce$condition)
if (length(conds) >= 2) {
  wt <- tryCatch(wilcox.test(conds[[1]], conds[[2]]), error = function(e) NULL)
  cat("\nPseudotime by condition (Wilcoxon):\n")
  cat(sprintf("  p = %.4g\n", if (!is.null(wt)) wt$p.value else NA))
  for (nm in names(conds))
    cat(sprintf("  Median %s: %.2f\n", nm, median(conds[[nm]])))
}

# ---- Gene-pseudotime correlations (Spearman, BH) -----------------------
# Compute t-approximation instead of per-gene cor.test (much faster).
message("Computing gene-pseudotime correlations (Spearman)...")
# Filter on sparse matrix to avoid materializing a 83K-cell × 33K-gene dense matrix (~21 GB)
g_keep     <- rowSums(logcounts(sce) > 0) >= 10
gene_keep  <- rownames(sce)[g_keep]
pt_vec     <- sce$slingshot_pseudotime
n_c        <- ncol(sce)
# Process 2000 genes at a time so peak memory per chunk is ~1.5 GB
chunk_size  <- 2000L
gene_chunks <- split(gene_keep, ceiling(seq_along(gene_keep) / chunk_size))
rho_vec <- unlist(lapply(gene_chunks, function(gns) {
  mat_chunk <- as.matrix(logcounts(sce)[gns, , drop = FALSE])
  apply(mat_chunk, 1, function(x) cor(x, pt_vec, method = "spearman"))
}))
invisible(gc())
t_stat  <- rho_vec * sqrt((n_c - 2) / pmax(1e-12, 1 - rho_vec^2))
p_vec   <- 2 * pt(abs(t_stat), df = n_c - 2, lower.tail = FALSE)
fdr_vec <- p.adjust(p_vec, method = "BH")
gene_df <- data.frame(gene = names(rho_vec), rho = round(rho_vec, 4),
                      p_value = p_vec, FDR = fdr_vec, stringsAsFactors = FALSE)
gene_df <- gene_df[order(gene_df$FDR, -abs(gene_df$rho)), ]
n_sig   <- sum(gene_df$FDR < 0.05, na.rm = TRUE)
cat(sprintf("\n%d trajectory-associated genes (FDR < 0.05)\n", n_sig))
cat("Top genes:\n"); print(head(gene_df[gene_df$FDR < 0.05, ], 10))
write.csv(gene_df[gene_df$FDR < 0.05, ],
          out_path(cfg, "tables", "spoke_slingshot_pseudotime_genes.csv"),
          row.names = FALSE)
invisible(gc())

# ---- Plots --------------------------------------------------------------
umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype   <- sce$celltype
umap_df$pseudotime <- sce$slingshot_pseudotime
umap_df$lineage    <- sce$slingshot_lineage

# Smooth curve coordinates per lineage
curves   <- slingCurves(sce)
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  pts <- curves[[i]]$s[curves[[i]]$ord, , drop = FALSE]
  data.frame(UMAP_1 = pts[, 1], UMAP_2 = pts[, 2],
             lineage = paste0("Lineage_", i))
}))

ct_pal <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                     "#FF7F00","#A65628","#F781BF","#999999"))(length(tnk_types)),
  tnk_types)

p_traj <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP_1, UMAP_2, color = celltype)) +
  ggplot2::geom_point(size = 0.4, alpha = 0.5) +
  ggplot2::geom_path(data = curve_df,
                     ggplot2::aes(group = lineage, color = NULL),
                     color = "black", linewidth = 1.5) +
  ggplot2::scale_color_manual(values = ct_pal) +
  ggplot2::labs(title = "Slingshot T/NK trajectories", color = NULL) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))

p_pt <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP_1, UMAP_2, color = pseudotime)) +
  ggplot2::geom_point(size = 0.4, alpha = 0.7) +
  ggplot2::scale_color_viridis_c(option = "plasma") +
  ggplot2::labs(title = "Pseudotime") +
  ggplot2::theme_classic(base_size = 12)

p_lin <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP_1, UMAP_2, color = lineage)) +
  ggplot2::geom_point(size = 0.4, alpha = 0.7) +
  ggplot2::labs(title = "Primary lineage") +
  ggplot2::theme_classic(base_size = 12)

ggsave(out_path(cfg, "plots", "spoke_slingshot_overview.png"),
       wrap_plots(p_traj, p_pt, p_lin, ncol = 3), width = 18, height = 6)
ggsave(out_path(cfg, "plots", "spoke_slingshot_overview.pdf"),
       wrap_plots(p_traj, p_pt, p_lin, ncol = 3), width = 18, height = 6)

# Pseudotime violin by cell type
vln_df <- data.frame(celltype = sce$celltype, pseudotime = sce$slingshot_pseudotime)
p_vln  <- ggplot2::ggplot(vln_df, ggplot2::aes(celltype, pseudotime, fill = celltype)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE) +
  ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  ggplot2::scale_fill_manual(values = ct_pal) +
  ggplot2::labs(x = NULL, y = "Pseudotime", title = "Pseudotime by cell type") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 plot.title  = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "spoke_slingshot_pseudotime_violin.png"),
       p_vln, width = 8, height = 5)
ggsave(out_path(cfg, "plots", "spoke_slingshot_pseudotime_violin.pdf"),
       p_vln, width = 8, height = 5)

# Pseudotime density by condition per lineage
cond_df <- data.frame(condition  = sce$condition,
                      pseudotime = sce$slingshot_pseudotime,
                      lineage    = sce$slingshot_lineage)
p_dens  <- ggplot2::ggplot(cond_df, ggplot2::aes(pseudotime, fill = condition)) +
  ggplot2::geom_density(alpha = 0.55) +
  ggplot2::facet_wrap(~ lineage, scales = "free_y") +
  ggplot2::labs(x = "Pseudotime", y = "Density",
                title = "Pseudotime distribution by condition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "spoke_slingshot_pseudotime_condition.png"),
       p_dens, width = 10, height = 5)
ggsave(out_path(cfg, "plots", "spoke_slingshot_pseudotime_condition.pdf"),
       p_dens, width = 10, height = 5)
rm(p_traj, p_pt, p_lin, p_vln, p_dens); invisible(gc())

# ---- Add pseudotime back to the full Seurat object and plot ------------
obj@meta.data[["slingshot_pseudotime"]] <- NA_real_
obj@meta.data[["slingshot_lineage"]]    <- NA_character_
m <- match(colnames(sub), colnames(obj))
obj@meta.data[["slingshot_pseudotime"]][m] <- sce$slingshot_pseudotime
obj@meta.data[["slingshot_lineage"]][m]    <- sce$slingshot_lineage

p_feat <- FeaturePlot(obj, "slingshot_pseudotime", pt.size = 0.3, order = TRUE) +
          ggplot2::scale_color_viridis_c(option = "plasma", na.value = "grey90") +
          ggplot2::ggtitle("Slingshot pseudotime (T/NK)")
p_dmp  <- DimPlot(obj, group.by = "slingshot_lineage", pt.size = 0.3) +
          ggplot2::ggtitle("Slingshot lineage (T/NK)")
ggsave(out_path(cfg, "plots", "spoke_slingshot_in_seurat.png"),
       p_feat + p_dmp, width = 14, height = 6)
ggsave(out_path(cfg, "plots", "spoke_slingshot_in_seurat.pdf"),
       p_feat + p_dmp, width = 14, height = 6)

saveRDS(sce, out_path(cfg, "objects", "spoke_slingshot.rds"))
log_session(cfg, "spoke_slingshot")
