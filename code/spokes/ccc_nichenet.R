#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# ccc_nichenet.R — ligand-target prediction (NicheNet).          (Part 10)
#
# NOTE: nichenetr is NOT in the current scrna_r.sif. Rebuild with:
#   devtools::install_github("saeyslab/nichenetr")
#   BiocManager::install("ComplexHeatmap")
#   install.packages(c("tidyverse","circlize","ggrepel","igraph","ggraph"))
#
# Prior models are downloaded from Zenodo on the first run and cached to
# results/references/nichenet/ (on xdisk, ~1.2 GB total).
# Set options(timeout=600) before this script if network is slow.
#
# In : results/objects/05_annotated.rds
# Out: results/tables/spoke_nichenet_ligands.csv
#      results/tables/spoke_nichenet_prioritized_lr.csv
#      results/tables/spoke_nichenet_sender_summary.csv
#      results/plots/nichenet/
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(nichenetr); library(Seurat); library(tidyverse)
  library(ComplexHeatmap); library(circlize)
  library(patchwork); library(ggrepel)
})
options(timeout = 600)
cfg <- load_config(); set.seed(cfg$seed)
sp  <- if (cfg$project$organism == "mouse") "mouse" else "human"

# ---- Step 1: Prior models (download once, cache to xdisk) ---------------
nn_dir   <- file.path(cfg$project$outdir, "references", "nichenet")
dir.create(nn_dir, recursive = TRUE, showWarnings = FALSE)

zenodo   <- "https://zenodo.org/record/7074291/files"
nn_files <- list(
  lr_network    = "lr_network_human_21122021.rds",
  ligand_target = "ligand_target_matrix_nsga2r_final.rds",
  weighted_nets = "weighted_networks_nsga2r_final.rds",
  ligand_tf     = "ligand_tf_matrix_nsga2r_final.rds"
)
if (sp == "mouse") {
  nn_files$lr_network    <- "lr_network_mouse_21122021.rds"
  nn_files$ligand_target <- "ligand_target_matrix_nsga2r_final_mouse.rds"
  nn_files$weighted_nets <- "weighted_networks_nsga2r_final_mouse.rds"
  nn_files$ligand_tf     <- "ligand_tf_matrix_nsga2r_final_mouse.rds"
}
for (nm in names(nn_files)) {
  lp <- file.path(nn_dir, nn_files[[nm]])
  if (!file.exists(lp)) {
    message("Downloading ", nn_files[[nm]], " ...")
    download.file(paste0(zenodo, "/", nn_files[[nm]]), lp, mode = "wb", quiet = TRUE)
  }
}
lr_network    <- readRDS(file.path(nn_dir, nn_files$lr_network)) %>%
                   distinct(from, to)
ligand_target <- readRDS(file.path(nn_dir, nn_files$ligand_target))
weighted_nets <- readRDS(file.path(nn_dir, nn_files$weighted_nets))
# ligand_tf loaded later only if signaling path step is run

obj <- readRDS(obj_path(cfg, "05_annotated"))
Idents(obj) <- "celltype"
cat("Cell types:\n"); print(sort(table(obj$celltype)))

# ---- Step 2: Receiver and sender cell types -----------------------------
# Receiver: T cells (biologically the primary NicheNet target in PBMCs)
all_ct   <- unique(obj$celltype)
receiver <- grep("T_cell|T cell|CD4|CD8", all_ct, value = TRUE, ignore.case = TRUE)[1]
if (is.na(receiver)) receiver <- all_ct[1]
senders  <- setdiff(all_ct, receiver)
cat(sprintf("\nReceiver: %s\nSenders : %s\n", receiver, paste(senders, collapse = ", ")))

# ---- Step 3: Expressed genes --------------------------------------------
# pct = 0.10: a gene is "expressed" if detected in ≥10% of cells
expressed_genes_receiver <- get_expressed_genes(
  receiver, obj, celltype_colname = "celltype", pct = 0.10)
expressed_genes_sender   <- unique(unlist(lapply(senders, get_expressed_genes,
  seurat_obj = obj, celltype_colname = "celltype", pct = 0.10)))

background_expressed_genes <- expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target)]
cat(sprintf("Background genes: %d\n", length(background_expressed_genes)))

# ---- Step 4: Gene set of interest (receiver DEGs, Post vs Healthy) ------
primary_grp <- cfg$de_contrasts[[1]]$group    # e.g. "Post"
primary_ref <- cfg$de_contrasts[[1]]$ref      # e.g. "Healthy"
cat(sprintf("\nDE contrast for gene set: %s vs %s\n", primary_grp, primary_ref))

sub_receiver <- subset(obj, idents = receiver)
Idents(sub_receiver) <- "condition"
de_receiver  <- FindMarkers(sub_receiver,
                            ident.1 = primary_grp,
                            ident.2 = primary_ref,
                            min.pct  = 0.10,
                            verbose  = FALSE)
geneset_oi <- de_receiver %>%
  filter(p_val_adj <= 0.05, avg_log2FC >= 0.25) %>%
  rownames() %>%
  intersect(rownames(ligand_target))
cat(sprintf("Gene set of interest: %d genes\n", length(geneset_oi)))
if (length(geneset_oi) < 10)
  warning("Gene set has <10 genes — NicheNet results may be unreliable. ",
          "Check the contrast or lower the FC threshold.")

# ---- Step 5: Potential ligands ------------------------------------------
expressed_receptors <- filter(lr_network, to %in% expressed_genes_receiver)$from %>%
                         intersect(expressed_genes_receiver)   # receptors that are expressed
# Agnostic: any ligand from lr_network with an expressed receptor
potential_ligands_agnostic <- filter(lr_network, to %in% expressed_genes_receiver)$from %>%
  unique() %>%
  intersect(rownames(obj))
# Focused: only sender-expressed ligands with expressed receptors
potential_ligands_focused <- potential_ligands_agnostic %>%
  intersect(expressed_genes_sender)

cat(sprintf("Potential ligands — agnostic: %d  focused: %d\n",
            length(potential_ligands_agnostic), length(potential_ligands_focused)))

# ---- Step 6: Ligand activity prediction ---------------------------------
message("Predicting ligand activities (agnostic)...")
ligand_activities_agnostic <- predict_ligand_activities(
  geneset            = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target,
  potential_ligands  = potential_ligands_agnostic)

message("Predicting ligand activities (focused)...")
ligand_activities_focused <- predict_ligand_activities(
  geneset            = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target,
  potential_ligands  = potential_ligands_focused)

# ---- Step 7: Top 20 ligands per approach, take union --------------------
top_n <- 20L
top_agnostic <- ligand_activities_agnostic %>%
  arrange(desc(aupr_corrected)) %>%
  slice_head(n = top_n) %>%
  pull(test_ligand)
top_focused  <- ligand_activities_focused %>%
  arrange(desc(aupr_corrected)) %>%
  slice_head(n = top_n) %>%
  pull(test_ligand)
best_ligands <- union(top_agnostic, top_focused)
cat(sprintf("\nTop ligands — agnostic: %d  focused: %d  union: %d\n",
            length(top_agnostic), length(top_focused), length(best_ligands)))
cat("Best upstream ligands:", paste(head(best_ligands, 15), collapse = ", "), "\n")

# Save ligand activities
write.csv(ligand_activities_focused %>% arrange(desc(aupr_corrected)),
          out_path(cfg, "tables", "spoke_nichenet_ligands.csv"), row.names = FALSE)

# ---- Plots --------------------------------------------------------------
plot_dir <- file.path(cfg$project$outdir, "plots", "nichenet")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Helper to build ligand activity barplot (AUPR)
ligand_activity_barplot <- function(activities, top_ligands, title) {
  df <- activities %>%
    filter(test_ligand %in% top_ligands) %>%
    mutate(test_ligand = factor(test_ligand,
                                levels = rev(top_ligands[top_ligands %in% test_ligand])))
  ggplot(df, aes(x = test_ligand, y = aupr_corrected)) +
    geom_col(fill = "#2196F3") +
    coord_flip() +
    labs(x = NULL, y = "AUPR (corrected)", title = title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))
}

# ---- Step 8: Activity barplots ------------------------------------------
p_agn <- ligand_activity_barplot(ligand_activities_agnostic, top_agnostic, "Agnostic")
p_foc <- ligand_activity_barplot(ligand_activities_focused,  top_focused,  "Focused")
ggsave(file.path(plot_dir, "ligand_activity_barplot.pdf"),
       p_agn + p_foc, width = 12, height = 6)
ggsave(file.path(plot_dir, "ligand_activity_barplot.png"),
       p_agn + p_foc, width = 12, height = 6)
rm(p_agn, p_foc); invisible(gc())

# ---- Step 9: Ligand-target heatmap --------------------------------------
tryCatch({
  lt_links <- lapply(best_ligands, function(lig)
    get_weighted_ligand_target_links(
      lig, geneset_oi, ligand_target, n = 250)) %>%
    bind_rows() %>%
    drop_na()
  lt_vis <- prepare_ligand_target_visualization(lt_links,
    ligand_activities_focused, cutoff = 1/3)
  order_ligands <- intersect(rev(best_ligands), rownames(lt_vis))
  order_targets <- colnames(lt_vis)[order(colSums(lt_vis), decreasing = TRUE)]
  p_lt <- make_heatmap_ggplot(
    t(lt_vis[order_ligands, order_targets]),
    y_name = "target", x_name = "ligand",
    color = "purple", legend_title = "Regulatory\npotential")
  ggsave(file.path(plot_dir, "ligand_target_heatmap.pdf"), p_lt, width = 10, height = 8)
  ggsave(file.path(plot_dir, "ligand_target_heatmap.png"), p_lt, width = 10, height = 8)
  rm(lt_links, lt_vis, p_lt); invisible(gc())
}, error = function(e)
  message("[skip] ligand-target heatmap: ", conditionMessage(e)))

# ---- Step 10: Ligand-receptor heatmap -----------------------------------
tryCatch({
  lr_links <- get_weighted_ligand_receptor_links(
    best_ligands, expressed_receptors,
    lr_network, weighted_nets$lr_sig)
  lr_vis <- prepare_ligand_receptor_visualization(lr_links,
    best_ligands, order_hclust = "both")
  p_lr <- make_heatmap_ggplot(
    t(lr_vis),
    y_name = "receptor", x_name = "ligand",
    color = "mediumvioletred", legend_title = "Prior\ninteraction")
  ggsave(file.path(plot_dir, "ligand_receptor_heatmap.pdf"), p_lr, width = 10, height = 7)
  ggsave(file.path(plot_dir, "ligand_receptor_heatmap.png"), p_lr, width = 10, height = 7)
  rm(lr_links, lr_vis, p_lr); invisible(gc())
}, error = function(e)
  message("[skip] ligand-receptor heatmap: ", conditionMessage(e)))

# ---- Step 11: Sender expression heatmap ---------------------------------
tryCatch({
  sub_senders <- subset(obj, idents = c(senders, receiver))
  Idents(sub_senders) <- "celltype"
  # Seurat 5: layer= not slot=
  avg_expr <- tryCatch(
    AverageExpression(sub_senders, features = best_ligands, assay = "RNA",
                      layer = "data")$RNA,
    error = function(e)
      AverageExpression(sub_senders, features = best_ligands, assay = "RNA",
                        slot  = "data")$RNA)
  avg_expr <- avg_expr[intersect(best_ligands, rownames(avg_expr)), , drop = FALSE]
  avg_scaled <- t(scale(t(avg_expr)))

  # Try pheatmap; fall back to ComplexHeatmap
  sender_hm <- tryCatch({
    suppressPackageStartupMessages(library(pheatmap))
    png(file.path(plot_dir, "ligand_sender_heatmap.png"), width = 1000, height = 800)
    pheatmap(avg_scaled,
             color           = colorRampPalette(c("navy","white","firebrick3"))(50),
             cluster_rows    = TRUE,
             cluster_cols    = TRUE,
             border_color    = NA,
             fontsize_row    = 8,
             fontsize_col    = 9)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    # ComplexHeatmap fallback
    tryCatch({
      col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
      png(file.path(plot_dir, "ligand_sender_heatmap.png"), width = 1000, height = 800)
      draw(Heatmap(avg_scaled, col = col_fun, name = "Scaled\nexpr",
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 9)))
      dev.off()
    }, error = function(e2) {
      if (dev.cur() > 1) dev.off()
      message("[skip] sender heatmap: ", conditionMessage(e2))
    })
  })
  rm(sub_senders, avg_expr, avg_scaled); invisible(gc())
}, error = function(e)
  message("[skip] sender expression heatmap: ", conditionMessage(e)))

# ---- Step 12: Extended prioritization (optional) -------------------------
tryCatch({
  message("\nRunning extended prioritization...")
  lr_filtered <- filter(lr_network,
                        from %in% best_ligands,
                        to   %in% expressed_genes_receiver)
  lr_features <- union(lr_filtered$from, lr_filtered$to) %>%
                   intersect(rownames(obj))

  sub_oi <- subset(obj, condition == primary_grp)
  Idents(sub_oi) <- "celltype"
  fm_all <- FindAllMarkers(sub_oi, features = lr_features, min.pct = 0,
                           logfc.threshold = 0, return.thresh = 1,
                           verbose = FALSE)

  avg_oi <- get_exprs_avg(sub_oi, celltype_colname = "celltype",
                          condition_colname = "condition",
                          condition_oi = primary_grp,
                          features = lr_features)
  de_cond <- FindMarkers(sub_receiver, ident.1 = primary_grp, ident.2 = primary_ref,
                         features = lr_features, min.pct = 0,
                         logfc.threshold = 0, verbose = FALSE) %>%
             rownames_to_column("gene")

  ic_sender   <- process_table_to_ic(fm_all,  table_type = "celltype_DE",
                                     lr_network = lr_filtered,
                                     senders_oi = senders, receivers_oi = receiver)
  ic_expr     <- process_table_to_ic(avg_oi,  table_type = "expression",
                                     lr_network = lr_filtered)
  ic_group    <- process_table_to_ic(de_cond, table_type = "group_DE",
                                     lr_network = lr_filtered)
  prio_table  <- generate_prioritization_tables(
    sender_receiver_info            = ic_sender,
    ligand_activities               = ligand_activities_focused,
    sender_receiver_de              = ic_group,
    ligand_expression_de            = ic_expr,
    scenario                        = "case_control")
  write.csv(prio_table,
            out_path(cfg, "tables", "spoke_nichenet_prioritized_lr.csv"),
            row.names = FALSE)
  cat("\nTop 10 prioritized L-R pairs:\n")
  print(head(prio_table[, c("sender","receiver","ligand","receptor","prioritization_score")], 10))
  rm(sub_oi, fm_all, avg_oi, de_cond, ic_sender, ic_expr, ic_group); invisible(gc())
}, error = function(e)
  message("[skip] extended prioritization: ", conditionMessage(e)))

# ---- Step 13: Signaling path (top ligand → top target) ------------------
tryCatch({
  suppressPackageStartupMessages({ library(igraph); library(ggraph) })
  lp <- file.path(nn_dir, nn_files$ligand_tf)
  if (file.exists(lp)) {
    ligand_tf <- readRDS(lp)
    top_lig   <- best_ligands[1]
    top_tgt   <- names(which.max(ligand_target[, top_lig, drop = FALSE][geneset_oi, ]))
    if (length(top_tgt) == 0)
      top_tgt <- geneset_oi[which.max(ligand_target[geneset_oi, top_lig])]
    sp_obj <- get_ligand_signaling_path(
      ligand_tf_matrix   = ligand_tf,
      ligands_all        = top_lig,
      targets_all        = top_tgt,
      top_n_regulators   = 4,
      weighted_networks  = weighted_nets,
      minmax_scaling     = TRUE)
    edges <- bind_rows(sp_obj$sig %>% mutate(type = "signaling"),
                       sp_obj$gr  %>% mutate(type = "GRN"))
    g <- graph_from_data_frame(edges, directed = TRUE)
    p_sp <- ggraph(g, layout = "kk") +
      geom_edge_link(aes(width = weight, alpha = weight),
                     arrow = arrow(length = unit(3, "mm"), type = "closed"),
                     end_cap = circle(3, "mm"), color = "#555555") +
      scale_edge_width(range = c(0.5, 2)) +
      geom_node_point(size = 5, color = "#E41A1C") +
      geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
      theme_void() +
      labs(title = paste0(top_lig, " → ", top_tgt, " signaling path"))
    ggsave(file.path(plot_dir, "signaling_path.pdf"), p_sp, width = 8, height = 6)
    ggsave(file.path(plot_dir, "signaling_path.png"), p_sp, width = 8, height = 6)
    rm(ligand_tf, sp_obj, edges, g, p_sp); invisible(gc())
  } else {
    message("[skip] signaling path: ligand_tf matrix not downloaded yet")
  }
}, error = function(e)
  message("[skip] signaling path: ", conditionMessage(e)))

# ---- Step 14: Sender assignment -----------------------------------------
tryCatch({
  sub_senders2 <- subset(obj, idents = senders)
  Idents(sub_senders2) <- "celltype"
  avg_snd <- tryCatch(
    AverageExpression(sub_senders2, features = best_ligands, assay = "RNA",
                      layer = "data")$RNA,
    error = function(e)
      AverageExpression(sub_senders2, features = best_ligands, assay = "RNA",
                        slot  = "data")$RNA)
  avg_snd <- avg_snd[intersect(best_ligands, rownames(avg_snd)), , drop = FALSE]
  best_sender <- apply(avg_snd, 1, function(x) colnames(avg_snd)[which.max(x)])
  sender_summary <- tibble(
    ligand           = names(best_sender),
    best_sender      = unname(best_sender)) %>%
    left_join(
      ligand_activities_focused %>%
        select(test_ligand, aupr_corrected, auroc, pearson) %>%
        rename(ligand = test_ligand),
      by = "ligand") %>%
    arrange(desc(aupr_corrected))
  cat("\nSender assignment summary:\n"); print(sender_summary)
  write.csv(sender_summary,
            out_path(cfg, "tables", "spoke_nichenet_sender_summary.csv"),
            row.names = FALSE)
  rm(sub_senders2, avg_snd); invisible(gc())
}, error = function(e)
  message("[skip] sender assignment: ", conditionMessage(e)))

# ---- Step 15: Circos plot -----------------------------------------------
tryCatch({
  # For each top ligand, take its 5 highest-weight targets in geneset_oi
  circos_links <- lapply(head(best_ligands, 15), function(lig) {
    wts <- ligand_target[geneset_oi, lig, drop = FALSE]
    top5 <- names(sort(wts[, 1], decreasing = TRUE)[seq_len(min(5, nrow(wts)))])
    tibble(sender = lig, target = top5, weight = wts[top5, 1])
  }) %>% bind_rows()

  # Color senders (ligands) distinctly; targets in grey
  n_lig    <- length(unique(circos_links$sender))
  lig_cols <- setNames(colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                                          "#FF7F00","#A65628"))(n_lig),
                       unique(circos_links$sender))
  all_nodes <- union(circos_links$sender, circos_links$target)
  node_cols <- setNames(
    ifelse(all_nodes %in% names(lig_cols), lig_cols[all_nodes], "grey70"),
    all_nodes)

  pdf(file.path(plot_dir, "nichenet_circos.pdf"), width = 10, height = 10)
  circos.clear()
  chordDiagram(circos_links[, c("sender","target","weight")],
               grid.col      = node_cols,
               transparency  = 0.4,
               annotationTrack = "grid",
               preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim   <- get.cell.meta.data("xlim")
    xplot  <- get.cell.meta.data("xplot")
    ylim   <- get.cell.meta.data("ylim")
    sector <- get.cell.meta.data("sector.index")
    if (abs(xplot[2] - xplot[1]) > 5)
      circos.text(mean(xlim), ylim[1], sector, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
  }, bg.border = NA)
  title(paste("NicheNet: sender → target ligand network\n",
              primary_grp, "vs", primary_ref, "—", receiver))
  circos.clear()
  dev.off()
  message("Circos plot saved.")
}, error = function(e) {
  if (dev.cur() > 1) { circos.clear(); dev.off() }
  message("[skip] circos: ", conditionMessage(e))
})

log_session(cfg, "spoke_nichenet")
