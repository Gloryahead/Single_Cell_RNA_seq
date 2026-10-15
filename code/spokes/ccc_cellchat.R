#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# ccc_cellchat.R — cell-cell communication via CellChat.          (Part 9)
#
# NOTE: CellChat is NOT in the current scrna_r.sif. Rebuild with:
#   devtools::install_github("jinworks/CellChat")
#   BiocManager::install(c("BiocNeighbors", "ComplexHeatmap"))
#   install.packages(c("NMF", "circlize", "ggalluvial"))
# before enabling this rule in Snakefile rule all.
#
# Strategy: build one CellChat object per condition, then compare
# the primary contrast (Healthy vs Post) with mergeCellChat().
# Per-condition objects are saved to objects/cellchat_<condition>.rds.
# The merged Healthy+Post object is the declared Snakemake output.
#
# In : results/objects/05_annotated.rds
# Out: results/objects/spoke_cellchat.rds   (merged Healthy+Post)
#      results/objects/cellchat_<condition>.rds  (per condition)
#      results/plots/cellchat/<condition>/  (per-condition plots)
#      results/plots/cellchat/comparison/   (differential plots)
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(CellChat); library(patchwork); library(NMF); library(ggalluvial)
})
options(stringsAsFactors = FALSE)
registerDoSEQ()   # suppress %dopar% warnings when running sequentially

cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))
sp  <- if (cfg$project$organism == "mouse") "mouse" else "human"
db  <- if (sp == "mouse") CellChatDB.mouse else CellChatDB.human

# ---- helpers ------------------------------------------------------------

# Safety wrapper: opens a PNG device, runs plot_fn(), closes it.
# Ensures dev.off() is always called even if the plot fails.
png_safe <- function(path, plot_fn, w = 800, h = 800, res = 100) {
  tryCatch({
    png(path, width = w, height = h, res = res)
    plot_fn()
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("  [skip] ", basename(path), ": ", conditionMessage(e))
  })
}

gg_safe <- function(p, path_base, w = 10, h = 8) {
  tryCatch({
    ggsave(paste0(path_base, ".png"), p, width = w, height = h)
    ggsave(paste0(path_base, ".pdf"), p, width = w, height = h)
  }, error = function(e)
    message("  [skip] ", basename(path_base), ": ", conditionMessage(e)))
}

# ---- per-condition CellChat pipeline ------------------------------------
build_cellchat <- function(sub_obj, cond_name) {
  message("\n=== Building CellChat: ", cond_name, " (", ncol(sub_obj), " cells) ===")
  samples_vec <- factor(sub_obj$sample)

  cc <- createCellChat(object   = sub_obj,
                       group.by = "celltype",
                       assay    = "RNA",
                       samples  = samples_vec)
  cc@DB <- db
  cc <- subsetData(cc)
  # do.fast = FALSE: exact statistical test (slower but more accurate)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", nboot = 100)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  cat(sprintf("  %s: %d L-R pairs, %d signaling pathways\n",
              cond_name,
              nrow(subsetCommunication(cc)),
              length(cc@netP$pathways)))
  cc
}

# ---- per-condition visualisation ----------------------------------------
plot_condition <- function(cc, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ct_levels <- levels(cc@idents)

  # Dynamic T-cell receiver indices for hierarchy plots
  t_idx <- grep("T_cell|T cell|CD4|CD8", ct_levels, ignore.case = TRUE)
  if (length(t_idx) == 0) t_idx <- seq_len(min(3L, length(ct_levels)))

  # Top signaling pathway (by total probability mass)
  path_totals <- sapply(cc@netP$pathways, function(p)
    sum(cc@netP$prob[, , p], na.rm = TRUE))
  top_path <- if (length(path_totals) > 0) names(which.max(path_totals)) else NULL

  # 01 Circle — interaction count
  png_safe(file.path(out_dir, "01_circle_count.png"), function()
    netVisual_circle(cc@net$count, weight.scale = TRUE, label.edge = FALSE,
                     title.name = "Number of interactions"), w = 800, h = 800)

  # 02 Circle — interaction strength
  png_safe(file.path(out_dir, "02_circle_strength.png"), function()
    netVisual_circle(cc@net$weight, weight.scale = TRUE, label.edge = FALSE,
                     title.name = "Interaction strength"), w = 800, h = 800)

  # 03-04 Top pathway: hierarchy + circle
  if (!is.null(top_path)) {
    png_safe(file.path(out_dir, "03_top_pathway_hierarchy.png"), function()
      netVisual_aggregate(cc, signaling = top_path, layout = "hierarchy",
                          vertex.receiver = t_idx), w = 1000, h = 800)
    png_safe(file.path(out_dir, "04_top_pathway_circle.png"), function()
      netVisual_aggregate(cc, signaling = top_path, layout = "circle"), w = 800, h = 800)
    # L-R contribution to top pathway
    tryCatch({
      enriched_lr <- extractEnrichedLR(cc, signaling = top_path, geneLR.return = FALSE)
      p_contrib <- netAnalysis_contribution(cc, signaling.name = top_path)
      gg_safe(p_contrib, file.path(out_dir, "05_top_pathway_LR_contribution"), w = 8, h = 5)
    }, error = function(e)
      message("  [skip] LR contribution: ", conditionMessage(e)))
  }

  # 06 Bubble — all L-R pairs
  tryCatch({
    p_bub <- netVisual_bubble(cc, remove.isolate = FALSE, angle.x = 45)
    gg_safe(p_bub, file.path(out_dir, "06_bubble_all_LR"), w = 14, h = 10)
  }, error = function(e) message("  [skip] bubble: ", conditionMessage(e)))

  # 07 Signaling role heatmap — outgoing
  png_safe(file.path(out_dir, "07_role_heatmap_outgoing.png"), function()
    netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing", height = 12),
    w = 1000, h = 1200)

  # 08 Signaling role heatmap — incoming
  png_safe(file.path(out_dir, "08_role_heatmap_incoming.png"), function()
    netAnalysis_signalingRole_heatmap(cc, pattern = "incoming", height = 12),
    w = 1000, h = 1200)

  # 09 Signaling role scatter
  tryCatch({
    p_scat <- netAnalysis_signalingRole_scatter(cc)
    gg_safe(p_scat, file.path(out_dir, "09_role_scatter"), w = 8, h = 6)
  }, error = function(e) message("  [skip] role scatter: ", conditionMessage(e)))

  # 10-12 Communication patterns — outgoing (selectK + river + dot)
  tryCatch({
    png_safe(file.path(out_dir, "10_selectK_outgoing.png"), function()
      selectK(cc, pattern = "outgoing"), w = 900, h = 600)
    cc_pat <- identifyCommunicationPatterns(cc, pattern = "outgoing",
                                            k = 4L, verbose = FALSE)
    tryCatch({
      p_riv <- netAnalysis_river(cc_pat, pattern = "outgoing")
      gg_safe(p_riv, file.path(out_dir, "11_pattern_river_outgoing"), w = 10, h = 8)
    }, error = function(e) message("  [skip] outgoing river: ", conditionMessage(e)))
    tryCatch({
      p_dot <- netAnalysis_dot(cc_pat, pattern = "outgoing")
      gg_safe(p_dot, file.path(out_dir, "12_pattern_dot_outgoing"), w = 9, h = 6)
    }, error = function(e) message("  [skip] outgoing dot: ", conditionMessage(e)))
  }, error = function(e)
    message("  [skip] outgoing patterns: ", conditionMessage(e)))

  # 13-15 Communication patterns — incoming
  tryCatch({
    cc_pat2 <- identifyCommunicationPatterns(cc, pattern = "incoming",
                                             k = 4L, verbose = FALSE)
    tryCatch({
      p_riv2 <- netAnalysis_river(cc_pat2, pattern = "incoming")
      gg_safe(p_riv2, file.path(out_dir, "13_pattern_river_incoming"), w = 10, h = 8)
    }, error = function(e) message("  [skip] incoming river: ", conditionMessage(e)))
    tryCatch({
      p_dot2 <- netAnalysis_dot(cc_pat2, pattern = "incoming")
      gg_safe(p_dot2, file.path(out_dir, "14_pattern_dot_incoming"), w = 9, h = 6)
    }, error = function(e) message("  [skip] incoming dot: ", conditionMessage(e)))
  }, error = function(e)
    message("  [skip] incoming patterns: ", conditionMessage(e)))

  # 15-16 Pathway functional similarity embedding (needs umap-learn via reticulate)
  tryCatch({
    cc <- computeNetSimilarity(cc, type = "functional")
    cc <- netEmbedding(cc, type = "functional")
    cc <- netClustering(cc, type = "functional")
    png_safe(file.path(out_dir, "15_pathway_embedding.png"), function()
      netVisual_embedding(cc, type = "functional", label.size = 3.5), w = 1000, h = 800)
    png_safe(file.path(out_dir, "16_pathway_embedding_zoom.png"), function()
      netVisual_embeddingZoomIn(cc, type = "functional", nCol = 2), w = 1200, h = 1000)
  }, error = function(e)
    message("  [skip] pathway embedding (umap-learn?): ", conditionMessage(e)))

  invisible(cc)
}

# ---- comparison plots (merged CellChat object) -------------------------
plot_comparison <- function(cc_merged, cond_names, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  object.list <- cc_merged@meta.list
  n           <- length(object.list)
  grp         <- seq_len(n)
  ref_name    <- cond_names[1]   # Healthy (reference)
  cmp_name    <- cond_names[2]   # Post

  # 01 Interaction count + strength bar chart
  tryCatch({
    p_cnt <- compareInteractions(cc_merged, show.legend = FALSE, group = grp,
                                 measure = "count")
    p_str <- compareInteractions(cc_merged, show.legend = FALSE, group = grp,
                                 measure = "weight")
    gg_safe(p_cnt + p_str, file.path(out_dir, "01_interaction_count_strength"), w = 10, h = 5)
  }, error = function(e) message("  [skip] compare bar: ", conditionMessage(e)))

  # 02-03 Circle per condition (loop) + diff interaction circles
  for (i in seq_len(n)) {
    png_safe(file.path(out_dir, sprintf("02_circle_%s.png", cond_names[i])),
             function() netVisual_circle(object.list[[i]]@net$count,
                                         weight.scale = TRUE, label.edge = FALSE,
                                         title.name = cond_names[i]), w = 800, h = 800)
  }
  png_safe(file.path(out_dir, "03a_diff_count.png"), function()
    netVisual_diffInteraction(cc_merged, weight.scale = TRUE), w = 800, h = 800)
  png_safe(file.path(out_dir, "03b_diff_strength.png"), function()
    netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight"),
    w = 800, h = 800)

  # 04 Differential heatmaps
  tryCatch({
    png_safe(file.path(out_dir, "04a_diff_heatmap_count.png"), function()
      netVisual_heatmap(cc_merged), w = 900, h = 900)
    png_safe(file.path(out_dir, "04b_diff_heatmap_strength.png"), function()
      netVisual_heatmap(cc_merged, measure = "weight"), w = 900, h = 900)
  }, error = function(e) message("  [skip] diff heatmap: ", conditionMessage(e)))

  # 05-06 rankNet: pathway-level strength comparison
  tryCatch({
    p_stk <- rankNet(cc_merged, mode = "comparison", stacked = TRUE,  do.stat = TRUE)
    p_sbs <- rankNet(cc_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
    gg_safe(p_stk, file.path(out_dir, "05_rankNet_stacked"),   w = 10, h = 8)
    gg_safe(p_sbs, file.path(out_dir, "06_rankNet_sidebyside"), w = 10, h = 8)
  }, error = function(e) message("  [skip] rankNet: ", conditionMessage(e)))

  # 07-10 Signaling role heatmaps per condition, outgoing + incoming
  for (i in seq_len(n)) {
    for (pat in c("outgoing", "incoming")) {
      png_safe(file.path(out_dir,
                         sprintf("%02d_role_%s_%s.png", 6 + 2*(i-1) + (pat=="incoming"), pat, cond_names[i])),
               function() netAnalysis_signalingRole_heatmap(
                 object.list[[i]], pattern = pat, title = paste(cond_names[i], pat),
                 height = 12, color.heatmap = if (pat == "outgoing") "GnBu" else "OrRd"),
               w = 1000, h = 1200)
    }
  }

  # 11 Bubble: monocyte sender → T/NK receivers, both conditions
  tryCatch({
    all_ct  <- levels(object.list[[1]]@idents)
    mono_src <- grep("Mono|monocyte|Macrophage|DC", all_ct, value = TRUE, ignore.case = TRUE)
    tnk_tgt  <- grep("T_cell|T cell|NK", all_ct, value = TRUE, ignore.case = TRUE)
    if (length(mono_src) > 0 && length(tnk_tgt) > 0) {
      p_bub2 <- netVisual_bubble(cc_merged,
                                  sources.use  = mono_src[1],
                                  targets.use  = tnk_tgt,
                                  comparison   = grp,
                                  angle.x = 45,
                                  remove.isolate = TRUE)
      gg_safe(p_bub2, file.path(out_dir, "11_bubble_mono_to_TNK"), w = 12, h = 8)
    }
  }, error = function(e) message("  [skip] mono→TNK bubble: ", conditionMessage(e)))

  # 12-13 Chord gene plot: monocyte→T/NK per condition
  tryCatch({
    all_ct  <- levels(object.list[[1]]@idents)
    mono_src <- grep("Mono|monocyte|Macrophage|DC", all_ct, value = TRUE, ignore.case = TRUE)
    tnk_tgt  <- grep("T_cell|T cell|NK", all_ct, value = TRUE, ignore.case = TRUE)
    if (length(mono_src) > 0 && length(tnk_tgt) > 0) {
      for (i in seq_len(n)) {
        png_safe(file.path(out_dir, sprintf("12_chord_mono_TNK_%s.png", cond_names[i])),
                 function() netVisual_chord_gene(object.list[[i]],
                                                  sources.use = mono_src[1],
                                                  targets.use = tnk_tgt,
                                                  lab.cex = 0.5,
                                                  title.name = cond_names[i]),
                 w = 900, h = 900)
      }
    }
  }, error = function(e) message("  [skip] chord: ", conditionMessage(e)))

  # 13-15 Pairwise functional embedding comparison
  tryCatch({
    cc_merged <- computeNetSimilarityPairwise(cc_merged, type = "functional")
    cc_merged <- netEmbedding(cc_merged, type = "functional")
    cc_merged <- netClustering(cc_merged, type = "functional")
    png_safe(file.path(out_dir, "13_pairwise_embedding.png"), function()
      netVisual_embeddingPairwise(cc_merged, type = "functional",
                                  label.size = 3.5), w = 1200, h = 900)
    png_safe(file.path(out_dir, "14_pairwise_embedding_zoom.png"), function()
      netVisual_embeddingPairwiseZoomIn(cc_merged, type = "functional",
                                        nCol = 2), w = 1200, h = 1000)
    tryCatch({
      p_sim <- rankSimilarity(cc_merged, type = "functional")
      gg_safe(p_sim, file.path(out_dir, "15_rankSimilarity"), w = 8, h = 6)
    }, error = function(e) message("  [skip] rankSimilarity: ", conditionMessage(e)))
  }, error = function(e)
    message("  [skip] pairwise embedding: ", conditionMessage(e)))

  invisible(cc_merged)
}

# =========================================================================
# Main
# =========================================================================
conditions <- sort(unique(obj$condition))
cat("Conditions found:", paste(conditions, collapse = ", "), "\n")

cc_by_cond <- list()
for (cond in conditions) {
  sub      <- obj[, obj$condition == cond]
  cc_cond  <- build_cellchat(sub, cond)
  out_dir  <- file.path(cfg$project$outdir, "plots", "cellchat", cond)
  plot_condition(cc_cond, out_dir)
  saveRDS(cc_cond, out_path(cfg, "objects", paste0("cellchat_", cond, ".rds")))
  cc_by_cond[[cond]] <- cc_cond
  rm(sub, cc_cond); invisible(gc())
}

# ---- Primary comparison: Healthy vs Post --------------------------------
# Identify reference and comparison from the first de_contrast in config.
primary_ref <- cfg$de_contrasts[[1]]$ref    # e.g. "Healthy"
primary_grp <- cfg$de_contrasts[[1]]$group  # e.g. "Post"

pair <- c(primary_ref, primary_grp)
pair <- pair[pair %in% names(cc_by_cond)]   # keep only available conditions

cc_merged <- NULL
if (length(pair) == 2) {
  message("\n=== Merging CellChat: ", paste(pair, collapse = " vs "), " ===")
  cc_merged <- mergeCellChat(
    list(cc_by_cond[[pair[1]]], cc_by_cond[[pair[2]]]),
    add.names = pair,
    cell.level = FALSE
  )
  cmp_dir <- file.path(cfg$project$outdir, "plots", "cellchat", "comparison")
  plot_comparison(cc_merged, pair, cmp_dir)
} else {
  message("WARNING: primary contrast conditions not both available for comparison. ",
          "Found: ", paste(names(cc_by_cond), collapse = ", "))
  # Fall back: save the first available condition object as the output
  cc_merged <- cc_by_cond[[1]]
}

saveRDS(cc_merged, out_path(cfg, "objects", "spoke_cellchat.rds"))
log_session(cfg, "spoke_cellchat")
