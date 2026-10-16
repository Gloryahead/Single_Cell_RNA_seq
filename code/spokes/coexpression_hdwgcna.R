#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# coexpression_hdwgcna.R — gene co-expression networks (hdWGCNA). (Part 12)
#
# NOTE: hdWGCNA is NOT in the current scrna_r.sif. Rebuild with:
#   devtools::install_github("smorabit/hdWGCNA", ref = "dev")
#   BiocManager::install(c("WGCNA", "UCell", "GeneOverlap"))
#   devtools::install_github("wjawaid/enrichR")   # enrichR GitHub version
#   install.packages(c("tidyverse","cowplot","igraph"))
#
# Memory strategy: subset to the target cell type BEFORE the pipeline.
# ModuleEigengenes densifies a cells × genes matrix; keeping all non-target
# cells resident would OOM at 90 GB. Subsetting drops ~35 GB of unused cells.
# MetacellsByGroups only uses target_ct metacells anyway, so results are identical.
#
# In : results/objects/05_annotated.rds
# Out: results/objects/spoke_hdwgcna.rds
#      results/tables/spoke_hdwgcna_*.csv
#      results/plots/spoke_hdwgcna_*.{pdf,png}
#      results/plots/hdwgcna_enrichr/
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(hdWGCNA); library(WGCNA); library(Seurat)
  library(tidyverse); library(patchwork); library(cowplot)
})
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
enableWGCNAThreads(nThreads = n_threads)
options(future.globals.maxSize = 8 * 1024^3)   # 8 GB for HubGeneNetworkPlot
cfg <- load_config(); set.seed(cfg$seed)

obj <- readRDS(obj_path(cfg, "05_annotated"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)

# ---- Step 1: Target cell type (most abundant → best metacell coverage) --
ct_counts <- sort(table(obj$celltype), decreasing = TRUE)
target_ct <- names(ct_counts)[1]
message("hdWGCNA target cell type: ", target_ct, " (", ct_counts[1], " cells)")

# Subset before pipeline — see memory note in header
obj <- subset(obj, subset = celltype == target_ct)
cat(sprintf("After subset: %d cells\n", ncol(obj)))

wgcna_name <- paste0(gsub("[^A-Za-z0-9]", "_", target_ct), "_WGCNA")

# ---- Step 2: Setup & metacells ------------------------------------------
obj <- SetupForWGCNA(obj,
                     gene_select = "fraction",
                     fraction    = 0.05,
                     wgcna_name  = wgcna_name)
cat(sprintf("WGCNA genes selected: %d\n", length(GetWGCNAGenes(obj))))

# group.by uses sample (orig.ident equivalent) not condition — gives better
# within-group metacell resolution for a 4-sample-per-condition design.
# max_shared = 10 limits cell reuse across metacells.
obj <- MetacellsByGroups(obj,
                         group.by   = c("celltype", "sample"),
                         k          = 20,
                         max_shared = 10,
                         reduction  = "harmony",
                         ident.group = "celltype")
obj <- NormalizeMetacells(obj)

# ---- Step 3: Expression data for WGCNA ----------------------------------
obj <- SetDatExpr(obj,
                  group_name = target_ct,
                  group.by   = "celltype",
                  assay      = "RNA",
                  layer      = "data")

# ---- Step 4: Soft-thresholding power ------------------------------------
obj <- TestSoftPowers(obj, networkType = "signed")

# Auto-select: lowest power where signed R² ≥ 0.85; fall back to 6
power_tbl  <- GetPowerTable(obj)
best_power <- power_tbl %>%
  filter(SFT.R.sq >= 0.85) %>%
  arrange(Power) %>%
  slice(1) %>%
  pull(Power)
if (length(best_power) == 0 || is.na(best_power)) {
  best_power <- 6L
  message("R² ≥ 0.85 not reached — using fallback soft_power = 6")
}
message("Selected soft power: ", best_power)

plot_list_sp <- PlotSoftPowers(obj)
pdf(out_path(cfg, "plots", "spoke_hdwgcna_soft_power.pdf"), width = 10, height = 5)
print(wrap_plots(plot_list_sp, ncol = 2))
dev.off()
ggsave(out_path(cfg, "plots", "spoke_hdwgcna_soft_power.png"),
       wrap_plots(plot_list_sp, ncol = 2), width = 10, height = 5)
rm(plot_list_sp); invisible(gc())

# ---- Step 5: Network construction ---------------------------------------
# overwrite_tom = TRUE: re-runs don't block on a stale TOM file from a prior
# failed attempt (ConstructNetwork writes TOM to disk before ModuleEigengenes).
obj <- ConstructNetwork(obj,
                        soft_power     = best_power,
                        setDatExpr     = FALSE,
                        networkType    = "signed",
                        TOMType        = "signed",
                        minModuleSize  = 30,
                        mergeCutHeight = 0.25,
                        overwrite_tom  = TRUE)

pdf(out_path(cfg, "plots", "spoke_hdwgcna_dendrogram.pdf"), width = 10, height = 5)
PlotDendrogram(obj, main = paste0(target_ct, " hdWGCNA dendrogram"))
dev.off()
png(out_path(cfg, "plots", "spoke_hdwgcna_dendrogram.png"), width = 1000, height = 500)
PlotDendrogram(obj, main = paste0(target_ct, " hdWGCNA dendrogram"))
dev.off()

# ---- Step 6: Module eigengenes ------------------------------------------
# ScaleData on WGCNA genes before ModuleEigengenes (needed for PCA within fn)
wgcna_genes <- GetWGCNAGenes(obj)
obj <- ScaleData(obj, features = wgcna_genes, verbose = FALSE)

# group.by.vars = "condition": harmonise metacell eigengenes across conditions
obj <- ModuleEigengenes(obj, group.by.vars = "condition")

# ---- Step 7: Connectivity & hub genes -----------------------------------
obj <- ModuleConnectivity(obj,
                          group.by   = "celltype",
                          group_name = target_ct)

modules  <- GetModules(obj, wgcna_name = wgcna_name)
hub_genes <- GetHubGenes(obj, n_hubs = 10)
write.csv(modules,   out_path(cfg, "tables", "spoke_hdwgcna_modules.csv"),   row.names = FALSE)
write.csv(hub_genes, out_path(cfg, "tables", "spoke_hdwgcna_hub_genes.csv"), row.names = FALSE)
cat(sprintf("\nModules: %d  Hub genes (top 10 each): %d\n",
            length(unique(modules$module)), nrow(hub_genes)))

# KME distribution per module
n_mods  <- length(unique(modules$module[modules$module != "grey"]))
kme_ncol <- min(5L, n_mods)
tryCatch({
  p_kme <- PlotKMEs(obj, ncol = kme_ncol)
  ggsave(out_path(cfg, "plots", "spoke_hdwgcna_kme_distributions.pdf"),
         p_kme, width = 14, height = 8)
  ggsave(out_path(cfg, "plots", "spoke_hdwgcna_kme_distributions.png"),
         p_kme, width = 14, height = 8)
  rm(p_kme); invisible(gc())
}, error = function(e) message("[skip] KME plot: ", conditionMessage(e)))

# ---- Step 8: Module eigengene UMAP plots --------------------------------
tryCatch({
  plot_list_hme <- ModuleFeaturePlot(obj, features = "hMEs",
                                     order = TRUE, reduction = "umap")
  if (length(plot_list_hme) > 0) {
    p_hme <- wrap_plots(plot_list_hme, ncol = min(4L, length(plot_list_hme)))
    ggsave(out_path(cfg, "plots", "spoke_hdwgcna_module_eigengene_umap.pdf"),
           p_hme, width = 16, height = 12)
    ggsave(out_path(cfg, "plots", "spoke_hdwgcna_module_eigengene_umap.png"),
           p_hme, width = 16, height = 12)
    rm(p_hme)
  }
  rm(plot_list_hme); invisible(gc())
}, error = function(e) message("[skip] module eigengene UMAP: ", conditionMessage(e)))

# Module correlogram
tryCatch({
  pdf(out_path(cfg, "plots", "spoke_hdwgcna_module_correlogram.pdf"), width = 8, height = 7)
  ModuleCorrelogram(obj)
  dev.off()
  png(out_path(cfg, "plots", "spoke_hdwgcna_module_correlogram.png"), width = 800, height = 700)
  ModuleCorrelogram(obj)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message("[skip] correlogram: ", conditionMessage(e))
})

# ---- Step 9: Single-cell module scores (UCell) --------------------------
tryCatch({
  obj <- ModuleExprScore(obj, n_genes = 25, method = "UCell")
  plot_list_sc <- ModuleFeaturePlot(obj, features = "scores",
                                    order = TRUE, ucell = TRUE,
                                    reduction = "umap")
  if (length(plot_list_sc) > 0) {
    p_sc <- wrap_plots(plot_list_sc, ncol = min(4L, length(plot_list_sc)))
    ggsave(out_path(cfg, "plots", "spoke_hdwgcna_module_scores_umap.pdf"),
           p_sc, width = 16, height = 12)
    ggsave(out_path(cfg, "plots", "spoke_hdwgcna_module_scores_umap.png"),
           p_sc, width = 16, height = 12)
    rm(p_sc)
  }
  rm(plot_list_sc); invisible(gc())
}, error = function(e) message("[skip] UCell scores: ", conditionMessage(e)))

# ---- Step 10: Enrichment (enrichR) --------------------------------------
tryCatch({
  dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
  obj <- RunEnrichr(obj, dbs = dbs, max_genes = 100)
  enrich_tbl <- GetEnrichrTable(obj)
  write.csv(enrich_tbl,
            out_path(cfg, "tables", "spoke_hdwgcna_enrichment.csv"),
            row.names = FALSE)
  enrichr_dir <- file.path(cfg$project$outdir, "plots", "hdwgcna_enrichr")
  dir.create(enrichr_dir, recursive = TRUE, showWarnings = FALSE)
  EnrichrBarPlot(obj, outdir = enrichr_dir, n_terms = 10,
                 plot_size = c(5, 7), logscale = TRUE)
  message("Enrichment: ", nrow(enrich_tbl), " terms across all modules")
}, error = function(e) message("[skip] enrichR: ", conditionMessage(e)))

# ---- Step 11: Differential module eigengenes (DME) ----------------------
tryCatch({
  primary_grp <- cfg$de_contrasts[[1]]$group   # e.g. "Post"
  primary_ref <- cfg$de_contrasts[[1]]$ref     # e.g. "Healthy"
  bcs_grp <- colnames(obj)[obj$condition == primary_grp]
  bcs_ref <- colnames(obj)[obj$condition == primary_ref]
  if (length(bcs_grp) >= 10 && length(bcs_ref) >= 10) {
    DMEs <- FindDMEs(obj,
                     barcodes1  = bcs_grp,
                     barcodes2  = bcs_ref,
                     test.use   = "wilcox",
                     wgcna_name = wgcna_name)
    write.csv(DMEs,
              out_path(cfg, "tables", "spoke_hdwgcna_DMEs.csv"),
              row.names = FALSE)
    cat(sprintf("\nDME results (%s vs %s):\n", primary_grp, primary_ref))
    print(DMEs[order(DMEs$p_val_adj), ])

    p_dme <- PlotDMEsVolcano(obj, DMEs = DMEs, wgcna_name = wgcna_name)
    ggsave(out_path(cfg, "plots", "spoke_hdwgcna_DME_volcano.pdf"),
           p_dme, width = 8, height = 6)
    ggsave(out_path(cfg, "plots", "spoke_hdwgcna_DME_volcano.png"),
           p_dme, width = 8, height = 6)
    rm(p_dme); invisible(gc())
  } else {
    message("DME skipped: insufficient cells in one or both conditions")
  }
}, error = function(e) message("[skip] DME analysis: ", conditionMessage(e)))

# ---- Step 12: Network visualization -------------------------------------
# Hub gene network (memory-intensive: 8 GB limit set above)
tryCatch({
  pdf(out_path(cfg, "plots", "spoke_hdwgcna_hub_gene_network.pdf"),
      width = 10, height = 10)
  HubGeneNetworkPlot(obj, n_hubs = 3, n_other = 5,
                     edge_prop = 0.75, mods = "all")
  dev.off()
  png(out_path(cfg, "plots", "spoke_hdwgcna_hub_gene_network.png"),
      width = 1000, height = 1000)
  HubGeneNetworkPlot(obj, n_hubs = 3, n_other = 5,
                     edge_prop = 0.75, mods = "all")
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message("[skip] hub network: ", conditionMessage(e))
})

# Gene network UMAP
tryCatch({
  obj <- RunModuleUMAP(obj, n_hubs = 10, n_neighbors = 15, min_dist = 0.1)
  pdf(out_path(cfg, "plots", "spoke_hdwgcna_gene_network_umap.pdf"),
      width = 10, height = 9)
  ModuleUMAPPlot(obj, edge_prop = 0.05, label_hubs = 5,
                 keep_grey_edges = FALSE)
  dev.off()
  png(out_path(cfg, "plots", "spoke_hdwgcna_gene_network_umap.png"),
      width = 1000, height = 900)
  ModuleUMAPPlot(obj, edge_prop = 0.05, label_hubs = 5,
                 keep_grey_edges = FALSE)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message("[skip] gene network UMAP: ", conditionMessage(e))
})

# ---- Save ---------------------------------------------------------------
saveRDS(obj, out_path(cfg, "objects", "spoke_hdwgcna.rds"))
log_session(cfg, "spoke_hdwgcna")
