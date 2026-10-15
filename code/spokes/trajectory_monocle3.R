#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# trajectory_monocle3.R — trajectory & pseudotime via Monocle 3. (Part 7)
#
# NOTE: monocle3 and SeuratWrappers are NOT in the current scrna_r.sif.
#       Rebuild the container with:
#           BiocManager::install("monocle3")
#           devtools::install_github("satijalab/seurat-wrappers")
#       before enabling this rule in Snakefile rule all.
#
# HPC NOTE: order_cells() without root_cells opens an interactive
#           browser — NEVER call it that way on a headless HPC node.
#           This script always selects roots programmatically.
#
# In : results/objects/05_annotated.rds
# Out: results/objects/spoke_monocle3.rds
#      results/tables/spoke_monocle3_trajectory_genes.csv
#      results/plots/spoke_monocle3_*.{png,pdf}
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(monocle3); library(SeuratWrappers); library(patchwork)
})
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# ---- Optional downsample for speed (keeps proportions intact) -----------
max_cells <- 50000L
if (ncol(obj) > max_cells) {
  message(sprintf("Downsampling from %d to %d cells (stratified by cell type)",
                  ncol(obj), max_cells))
  set.seed(cfg$seed)
  ct_vec  <- obj$celltype
  keep    <- unlist(lapply(unique(ct_vec), function(ct) {
    idx <- which(ct_vec == ct)
    n   <- max(50L, round(length(idx) * max_cells / ncol(obj)))
    sample(idx, min(n, length(idx)))
  }))
  obj <- obj[, keep]
}
cat(sprintf("Monocle3 running on %d cells\n", ncol(obj)))

# ---- Seurat -> CDS conversion, copy Seurat UMAP -------------------------
cds <- as.cell_data_set(obj)

# Transfer Seurat UMAP into CDS so we analyse in the same embedding
umap_mat <- Embeddings(obj, "umap")
cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- umap_mat

# ---- Monocle3 preprocessing ---------------------------------------------
# norm_method="log" — skip normalisation since Seurat already normalised.
# num_dim=50 matches the Harmony PC space used for UMAP.
cds <- preprocess_cds(cds, num_dim = 50, norm_method = "log")

# ---- Clustering (UMAP-based, louvain) -----------------------------------
# k=20 gives stable community detection across PBMC datasets.
cds <- cluster_cells(cds, reduction_method = "UMAP",
                     k = 20, cluster_method = "louvain")
message("Monocle3 partitions: ", length(table(partitions(cds))))
message("Monocle3 clusters:   ", length(table(clusters(cds))))

# ---- Learn trajectory graph ---------------------------------------------
# use_partition = TRUE keeps biologically separate compartments
# (T cells vs monocytes vs B cells) in separate lineages.
# ncenter = 500 controls trajectory graph resolution.
cds <- learn_graph(cds, use_partition = TRUE,
                   close_loop = FALSE, ncenter = 500)

# ---- Programmatic root selection per partition --------------------------
# Priority: T_cells (naive) > T_cells > NK > first N cells in partition.
# Root selection is REQUIRED for headless HPC execution (no interactive UI).
get_roots <- function(cds, root_priority, n_cand = 50L) {
  parts   <- partitions(cds)
  types   <- colData(cds)$celltype
  roots   <- character(0)
  for (p in unique(parts)) {
    p_cells <- names(parts)[parts == p]
    found   <- NULL
    for (pat in root_priority) {
      hits <- p_cells[grepl(pat, types[p_cells], ignore.case = TRUE)]
      if (length(hits) >= 1) {
        found <- sample(hits, min(n_cand, length(hits)))
        break
      }
    }
    if (is.null(found)) found <- head(p_cells, n_cand)
    message(sprintf("  Partition %s root: %d cells (type pattern matched)", p, length(found)))
    roots <- c(roots, found)
  }
  roots
}

root_patterns <- c("naive.*T", "CCR7.*T", "T_cell", "T cell", "CD4", "CD8",
                   "NK_cell", "NK cell", "B_cell", "Monocyte")
message("Selecting trajectory roots:")
root_cells <- get_roots(cds, root_patterns)
cds <- order_cells(cds, root_cells = root_cells)

cat("\nPseudotime range by cell type:\n")
pt_by_ct <- tapply(cds$monocle3_pseudotime, colData(cds)$celltype,
                   function(x) round(median(x[is.finite(x)]), 2))
print(sort(pt_by_ct))

# ---- Trajectory-associated genes (graph test) --------------------------
message("Running graph_test to find trajectory-associated genes...")
pr_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
pr_res <- pr_res[!is.na(pr_res$q_value), ]
pr_res <- pr_res[order(pr_res$q_value, -pr_res$morans_I), ]
n_sig  <- sum(pr_res$q_value < 0.05, na.rm = TRUE)
cat(sprintf("\n%d trajectory-associated genes (q < 0.05)\n", n_sig))
cat("Top genes:\n"); print(head(pr_res[pr_res$q_value < 0.05,
                                       c("gene_short_name","morans_I","q_value")], 10))
write.csv(pr_res[pr_res$q_value < 0.05, ],
          out_path(cfg, "tables", "spoke_monocle3_trajectory_genes.csv"),
          row.names = FALSE)
top_genes <- head(pr_res$gene_short_name[pr_res$q_value < 0.05], 6)

# ---- Plots --------------------------------------------------------------
# Trajectory colored by cell type
png_w <- 12; png_h <- 6

p_traj <- plot_cells(cds,
                     color_cells_by   = "celltype",
                     label_groups_by_cluster = FALSE,
                     label_leaves     = FALSE,
                     label_branch_points = FALSE,
                     graph_label_size = 3,
                     cell_size        = 0.4,
                     trajectory_graph_color = "black") +
          ggplot2::ggtitle("Monocle3 trajectory — cell type") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

p_pt <- plot_cells(cds,
                   color_cells_by   = "monocle3_pseudotime",
                   label_cell_groups = FALSE,
                   label_leaves      = FALSE,
                   label_branch_points = FALSE,
                   cell_size         = 0.4,
                   trajectory_graph_color = "black") +
        ggplot2::scale_color_viridis_c(option = "plasma") +
        ggplot2::ggtitle("Monocle3 pseudotime") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggsave(out_path(cfg, "plots", "spoke_monocle3_trajectory.png"),
       p_traj + p_pt, width = png_w, height = png_h)
ggsave(out_path(cfg, "plots", "spoke_monocle3_trajectory.pdf"),
       p_traj + p_pt, width = png_w, height = png_h)

# Condition split
p_cond <- plot_cells(cds,
                     color_cells_by = "condition",
                     label_cell_groups = FALSE,
                     label_leaves   = FALSE,
                     label_branch_points = FALSE,
                     cell_size      = 0.4) +
          ggplot2::ggtitle("Monocle3 trajectory — condition") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "spoke_monocle3_trajectory_condition.png"),
       p_cond, width = 8, height = 6)
ggsave(out_path(cfg, "plots", "spoke_monocle3_trajectory_condition.pdf"),
       p_cond, width = 8, height = 6)

# Top trajectory genes on UMAP
if (length(top_genes) > 0) {
  p_genes <- plot_cells(cds,
                        genes     = top_genes,
                        show_trajectory_graph = FALSE,
                        label_cell_groups = FALSE,
                        cell_size = 0.3)
  ggsave(out_path(cfg, "plots", "spoke_monocle3_top_trajectory_genes.png"),
         p_genes, width = 14, height = 10)
  ggsave(out_path(cfg, "plots", "spoke_monocle3_top_trajectory_genes.pdf"),
         p_genes, width = 14, height = 10)
}
rm(p_traj, p_pt, p_cond); invisible(gc())

# ---- Write pseudotime back to original Seurat object -------------------
sub_cells   <- colnames(obj)[colnames(obj) %in% rownames(colData(cds))]
obj@meta.data[["monocle3_pseudotime"]] <- NA_real_
obj@meta.data[["monocle3_cluster"]]    <- NA_character_
m <- match(sub_cells, colnames(obj))
obj@meta.data[["monocle3_pseudotime"]][m] <- cds$monocle3_pseudotime[sub_cells]
obj@meta.data[["monocle3_cluster"]][m]    <- as.character(clusters(cds)[sub_cells])

# FeaturePlot pseudotime in the Seurat UMAP embedding
p_feat <- FeaturePlot(obj, "monocle3_pseudotime", pt.size = 0.3, order = TRUE) +
          ggplot2::scale_color_viridis_c(option = "plasma", na.value = "grey90") +
          ggplot2::ggtitle("Monocle3 pseudotime")
ggsave(out_path(cfg, "plots", "spoke_monocle3_in_seurat.png"),  p_feat, width = 7, height = 6)
ggsave(out_path(cfg, "plots", "spoke_monocle3_in_seurat.pdf"),  p_feat, width = 7, height = 6)

saveRDS(cds, out_path(cfg, "objects", "spoke_monocle3.rds"))
log_session(cfg, "spoke_monocle3")
