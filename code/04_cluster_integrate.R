#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 04_cluster_integrate.R — filtered -> integrated, clustered object.
#                                                       (Tutorial Part 3)
# Inputs  : results/objects/03_filtered.rds
# Outputs : results/objects/04_clustered.rds + UMAPs in results/plots/
# Run     : Rscript code/04_cluster_integrate.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
cfg <- load_config(); set.seed(cfg$seed)
cc  <- cfg$cluster

obj <- readRDS(obj_path(cfg, "03_filtered"))

# ---- standard normalization -> HVG -> scale -> PCA ----
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = cc$n_hvg)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = cc$n_pcs)

# ---- batch integration (Part 3): only correct nuisance batch, never condition ----
reduction <- "pca"
if (cc$integration == "harmony") {
  suppressPackageStartupMessages(library(harmony))
  obj <- RunHarmony(obj, group.by.vars = cc$batch_key)
  reduction <- "harmony"
}

# ---- batch-mixing metric: kNN sample entropy (Harmony only) -------------
# For each cell, compute the normalised Shannon entropy of the sample
# distribution in its 20 nearest Harmony neighbours.
# Score 1.0 = perfect mixing across samples; 0.0 = all neighbours same sample.
# Scores below ~0.7 suggest residual batch structure worth investigating.
if (cc$integration == "harmony") {
  suppressPackageStartupMessages(library(FNN))
  emb_all   <- Embeddings(obj, "harmony")[, seq_len(cc$n_pcs), drop = FALSE]
  k_mix     <- 20L
  nn_idx    <- FNN::get.knnx(emb_all, emb_all, k = k_mix + 1L)$nn.index[, -1L]
  samp_vec  <- obj$sample
  n_samples <- length(unique(samp_vec))

  mixing_score <- apply(nn_idx, 1L, function(idx) {
    p <- as.numeric(table(samp_vec[idx])) / k_mix
    p <- p[p > 0]
    -sum(p * log(p)) / log(n_samples)  # normalised 0 – 1
  })
  obj$knn_mixing_score <- mixing_score

  mix_summ <- aggregate(knn_mixing_score ~ sample + condition,
                        data = obj@meta.data, FUN = mean)
  mix_summ$knn_mixing_score <- round(mix_summ$knn_mixing_score, 4)
  write.csv(mix_summ,
            out_path(cfg, "qc", "integration_mixing_score.csv"),
            row.names = FALSE)

  overall_mean <- round(mean(mixing_score), 4)
  cat(sprintf("\n=== Harmony batch-mixing (kNN entropy, k=%d) ===\n", k_mix))
  cat(sprintf("Overall mean mixing score: %.4f  (1.00 = perfect, <0.70 = investigate)\n",
              overall_mean))
  print(mix_summ[order(mix_summ$sample), ])

  p_mix <- ggplot2::ggplot(
    obj@meta.data,
    ggplot2::aes(x = sample, y = knn_mixing_score, fill = condition)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.3) +
    ggplot2::geom_hline(yintercept = overall_mean,
                        linetype = "dashed", color = "grey40") +
    ggplot2::annotate("text", x = Inf, y = overall_mean + 0.01,
                      label = paste0("mean=", overall_mean),
                      hjust = 1.1, size = 3, color = "grey30") +
    ggplot2::labs(x = NULL, y = "kNN mixing score (0–1)", fill = "Condition",
                  title = sprintf("Harmony integration quality (k=%d nearest neighbours)",
                                  k_mix)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title  = ggplot2::element_text(hjust = 0.5))
  ggsave(out_path(cfg, "plots", "04_harmony_mixing_score.png"),
         p_mix, width = 10, height = 5)
  ggsave(out_path(cfg, "plots", "04_harmony_mixing_score.pdf"),
         p_mix, width = 10, height = 5)
  rm(emb_all, nn_idx, mix_summ, p_mix); invisible(gc())
}

# ---- neighbors -> clusters -> UMAP, all on the integrated reduction ----
dims <- 1:cc$n_pcs
obj <- FindNeighbors(obj, reduction = reduction, dims = dims)

# ---- resolution sweep: select best by mean silhouette on 5000-cell subsample ----
suppressPackageStartupMessages(library(cluster))
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)
n_sil       <- min(5000L, ncol(obj))
set.seed(cfg$seed)
sil_cells <- sample(colnames(obj), n_sil)
emb_sub   <- Embeddings(obj, reduction)[sil_cells, seq_len(cc$n_pcs), drop = FALSE]
dist_sub  <- dist(emb_sub)

sil_rows <- list()
for (res in resolutions) {
  obj        <- FindClusters(obj, resolution = res, verbose = FALSE)
  labels_sub <- as.integer(Idents(obj)[sil_cells])
  n_cl       <- length(unique(labels_sub))
  if (n_cl < 2L) {
    sil_rows[[length(sil_rows) + 1L]] <- data.frame(resolution = res, n_clusters = n_cl,
                                                      mean_silhouette = NA_real_)
    next
  }
  sil <- silhouette(labels_sub, dist_sub)
  sil_rows[[length(sil_rows) + 1L]] <- data.frame(
    resolution      = res,
    n_clusters      = n_cl,
    mean_silhouette = round(mean(sil[, "sil_width"]), 4)
  )
}
rm(emb_sub, dist_sub); invisible(gc())

sil_df <- do.call(rbind, sil_rows)
write.csv(sil_df, out_path(cfg, "qc", "clustering_silhouette.csv"), row.names = FALSE)
cat("\n=== Resolution sweep (silhouette on", n_sil, "cells) ===\n"); print(sil_df)

best_idx <- which.max(sil_df$mean_silhouette)
best_res <- sil_df$resolution[best_idx]
message(sprintf("Selected resolution: %.1f  (%d clusters, silhouette = %.4f)",
                best_res, sil_df$n_clusters[best_idx], sil_df$mean_silhouette[best_idx]))
obj <- FindClusters(obj, resolution = best_res, verbose = FALSE)

obj <- RunUMAP(obj, reduction = reduction, dims = dims)
saveRDS(obj, obj_path(cfg, "04_clustered"))

# ---- cluster quality checks ----
clust_samp <- table(cluster = Idents(obj), sample = obj$sample)
clust_pct  <- round(100 * sweep(clust_samp, 1L, rowSums(clust_samp), "/"), 1)
write.csv(as.data.frame.matrix(clust_pct),
          out_path(cfg, "tables", "04_cluster_sample_pct.csv"))
dominated <- sort(apply(clust_pct, 1L, max), decreasing = TRUE)
cat("\n=== Cluster quality check ===\n")
cat("Max single-sample % per cluster (>70 = possible batch artefact):\n"); print(dominated)
small_n <- ncol(obj) * 0.01
flagged <- names(which(as.vector(table(Idents(obj))) < small_n))
cat(sprintf("Clusters with <1%% of cells (n < %.0f): %s\n", small_n,
            if (length(flagged)) paste(flagged, collapse = ", ") else "none"))

# ---- condition composition per cluster (stacked bar, PNG + PDF) ----
cond_levels  <- levels(factor(obj$condition))
cond_palette <- setNames(
  c("#2196F3", "#FF9800", "#F44336", "#9C27B0", "#4CAF50")[seq_along(cond_levels)],
  cond_levels
)
cond_df <- as.data.frame(table(cluster = Idents(obj), condition = obj$condition))
cond_df$pct <- ave(cond_df$Freq, cond_df$cluster, FUN = function(x) 100 * x / sum(x))
p_cond <- ggplot2::ggplot(cond_df, ggplot2::aes(x = cluster, y = pct, fill = condition)) +
  ggplot2::geom_col(position = "stack", width = 0.8) +
  ggplot2::scale_fill_manual(values = cond_palette) +
  ggplot2::labs(x = "Cluster", y = "% cells", fill = "Condition",
                title = "Condition composition per cluster") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 plot.title  = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "04_cluster_condition_dist.png"), p_cond, width = 12, height = 5)
ggsave(out_path(cfg, "plots", "04_cluster_condition_dist.pdf"), p_cond, width = 12, height = 5)

# sanity check: cells should mix by SAMPLE but separate by biology
p1 <- DimPlot(obj, label = TRUE, repel = TRUE, label.size = 3.5, pt.size = 0.3) +
      NoLegend() +
      ggtitle("Clusters") +
      theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(obj, group.by = "sample", pt.size = 0.3) +
      ggtitle("Sample mixing") +
      theme(plot.title = element_text(hjust = 0.5))
p_clusters <- p1 + p2
ggsave(out_path(cfg, "plots", "04_umap_clusters_samples.png"), p_clusters, width = 13, height = 5)
ggsave(out_path(cfg, "plots", "04_umap_clusters_samples.pdf"), p_clusters, width = 13, height = 5)
log_session(cfg, "04_cluster_integrate")
