#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 06_dge_pathway.R — DE + proportions + pathway enrichment.
#                                                       (Tutorial Part 5)
# Inputs  : results/objects/05_annotated.rds
# Outputs : results/tables/06_*.csv, results/plots/06_*.{png,pdf}
# Run     : Rscript code/06_dge_pathway.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(msigdbr)
  library(Matrix); library(speckle); library(ggrepel)
})
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))
sp  <- if (cfg$project$organism == "mouse") "Mus musculus" else "Homo sapiens"

celltypes <- sort(unique(obj$celltype))

# =========================================================================
# (A) Cell-type proportions
# =========================================================================
prop_counts <- as.data.frame(
  table(celltype = obj$celltype, sample = obj$sample, condition = obj$condition))
prop_counts <- prop_counts[prop_counts$Freq > 0, ]
prop_counts <- do.call(rbind, lapply(split(prop_counts, prop_counts$sample), function(x) {
  x$pct <- 100 * x$Freq / sum(x$Freq); x
}))
write.csv(prop_counts, out_path(cfg, "tables", "06_proportions_per_sample.csv"),
          row.names = FALSE)

prop_summ <- aggregate(pct ~ celltype + condition, prop_counts,
                       FUN = function(x) c(mean = round(mean(x), 2), sd = round(sd(x), 2)))
prop_summ <- do.call(data.frame, prop_summ)
names(prop_summ)[3:4] <- c("mean_pct", "sd_pct")

ct_pal <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                     "#A65628","#F781BF","#999999","#66C2A5","#FC8D62"))(
    length(celltypes)), celltypes)

# Stacked barplot: composition per sample
p_stack <- ggplot2::ggplot(prop_counts,
                           ggplot2::aes(x = sample, y = pct, fill = celltype)) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::scale_fill_manual(values = ct_pal) +
  ggplot2::labs(x = NULL, y = "% cells", fill = "Cell type",
                title = "Cell type composition per sample") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 plot.title  = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "06_stacked_composition.png"), p_stack, width = 10, height = 5)
ggsave(out_path(cfg, "plots", "06_stacked_composition.pdf"), p_stack, width = 10, height = 5)

# Grouped barplot: mean ± SD per condition
p_grp <- ggplot2::ggplot(prop_summ,
                         ggplot2::aes(x = celltype, y = mean_pct, fill = condition)) +
  ggplot2::geom_col(position = "dodge", width = 0.7) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, mean_pct - sd_pct), ymax = mean_pct + sd_pct),
    position = ggplot2::position_dodge(0.7), width = 0.25) +
  ggplot2::labs(x = NULL, y = "Mean % cells ± SD", fill = "Condition",
                title = "Cell type proportions by condition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 plot.title  = ggplot2::element_text(hjust = 0.5))
ggsave(out_path(cfg, "plots", "06_grouped_proportions.png"), p_grp, width = 11, height = 5)
ggsave(out_path(cfg, "plots", "06_grouped_proportions.pdf"), p_grp, width = 11, height = 5)
rm(p_stack, p_grp); invisible(gc())

# Propeller test: accounts for sample-level composition variability
prop_test <- propeller(clusters = obj$celltype, sample = obj$sample, group = obj$condition)
write.csv(prop_test, out_path(cfg, "tables", "06_proportions.csv"))
cat("\n=== Propeller proportion test ===\n"); print(head(prop_test, 10))

# =========================================================================
# (B) Cell-level DE with FindMarkers (Wilcoxon, per cell type)
# =========================================================================
# logfc.threshold = 0.1 (not 0) to keep runtime feasible.
# Use pseudobulk stats (section C) for fgsea — not these cell-level stats.
volcano_pb <- function(df, x_col, y_col, fc_col, title, outfile) {
  df <- df[!is.na(df[[y_col]]), ]
  df$sig   <- df[[y_col]] < 0.05 & abs(df[[fc_col]]) > 0.25
  rnk      <- rank(-abs(df[[fc_col]]), ties.method = "first")
  df$label <- ifelse(df$sig & rnk <= 10, df$gene, NA_character_)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x_col]], -log10(.data[[y_col]] + 1e-300),
                                         color = sig, label = label)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.5) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#E41A1C")) +
    ggplot2::geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey60") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    ggplot2::labs(title = title, x = x_col, y = paste0("-log10(", y_col, ")")) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 10))
  ggsave(paste0(outfile, ".png"), p, width = 6, height = 5)
  ggsave(paste0(outfile, ".pdf"), p, width = 6, height = 5)
}

Idents(obj) <- "condition"
for (cset in cfg$de_contrasts) {
  message("Cell-level DE (FindMarkers): ", cset$name)
  ct_res <- lapply(celltypes, function(ct) {
    cells <- which(obj$celltype == ct)
    if (sum(obj$condition[cells] == cset$group) < 3 ||
        sum(obj$condition[cells] == cset$ref)   < 3) return(NULL)
    sub <- obj[, cells]; Idents(sub) <- "condition"
    tryCatch(
      FindMarkers(sub, ident.1 = cset$group, ident.2 = cset$ref,
                  logfc.threshold = 0.1, min.pct = 0.1,
                  max.cells.per.ident = 1000, verbose = FALSE),
      error = function(e) { message("  skipped ", ct, ": ", conditionMessage(e)); NULL }
    )
  })
  names(ct_res) <- celltypes
  ct_res <- Filter(Negate(is.null), ct_res)
  if (length(ct_res) == 0) next

  for (ct in names(ct_res)) {
    df <- ct_res[[ct]]; df$gene <- rownames(df)
    volcano_pb(df, "avg_log2FC", "p_val_adj", "avg_log2FC",
               paste0(ct, " — ", cset$name),
               out_path(cfg, "plots",
                        paste0("06_volcano_cell_", cset$name, "_",
                               gsub("[^A-Za-z0-9]", "_", ct))))
  }
  out_df <- do.call(rbind, lapply(names(ct_res), function(ct) {
    df <- ct_res[[ct]]; df$gene <- rownames(df); df$celltype <- ct; df
  }))
  write.csv(out_df,
            out_path(cfg, "tables", paste0("06_findmarkers_", cset$name, ".csv")),
            row.names = FALSE)
  n_sig <- sum(out_df$p_val_adj < 0.05 & abs(out_df$avg_log2FC) > 0.25, na.rm = TRUE)
  message("  ", cset$name, ": ", n_sig, " significant cell-level hits (FDR<0.05, |FC|>0.25)")
  rm(ct_res, out_df); invisible(gc())
}
Idents(obj) <- "celltype"

# =========================================================================
# (C) Pseudobulk DE with DESeq2 + volcano + fgsea per cell type
# =========================================================================
obj$pb  <- paste(obj$celltype, obj$condition, obj$sample, sep = "..")
counts  <- GetAssayData(obj, layer = "counts")
grp     <- factor(obj$pb)
mm      <- sparse.model.matrix(~ 0 + grp); colnames(mm) <- levels(grp)
pb      <- counts %*% mm
meta    <- as.data.frame(do.call(rbind, strsplit(colnames(pb), "\\.\\.")))
colnames(meta) <- c("celltype", "condition", "sample"); rownames(meta) <- colnames(pb)
rm(counts, grp, mm); invisible(gc())

run_deseq2 <- function(ct, grp_lvl, ref_lvl) {
  keep <- meta$celltype == ct & meta$condition %in% c(grp_lvl, ref_lvl)
  cd   <- meta[keep, , drop = FALSE]
  if (min(table(factor(cd$condition, levels = c(grp_lvl, ref_lvl)))) < 2) return(NULL)
  m  <- round(as.matrix(pb[, rownames(cd)]))
  cd$condition <- relevel(factor(cd$condition), ref = ref_lvl)
  dds <- DESeqDataSetFromMatrix(m, cd, ~ condition)
  # keep genes with ≥10 counts in at least 3 samples
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  if (nrow(dds) < 10) return(NULL)
  res        <- as.data.frame(results(DESeq(dds, quiet = TRUE),
                                      contrast = c("condition", grp_lvl, ref_lvl)))
  res$gene   <- rownames(res); res$celltype <- ct; res
}

hm   <- msigdbr(species = sp, collection = "H")
pathways <- split(hm$gene_symbol, hm$gs_name)

for (cset in cfg$de_contrasts) {
  message("Pseudobulk DE: ", cset$name)
  de_list <- lapply(celltypes, run_deseq2, cset$group, cset$ref)
  names(de_list) <- celltypes
  de_list <- Filter(Negate(is.null), de_list)
  if (length(de_list) == 0) { message("  no estimable cell types"); next }

  de <- do.call(rbind, de_list)
  write.csv(de, out_path(cfg, "tables", paste0("06_de_", cset$name, ".csv")),
            row.names = FALSE)
  n_sig <- sum(de$padj < 0.05 & abs(de$log2FoldChange) > 0.25, na.rm = TRUE)
  message("  ", n_sig, " significant pseudobulk hits across all cell types")

  # Volcano per cell type
  for (ct in names(de_list)) {
    d <- de_list[[ct]]
    d$gene <- rownames(d)
    volcano_pb(d, "log2FoldChange", "padj", "log2FoldChange",
               paste0(ct, " — ", cset$name, " (pseudobulk)"),
               out_path(cfg, "plots",
                        paste0("06_volcano_pb_", cset$name, "_",
                               gsub("[^A-Za-z0-9]", "_", ct))))
  }

  # fgsea on every cell type with ≥50 finite stats
  fgsea_rows <- list()
  for (ct in names(de_list)) {
    d <- de_list[[ct]]
    rnk <- setNames(d$stat, d$gene)
    rnk <- rnk[is.finite(rnk)]
    if (length(rnk) < 50) next
    g <- tryCatch(
      fgsea(pathways, rnk, minSize = 10, maxSize = 500),
      error = function(e) NULL
    )
    if (is.null(g) || nrow(g) == 0) next
    g$celltype <- ct
    fgsea_rows[[ct]] <- g
  }
  if (length(fgsea_rows) > 0) {
    fgsea_out <- do.call(rbind, lapply(fgsea_rows, function(x)
      as.data.frame(lapply(x, function(col)
        if (is.list(col)) sapply(col, paste, collapse = ";") else col))))
    write.csv(fgsea_out,
              out_path(cfg, "tables", paste0("06_fgsea_", cset$name, ".csv")),
              row.names = FALSE)
    n_pw <- sum(fgsea_out$padj < 0.05, na.rm = TRUE)
    message("  fgsea: ", n_pw, " significant pathways across cell types")
  }
  rm(de_list, de, fgsea_rows); invisible(gc())
  message("Done contrast: ", cset$name)
}

# =========================================================================
# (D) Functional state scoring: HALLMARK_INFLAMMATORY_RESPONSE
# =========================================================================
inflam <- pathways[["HALLMARK_INFLAMMATORY_RESPONSE"]]
inflam <- intersect(inflam, rownames(obj))
if (length(inflam) > 0) {
  obj <- AddModuleScore(obj, features = list(inflam), name = "Inflam", nbin = 24, ctrl = 100)
  names(obj@meta.data)[ncol(obj@meta.data)] <- "Inflammatory_Response_Score"

  # Ridge plot: distribution per cell type × condition
  tryCatch({
    suppressPackageStartupMessages(library(ggridges))
    p_ridge <- ggplot2::ggplot(
      obj@meta.data,
      ggplot2::aes(x = Inflammatory_Response_Score, y = celltype,
                   fill = condition)) +
      ggridges::geom_density_ridges(alpha = 0.65, scale = 1.2) +
      ggplot2::labs(x = "Inflammatory Response Score", y = NULL,
                    fill = "Condition", title = "HALLMARK_INFLAMMATORY_RESPONSE") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    ggsave(out_path(cfg, "plots", "06_ridge_inflammatory_response.png"),
           p_ridge, width = 9, height = 6)
    ggsave(out_path(cfg, "plots", "06_ridge_inflammatory_response.pdf"),
           p_ridge, width = 9, height = 6)
    rm(p_ridge); invisible(gc())
  }, error = function(e) message("ggridges not available: ", conditionMessage(e)))

  # Pseudobulk t-test per cell type: mean score Post vs Healthy
  score_test <- do.call(rbind, lapply(celltypes, function(ct) {
    d <- obj@meta.data[obj$celltype == ct, ]
    h <- d$Inflammatory_Response_Score[d$condition == "Healthy"]
    p <- d$Inflammatory_Response_Score[d$condition == "Post"]
    if (length(h) < 3 || length(p) < 3) return(NULL)
    tt <- t.test(p, h)
    data.frame(celltype = ct, mean_healthy = round(mean(h), 4),
               mean_post = round(mean(p), 4),
               diff = round(mean(p) - mean(h), 4),
               p_value = tt$p.value, stringsAsFactors = FALSE)
  }))
  if (!is.null(score_test) && nrow(score_test) > 0) {
    score_test$FDR <- p.adjust(score_test$p_value, method = "BH")
    write.csv(score_test,
              out_path(cfg, "tables", "06_inflammatory_score_test.csv"),
              row.names = FALSE)
    cat("\n=== Inflammatory Response score: Post vs Healthy ===\n")
    print(score_test[order(score_test$FDR), ])
  }
  rm(inflam); invisible(gc())
}

log_session(cfg, "06_dge_pathway")
