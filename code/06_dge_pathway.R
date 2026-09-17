#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 06_dge_pathway.R — pseudobulk DE + proportion test + pathway enrichment.
#                                                       (Tutorial Part 5)
# Inputs  : results/objects/05_annotated.rds
# Outputs : results/tables/06_de_<contrast>.csv, 06_fgsea_<contrast>.csv,
#           06_proportions.csv
# Run     : Rscript code/06_dge_pathway.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(msigdbr); library(Matrix); library(speckle)
})
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# ===== (A) cell-type proportion testing across groups (Part 5) =====
# propeller correctly accounts for sample-level variation in composition.
prop <- propeller(clusters = obj$celltype, sample = obj$sample, group = obj$condition)
write.csv(prop, out_path(cfg, "tables", "06_proportions.csv"))

# ===== (B) pseudobulk DE (the rigorous way to compare GROUPS) =====
# Sum raw counts per celltype..condition..sample (manual = no name munging).
obj$pb <- paste(obj$celltype, obj$condition, obj$sample, sep = "..")
counts <- GetAssayData(obj, layer = "counts")
grp    <- factor(obj$pb)
mm     <- sparse.model.matrix(~ 0 + grp); colnames(mm) <- levels(grp)
pb     <- counts %*% mm
meta   <- as.data.frame(do.call(rbind, strsplit(colnames(pb), "\\.\\.")))
colnames(meta) <- c("celltype", "condition", "sample"); rownames(meta) <- colnames(pb)

run_contrast <- function(ct, grp_lvl, ref_lvl) {
  keep <- meta$celltype == ct & meta$condition %in% c(grp_lvl, ref_lvl)
  cd   <- meta[keep, , drop = FALSE]
  if (min(table(cd$condition)) < 2) return(NULL)        # need >=2 donors/side
  m  <- round(as.matrix(pb[, rownames(cd)]))
  cd$condition <- relevel(factor(cd$condition), ref = ref_lvl)
  dds <- DESeqDataSetFromMatrix(m, cd, ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  res <- as.data.frame(results(DESeq(dds, quiet = TRUE),
                               contrast = c("condition", grp_lvl, ref_lvl)))
  res$gene <- rownames(res); res$celltype <- ct; res
}

# ===== (C) pathway enrichment (fgsea on Hallmark) per contrast =====
sp <- if (cfg$project$organism == "mouse") "Mus musculus" else "Homo sapiens"
hm <- msigdbr(species = sp, category = "H")
pathways  <- split(hm$gene_symbol, hm$gs_name)
celltypes <- unique(obj$celltype)

for (cset in cfg$de_contrasts) {
  de <- do.call(rbind, lapply(celltypes, run_contrast, cset$group, cset$ref))
  if (is.null(de)) { message("No estimable cell types for ", cset$name); next }
  write.csv(de, out_path(cfg, "tables", paste0("06_de_", cset$name, ".csv")), row.names = FALSE)

  ct1 <- names(sort(table(de$celltype[!is.na(de$padj) & de$padj < 0.05]), decreasing = TRUE))[1]
  if (!is.na(ct1)) {
    d <- subset(de, celltype == ct1 & !is.na(stat))
    g <- fgsea(pathways, setNames(d$stat, d$gene), minSize = 10, maxSize = 500)
    write.csv(apply(g, 2, as.character),
              out_path(cfg, "tables", paste0("06_fgsea_", cset$name, ".csv")), row.names = FALSE)
  }
  message("Done contrast: ", cset$name)
}
log_session(cfg, "06_dge_pathway")
