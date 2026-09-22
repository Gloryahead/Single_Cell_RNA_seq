#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# diagnose_labels.R — settle the CD8-vs-"T cells" question with NUMBERS.
# Read-only: inspects 05_annotated.rds, changes nothing.
# Run: Rscript code/diagnose_labels.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
cfg <- load_config()
obj <- readRDS(obj_path(cfg, "05_annotated"))

# ---- 1) how many cells per label? tiny labels are usually artifacts ----
cat("\n=== cells per celltype ===\n"); print(sort(table(obj$celltype)))
cat("\n=== celltype x condition ===\n"); print(table(obj$celltype, obj$condition))

# ---- 2) mean expression of the DEFINITIVE lineage genes per label ----
# CD8A/CD8B  = cytotoxic T   |  CD4/IL7R = helper T
# GZMK/GNLY/NKG7 = cytotoxic effector  |  PPBP/PF4 = platelet
genes <- c("CD3D","CD3E","CD4","IL7R","CD8A","CD8B",
           "GZMK","GZMB","GNLY","NKG7","PPBP","PF4","ITGA2B")
genes <- genes[genes %in% rownames(obj)]           # drop any absent from the assay

avg <- AverageExpression(obj, features = genes, group.by = "celltype", layer = "data")$RNA
cat("\n=== mean normalized expression (rows=gene, cols=celltype) ===\n")
print(round(as.matrix(avg), 2))

# ---- 3) percent of cells expressing each gene, per label ----
pct <- sapply(levels(factor(obj$celltype)), function(ct) {
  cells <- colnames(obj)[obj$celltype == ct]
  m <- GetAssayData(obj, layer = "counts")[genes, cells, drop = FALSE]
  round(100 * Matrix::rowMeans(m > 0), 1)
})
cat("\n=== % of cells expressing (rows=gene, cols=celltype) ===\n")
print(pct)

cat("\n---------------- how to read this ----------------\n",
    "* A real CD8 T population: CD3D+ CD8A/CD8B high, CD4/IL7R low.\n",
    "* A real CD4 T population: CD3D+ IL7R/CD4 high, CD8A/CD8B low.\n",
    "* If your 'CD8+ T cells' show HIGH IL7R and LOW CD8A/CD8B, and your\n",
    "  generic 'T cells' show HIGH CD8A/GZMK/NKG7 -> the labels are SWAPPED.\n",
    "* 'Progenitors' with high PPBP/PF4/ITGA2B are PLATELETS, not progenitors.\n",
    "* Labels with <50 cells are usually artifacts (esp. Neutrophils in PBMC).\n")
