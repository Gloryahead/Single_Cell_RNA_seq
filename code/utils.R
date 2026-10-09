# ---------------------------------------------------------------------
# utils.R — shared helpers. Sourced by every R script after setwd(PROJECT).
# Keeping repeated logic + paths here is what makes the pipeline reusable.
# ---------------------------------------------------------------------
suppressPackageStartupMessages({ library(yaml); library(Seurat); library(ggplot2) })

# ---- config + path helpers -------------------------------------------
load_config <- function(path = "config.yaml") yaml::read_yaml(path)

# build a results path (absolute, from config) and ensure its folder exists
out_path <- function(cfg, ...) {
  p <- file.path(cfg$project$outdir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
obj_path <- function(cfg, stage) out_path(cfg, "objects", paste0(stage, ".rds"))
mito_pattern <- function(cfg) if (cfg$project$organism == "mouse") "^mt-" else "^MT-"

# write exact package versions for reproducibility — call at end of each script
log_session <- function(cfg, stage) {
  writeLines(capture.output(sessionInfo()),
             out_path(cfg, "logs", paste0(stage, "_sessionInfo.txt")))
}

# ---- QC helpers (Part 2) ---------------------------------------------
add_qc_metrics <- function(obj, cfg) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern(cfg))
  obj
}
apply_qc_filter <- function(obj, cfg) {
  q <- cfg$qc
  subset(obj, subset = nFeature_RNA > q$min_genes &
                       nFeature_RNA < q$max_genes &
                       percent.mt   < q$max_mito_pct)
}
# scDblFinder doublet flagging (Part 2-2) -> adds $doublet column
run_scdblfinder <- function(obj) {
  suppressPackageStartupMessages({ library(scDblFinder); library(SingleCellExperiment) })
  sce <- as.SingleCellExperiment(obj)
  has_s <- "sample" %in% colnames(SummarizedExperiment::colData(sce))
  sce <- scDblFinder(sce, samples = if (has_s) "sample" else NULL)
  obj$doublet <- sce$scDblFinder.class
  obj
}
# SoupX ambient-RNA correction (Part 2-2) -> corrected counts matrix
run_soupx <- function(filtered_dir, raw_dir) {
  suppressPackageStartupMessages({ library(SoupX) })
  filt <- Read10X(filtered_dir); raw <- Read10X(raw_dir)
  sc <- SoupChannel(raw, filt)
  tmp <- CreateSeuratObject(filt) |> NormalizeData() |> FindVariableFeatures() |>
           ScaleData() |> RunPCA() |> FindNeighbors() |> FindClusters() |> suppressWarnings()
  sc <- setClusters(sc, Idents(tmp)); sc <- autoEstCont(sc)
  adjustCounts(sc)
}

# ---- export hub metadata for the Python velocity/fate spokes -----------
# Tutorial approach (ngs101.com Part 13): write two small CSVs from R;
# Python builds AnnData from STARsolo matrices and joins these.
# Avoids sceasy/reticulate which duplicates the full object in memory.
export_meta_for_velocity <- function(obj, cfg) {
  meta <- obj@meta.data
  meta$barcode <- rownames(meta)
  write.csv(meta, out_path(cfg, "objects", "05_annotated_meta.csv"), row.names = FALSE)

  umap <- as.data.frame(Embeddings(obj, reduction = "umap"))
  umap$barcode <- rownames(umap)
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  write.csv(umap, out_path(cfg, "objects", "05_annotated_umap.csv"), row.names = FALSE)
}
