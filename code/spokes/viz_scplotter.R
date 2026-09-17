#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# viz_scplotter.R — publication-ready figures (scplotter).        (Part 15)
# Install: install.packages("scplotter")  or  remotes::install_github("pwwang/scplotter")
# In: results/objects/05_annotated.rds   Out: results/plots/spoke_scplotter_*.png
# Run: Rscript code/spokes/viz_scplotter.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(scplotter))
cfg <- load_config(); set.seed(cfg$seed)
obj <- readRDS(obj_path(cfg, "05_annotated"))

# one-line, polished versions of the plots you'd otherwise hand-build:
ggplot2::ggsave(out_path(cfg, "plots", "spoke_scplotter_umap.png"),
                CellDimPlot(obj, group_by = "celltype"), width = 7, height = 6)
ggplot2::ggsave(out_path(cfg, "plots", "spoke_scplotter_props.png"),
                CellStatPlot(obj, ident = "celltype", group_by = "condition",
                             frac = "group"), width = 7, height = 5)
log_session(cfg, "spoke_scplotter")
