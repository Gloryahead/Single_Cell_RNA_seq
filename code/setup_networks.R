#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# setup_networks.R — Pre-download OmniPath network files to the results cache.
#
# WHY: Puma compute nodes have no outbound HTTP. Parts 16 and 17 try to
# download CollecTRI and PROGENy at runtime; if that fails they fall back
# to cached CSVs. Run this script ONCE from a LOGIN NODE (which has HTTP)
# before submitting any pipeline jobs, and the CSVs will be ready.
#
# Usage (login node, inside container):
#   apptainer exec $DATA/containers/scrna_r.sif \
#     Rscript code/setup_networks.R
#
# Or with the active conda env:
#   conda run -n scrna Rscript code/setup_networks.R
# ---------------------------------------------------------------------
PROJECT <- Sys.getenv("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
setwd(PROJECT); source("code/utils.R")
suppressPackageStartupMessages(library(decoupleR))
cfg <- load_config()

ref_dir <- file.path(cfg$project$outdir, "references")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
cat("Network cache directory:", ref_dir, "\n\n")

# ---- CollecTRI (Part 16: TF activity via run_ulm) ------------------------
collectri_csv <- file.path(ref_dir, "collectri_human.csv")
cat("1/2  CollecTRI (human) ...\n")
tryCatch({
  net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
  write.csv(net, collectri_csv, row.names = FALSE)
  cat(sprintf("     OK  %s\n     %d interactions | %d TFs\n",
              collectri_csv, nrow(net), length(unique(net$source))))
}, error = function(e) {
  cat("     FAILED:", conditionMessage(e), "\n")
  cat("     Part 16 will fall back to GTRD (MSigDB C3, no repression sign).\n")
})

# ---- PROGENy (Part 17: pathway activity via run_mlm) ---------------------
progeny_csv <- file.path(ref_dir, "progeny_human.csv")
cat("\n2/2  PROGENy (human, top=500) ...\n")
tryCatch({
  net <- decoupleR::get_progeny(organism = "human", top = 500)
  write.csv(net, progeny_csv, row.names = FALSE)
  cat(sprintf("     OK  %s\n     %d gene-pathway links | %d pathways\n",
              progeny_csv, nrow(net), length(unique(net$source))))
}, error = function(e) {
  cat("     FAILED:", conditionMessage(e), "\n")
  cat("     Part 17 has no fallback — the job will fail if this step is skipped.\n")
  cat("     Try again later or provide a pre-built progeny_human.csv manually.\n")
})

cat("\nDone. Both CSVs are read automatically by Parts 16 and 17 at job start.\n")
cat("Re-run this script to refresh the networks (e.g. after OmniPath releases an update).\n")
