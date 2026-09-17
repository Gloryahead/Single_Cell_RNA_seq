# =====================================================================
# install_github_packages.R
# Run this INSIDE the activated `scrna` env, after `micromamba activate scrna`,
# to install the packages that are NOT reliably available on conda/bioconda.
#   R --no-save < install_github_packages.R
# remotes / devtools / BiocManager come from the conda env already.
# =====================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

gh <- function(repo, ...) remotes::install_github(repo, upgrade = "never", ...)

# ---- Part 3: Seurat integration wrappers ----
gh("satijalab/seurat-wrappers")

# ---- Part 2-2: DoubletFinder (scDblFinder is already in the conda env) ----
gh("chris-mcginnis-ucsf/DoubletFinder")

# ---- Part 7: Monocle 3 (github-only; leidenbase must come first) ----
gh("cole-trapnell-lab/leidenbase")
gh("cole-trapnell-lab/monocle3")

# ---- Part 9: CellChat ----
gh("jinworks/CellChat")

# ---- Part 10: NicheNet (network files still need downloading from Zenodo) ----
gh("saeyslab/nichenetr")
# Optional multi-sample extension:
# gh("saeyslab/multinichenetr")

# ---- Part 11: CopyKAT ----
gh("navinlabcode/copykat")

# ---- Part 12: hdWGCNA (dev branch) + up-to-date enrichR ----
gh("smorabit/hdWGCNA", ref = "dev")
# enrichR now comes from the conda env (r-enrichr). Only uncomment the next
# line if you specifically want the newer GitHub build the tutorial mentions:
# gh("wjawaid/enrichR")

# ---- Part 15: scplotter (try CRAN first, fall back to GitHub) ----
if (!requireNamespace("scplotter", quietly = TRUE)) {
  tryCatch(install.packages("scplotter"),
           error = function(e) gh("pwwang/scplotter"))
}

message("Done. Restart R and library() each package to verify.")
