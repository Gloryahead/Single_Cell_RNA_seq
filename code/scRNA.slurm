#!/bin/bash

# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------

#SBATCH --job-name=scRNA-seq
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=120:00:00

# --------------------------------------------------------------
### PART 2: Environment setup
# --------------------------------------------------------------

echo "========== Job started: $(date) =========="

# Code (scripts, config, envs) lives in home; data and results live on xdisk.
export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq
cd "$SCRNA_PROJECT"            # run from code root so ./code/... paths resolve
mkdir -p "$SCRNA_DATA/results/logs"

# Load micromamba shell functions, then switch to the haining group prefix
# where the `scrna` (R) and `scrna-velocity` (Python) environments live.
set +eu; source ~/.bashrc; set -eu
mamba-haining   # sets MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba

# Cell Ranger: installed once via code/00b_install_cellranger.sh
source "$SCRNA_DATA/software/activate_cellranger.sh" 2>/dev/null || \
  echo "WARNING: activate_cellranger.sh not found — run 00b_install_cellranger.sh first"

command -v cellranger >/dev/null && cellranger --version || \
  { echo "FATAL: cellranger not on PATH"; exit 1; }

# --------------------------------------------------------------
### PART 3: Helper functions
# --------------------------------------------------------------

run_if_exists() {
    local SCRIPT=$1
    if [[ -f "$SCRIPT" ]]; then
        echo "---------- Running: $(basename $SCRIPT) ($(date)) ----------"
        bash "$SCRIPT"
    else
        echo "WARNING: Script not found: $SCRIPT"
    fi
}

run_R_if_exists() {
    local SCRIPT=$1
    if [[ -f "$SCRIPT" ]]; then
        echo "---------- Running R: $(basename $SCRIPT) ($(date)) ----------"
        SCRNA_PROJECT="$SCRNA_PROJECT" SCRNA_DATA="$SCRNA_DATA" \
          micromamba run -n scrna Rscript "$SCRIPT"
    else
        echo "WARNING: Script not found: $SCRIPT"
    fi
}

run_PY_if_exists() {
    local SCRIPT=$1
    if [[ -f "$SCRIPT" ]]; then
        echo "---------- Running Python: $(basename $SCRIPT) ($(date)) ----------"
        SCRNA_PROJECT="$SCRNA_PROJECT" SCRNA_DATA="$SCRNA_DATA" \
          micromamba run -n scrna-velocity python "$SCRIPT"
    else
        echo "WARNING: Script not found: $SCRIPT"
    fi
}

# --------------------------------------------------------------
### PART 4: Run the pipeline
# --------------------------------------------------------------

# ── One-time setup (comment out after first run) ──────────────
# run_if_exists ./code/00_setup_dirs.sh

# ── Part 1: Download FASTQ + FastQC + Cell Ranger alignment ──
# These are compute-heavy; submit as separate sbatch jobs if preferred.
# run_if_exists ./code/01_download_fastq.sh
# run_if_exists ./code/01b_fastqc.sh
# run_if_exists ./code/02_align_cellranger.sh

# ── Part 2: QC filtering ──────────────────────────────────────
run_R_if_exists ./code/03_qc_filter.R

# ── Part 3: Integration & clustering ─────────────────────────
run_R_if_exists ./code/04_cluster_integrate.R

# ── Part 4: Cell type annotation ─────────────────────────────
run_R_if_exists ./code/05_annotate.R

# ── Part 5: DE + proportion testing + pathway enrichment ─────
run_R_if_exists ./code/06_dge_pathway.R

# ── Part 7: Trajectory (Monocle 3 + Slingshot) ───────────────
run_R_if_exists ./code/spokes/trajectory_monocle3.R
run_R_if_exists ./code/spokes/trajectory_slingshot.R

# ── Part 8: PDX deconvolution (skip if not PDX data) ─────────
# run_if_exists ./code/alt/08_pdx_human_mouse.sh

# ── Part 9: Cell-cell communication (CellChat) ───────────────
run_R_if_exists ./code/spokes/ccc_cellchat.R

# ── Part 10: Ligand-target prediction (NicheNet) ─────────────
run_R_if_exists ./code/spokes/ccc_nichenet.R

# ── Part 11: CNV inference (CopyKAT) ─────────────────────────
run_R_if_exists ./code/spokes/cnv_copykat.R

# ── Part 12: Co-expression networks (hdWGCNA) ────────────────
run_R_if_exists ./code/spokes/coexpression_hdwgcna.R

# ── Part 13: RNA velocity (STARsolo alignment, then scVelo) ──
run_if_exists    ./code/spokes/starsolo_velocity.sh
run_PY_if_exists ./code/spokes/velocity_scvelo.py

# ── Part 14: Cell fate (CellRank) ────────────────────────────
run_PY_if_exists ./code/spokes/fate_cellrank.py

# ── Part 15: Publication figures (scplotter) ─────────────────
run_R_if_exists ./code/spokes/viz_scplotter.R

# ── Part 16: TF activity (decoupleR / CollecTRI) ─────────────
run_R_if_exists ./code/spokes/grn_decoupler.R

# ── Part 17: Pathway activity (PROGENy) ──────────────────────
run_R_if_exists ./code/spokes/pathway_progeny.R

echo "========== Job finished: $(date) =========="
