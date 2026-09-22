#!/bin/bash
#SBATCH --job-name=scRNA_setup
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=4:00:00

echo "========== Job started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq
export THREADS=${SLURM_NTASKS}

# --output path above requires results/logs/ to already exist.
# Create it manually before the first sbatch:
#   mkdir -p /xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs
echo "Code root : $SCRNA_PROJECT"
echo "Data root : $SCRNA_DATA"

run_script() {
    local SCRIPT=$1
    if [[ -f "${SCRIPT}" ]]; then
        echo ""
        echo "=========================================="
        echo " Running: $(basename ${SCRIPT})"
        echo " Started: $(date)"
        echo "=========================================="
        bash "${SCRIPT}"
        echo " Finished: $(date)"
    else
        echo "WARNING: Script not found: ${SCRIPT}"
    fi
}

# --------------------------------------------------------------
### Step 1 — Create xdisk directory skeleton (idempotent)
# --------------------------------------------------------------
run_script ${SCRNA_PROJECT}/code/00_setup_dirs.sh

# --------------------------------------------------------------
### Step 2 — Install Cell Ranger into $SCRNA_DATA/software/
# NOTE: The signed CR_URL in 00b_install_cellranger.sh EXPIRES.
#       Regenerate it from https://www.10xgenomics.com/support/software/cell-ranger/downloads
#       before running this step if you see an HTTP 403 error.
# --------------------------------------------------------------
run_script ${SCRNA_PROJECT}/code/00b_install_cellranger.sh

# --------------------------------------------------------------
### Step 3 — Download Cell Ranger GRCh38 reference genome (~20 GB)
# --------------------------------------------------------------
echo ""
echo "=========================================="
echo " Downloading Cell Ranger GRCh38 reference"
echo " Started: $(date)"
echo "=========================================="
REF_DIR="${SCRNA_DATA}/references/cellranger"
mkdir -p "${REF_DIR}"
TARBALL="${REF_DIR}/refdata-gex-GRCh38-2024-A.tar.gz"

if [[ -d "${REF_DIR}/refdata-gex-GRCh38-2024-A" ]]; then
    echo "Reference already present — skipping download."
else
    curl -fL \
      "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" \
      -o "${TARBALL}"
    echo "Unpacking ..."
    tar -xzf "${TARBALL}" -C "${REF_DIR}"
    rm -f "${TARBALL}"
    echo "Reference ready: ${REF_DIR}/refdata-gex-GRCh38-2024-A"
fi
echo " Finished: $(date)"

echo ""
echo "========== Setup complete at: $(date) =========="
echo ""
echo "Next steps:"
echo "  1. Create conda envs (if not done):  see evns/csRNA.yml + evns/scRNA_python.yml"
echo "  2. Download FASTQs:  sbatch code/01_download.slurm   (or submit 01_download_fastq.sh)"
echo "  3. Align:            sbatch code/02_align.slurm      (after FASTQs finish)"
echo "  4. Analyse:          snakemake --profile profile/slurm -j 8 all"
