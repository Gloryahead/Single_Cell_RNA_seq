#!/bin/bash
#SBATCH --job-name=scRNA_align
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=24:00:00

echo "========== Job started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq

echo "Code root : $SCRNA_PROJECT"
echo "Data root : $SCRNA_DATA"
echo "Cores     : $SLURM_NTASKS"

# Activate haining micromamba prefix, then put Cell Ranger on PATH
set +eu
source ~/.bashrc
mamba-haining
source "${SCRNA_DATA}/software/activate_cellranger.sh"
set -eu

echo "Cell Ranger: $(cellranger --version 2>&1 | head -1)"

# --------------------------------------------------------------
### Align all 12 samples with Cell Ranger count (sequential)
# Healthy×4, Pre×4, Post×4 → results/counts/<ID>/outs/
# Cell Ranger uses --localcores=$SLURM_NTASKS and --localmem=90
# (~1-2 h per sample; 12 samples fit in 24 h wall time)
# --------------------------------------------------------------
echo ""
echo "=========================================="
echo " Running: 02_align_cellranger.sh"
echo " Started: $(date)"
echo "=========================================="
bash "${SCRNA_PROJECT}/code/02_align_cellranger.sh"
echo " Finished: $(date)"

echo ""
echo "========== Alignment complete at: $(date) =========="
echo ""
echo "Next steps:"
echo "  Check outputs : ls \$SCRNA_DATA/results/counts/"
echo "  Run pipeline  : sbatch code/03_pipeline.slurm   (or snakemake)"
