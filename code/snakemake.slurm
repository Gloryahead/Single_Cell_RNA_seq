#!/bin/bash
#SBATCH --job-name=scRNA_snakemake
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=120:00:00

# This job is a lightweight Snakemake orchestrator — it uses almost no CPU/RAM.
# It submits each pipeline rule as its own child SLURM job (up to 8 at a time)
# and monitors them to completion. Wall time is set high to cover the full pipeline.

echo "========== Snakemake started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq

echo "Code root : $SCRNA_PROJECT"
echo "Data root : $SCRNA_DATA"

cd "$SCRNA_PROJECT"

# Activate haining micromamba prefix, then the existing snakemake_env
set +eu
source ~/.bashrc
mamba-haining
micromamba activate snakemake_env
set -eu

echo "Snakemake: $(snakemake --version)"
echo "Containers:"
echo "  R   : $SCRNA_DATA/containers/scrna_r.sif"
echo "  Python: $SCRNA_DATA/containers/scrna_python.sif"

# --------------------------------------------------------------
### Launch full pipeline via SLURM executor
# Each rule runs as its own SLURM child job with per-rule resources.
# Up to 8 rules run simultaneously; failed jobs retry once automatically.
# --------------------------------------------------------------
snakemake \
  --profile profile/slurm \
  --use-singularity \
  --singularity-args "--bind $SCRNA_DATA" \
  -j 8 \
  all

echo ""
echo "========== Snakemake finished at: $(date) =========="
