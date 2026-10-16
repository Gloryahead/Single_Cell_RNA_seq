#!/bin/bash
#SBATCH --job-name=setup_networks
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=00:30:00

# Pre-downloads CollecTRI and PROGENy network CSVs to results/references/.
# Run this ONCE before submitting the main pipeline if you want the downloads
# to happen in a controlled job rather than during the Snakemake orchestrator.
#
#   sbatch code/setup_networks.slurm
#
# The Snakemake pipeline (snakemake.slurm) includes a localrule that does
# the same thing automatically — this script is an explicit alternative.

echo "========== setup_networks started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq

apptainer exec \
  --bind "$SCRNA_DATA" \
  "$SCRNA_DATA/containers/scrna_r.sif" \
  Rscript "$SCRNA_PROJECT/code/setup_networks.R"

echo "========== setup_networks finished at: $(date) =========="
