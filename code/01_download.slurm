#!/bin/bash
#SBATCH --job-name=scRNA_download
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00

echo "========== Job started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq
export THREADS=${SLURM_NTASKS}

echo "Code root : $SCRNA_PROJECT"
echo "Data root : $SCRNA_DATA"
echo "Threads   : $THREADS"

# Activate haining micromamba prefix, then the scrna env (provides sra-tools)
set +eu
source ~/.bashrc
mamba-haining
micromamba activate scrna
set -eu

echo "sra-tools: $(prefetch --version 2>&1 | head -1)"

# --------------------------------------------------------------
### Download all 12 GSE174609 SRA runs and rename for Cell Ranger
# SRR14575500-11: Healthy×4, Pre×4, Post×4
# Output: $SCRNA_DATA/data/fastq/{SampleID}/{ID}_S1_L001_R[12]_001.fastq.gz
# --------------------------------------------------------------
echo ""
echo "=========================================="
echo " Running: 01_download_fastq.sh"
echo " Started: $(date)"
echo "=========================================="
bash "${SCRNA_PROJECT}/code/01_download_fastq.sh"
echo " Finished: $(date)"

echo ""
echo "========== Download complete at: $(date) =========="
echo ""
echo "Next steps:"
echo "  Verify FASTQs : ls \$SCRNA_DATA/data/fastq/"
echo "  Align         : sbatch code/02_align.slurm"
