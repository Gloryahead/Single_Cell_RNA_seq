#!/bin/bash
#SBATCH --job-name=scRNA_containers
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=08:00:00

echo "========== Job started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq

echo "Code root : $SCRNA_PROJECT"
echo "Data root : $SCRNA_DATA"

mkdir -p "$SCRNA_DATA/containers"

# --------------------------------------------------------------
### Build R container (Seurat, Harmony, SingleR, Monocle3, etc.)
# --------------------------------------------------------------
echo ""
echo "=========================================="
echo " Building R container"
echo " Started: $(date)"
echo "=========================================="
apptainer build --fakeroot \
  "$SCRNA_DATA/containers/scrna_r.sif" \
  "$SCRNA_PROJECT/containers/scrna_r.def"
echo " Finished R container: $(date)"

# --------------------------------------------------------------
### Build Python container (scVelo, CellRank, anndata, etc.)
# --------------------------------------------------------------
echo ""
echo "=========================================="
echo " Building Python container"
echo " Started: $(date)"
echo "=========================================="
apptainer build --fakeroot \
  "$SCRNA_DATA/containers/scrna_python.sif" \
  "$SCRNA_PROJECT/containers/scrna_python.def"
echo " Finished Python container: $(date)"

echo ""
echo "========== Containers built at: $(date) =========="
echo ""
echo "Next step:"
echo "  sbatch code/snakemake.slurm"
