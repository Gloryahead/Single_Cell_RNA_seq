#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 02_align_cellranger.sh — Cell Ranger count per sample.   (Tutorial Part 1b)
# Inputs  : data/fastq/<ID>/   Outputs : results/counts/<ID>/outs/...
# Run     : bash code/02_align_cellranger.sh          # all 12, sequentially
#           bash code/02_align_cellranger.sh Healthy_1  # just one (optional arg)
# Setup   : module load cellranger/10.0.0   (NOT a conda package)
# Compute : ~hours + tens of GB per run.
# ---------------------------------------------------------------------
set -euo pipefail
# Activate haining group micromamba environments
set +eu; source ~/.bashrc; set -eu
mamba-haining
cd "${SCRNA_PROJECT:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"
mkdir -p results/logs results/counts

REF=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/references/cellranger/refdata-gex-GRCh38-2024-A
FASTQ_ROOT=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/data/fastq
CORES="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-16}}"

IDS=(Healthy_1 Healthy_2 Healthy_3 Healthy_4
     Pre_Patient_1 Pre_Patient_2 Pre_Patient_3 Pre_Patient_4
     Post_Patient_1 Post_Patient_2 Post_Patient_3 Post_Patient_4)

align_one () {
  local ID="$1"
  echo "[$ID] aligning on $(hostname)"
  # cellranger writes its output dir into the CWD, so run from results/counts.
  # It refuses to overwrite an existing --id dir; remove a stale partial first.
  ( cd results/counts
    rm -rf "$ID"
    cellranger count \
      --id="$ID" \
      --transcriptome="$REF" \
      --fastqs="$FASTQ_ROOT/$ID" \
      --sample="$ID" \
      --localcores="$CORES" --localmem=90 \
      --create-bam=true
  )
  echo "[$ID] done -> results/counts/$ID/outs/filtered_feature_bc_matrix/"
}

# ---- run: all samples, or just the one passed as $1 ----
for ID in "${IDS[@]}"; do
  if [ "$#" -eq 0 ] || [ "$1" = "$ID" ]; then
    align_one "$ID"
  fi
done
echo "All requested alignments complete."
