#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 01b_fastqc.sh — read-level QC on raw FASTQ (FastQC + MultiQC). (Part 1)
# Runs AFTER 01_download_fastq.sh, BEFORE 02_align_cellranger.sh.
# Inputs  : data/fastq/<ID>/*_R[12]_001.fastq.gz
# Outputs : results/qc/fastqc/*.html  +  results/qc/multiqc_report.html
# Run     : bash code/01b_fastqc.sh              # all samples
# Env     : micromamba activate scrna   (fastqc, multiqc)
# ---------------------------------------------------------------------
set -euo pipefail
# Activate haining group micromamba environments
set +eu; source ~/.bashrc; set -eu
mamba-haining
DATA="${SCRNA_DATA:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-6}}"
mkdir -p "$DATA/results/qc/fastqc"

# gather every FASTQ produced by the download step
mapfile -t FQS < <(find "$DATA/data/fastq" -name '*_001.fastq.gz' | sort)
if [ "${#FQS[@]}" -eq 0 ]; then
  echo "No FASTQ found under $DATA/data/fastq — run 01_download_fastq.sh first"; exit 1
fi
echo "Running FastQC on ${#FQS[@]} files with $THREADS threads"

# 1) per-file read-quality reports
fastqc -t "$THREADS" -o "$DATA/results/qc/fastqc" "${FQS[@]}"

# 2) aggregate all reports into ONE summary
multiqc "$DATA/results/qc/fastqc" -o "$DATA/results/qc" -n multiqc_report --force

echo "Done -> open results/qc/multiqc_report.html"
echo "For 10x: check R1 length (~28bp barcode+UMI) and R2 per-base quality."
