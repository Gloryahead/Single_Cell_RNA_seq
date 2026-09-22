#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 00_setup_dirs.sh — create the project skeleton (idempotent; safe to rerun).
# Project : GSE174609 PBMC scRNA-seq (periodontitis ctrl/pre/post treatment)
# Run     : bash code/00_setup_dirs.sh
# ---------------------------------------------------------------------
set -euo pipefail
DATA="${SCRNA_DATA:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"

mkdir -p "$DATA/data/fastq"
mkdir -p "$DATA/references/cellranger"
mkdir -p "$DATA/references/GRCh38"
mkdir -p "$DATA/references/star_index_GRCh38"
mkdir -p "$DATA/software"
mkdir -p "$DATA/containers"
mkdir -p "$DATA/results/counts"
mkdir -p "$DATA/results/qc"
mkdir -p "$DATA/results/objects"
mkdir -p "$DATA/results/tables"
mkdir -p "$DATA/results/plots"
mkdir -p "$DATA/results/logs"
mkdir -p "$DATA/results/velocity"

echo "Data skeleton ready under $DATA"
find "$DATA" -maxdepth 2 -type d | sort
