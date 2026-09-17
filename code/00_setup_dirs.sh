#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 00_setup_dirs.sh — create the project skeleton (idempotent; safe to rerun).
# Project : GSE174609 PBMC scRNA-seq (periodontitis ctrl/pre/post treatment)
# Run     : bash code/00_setup_dirs.sh
# ---------------------------------------------------------------------
set -euo pipefail
PROJECT=${SCRNA_PROJECT:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}
cd "$PROJECT"

# SCRNA_PROJECT = git repo root AND xdisk data root (clone the repo to xdisk)
mkdir -p data/fastq
mkdir -p references/cellranger
mkdir -p references/GRCh38
mkdir -p references/star_index_GRCh38
mkdir -p software
mkdir -p results/counts
mkdir -p results/qc
mkdir -p results/objects
mkdir -p results/tables
mkdir -p results/plots
mkdir -p results/logs
mkdir -p results/velocity

echo "Project skeleton ready under $PROJECT"
find "$PROJECT" -maxdepth 2 -type d | sort
