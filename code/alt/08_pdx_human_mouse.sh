#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 08_pdx_human_mouse.sh — separate human vs mouse reads for PDX samples. (Part 8)
# ALTERNATIVE entry point: use this INSTEAD of a plain human alignment when a
# sample is a patient-derived xenograft (human tumor + mouse stroma). Produces
# human-only FASTQs that then go through 02_align_cellranger.sh normally.
# Run (per PDX sample): bash code/alt/08_pdx_human_mouse.sh <ID>
# Tooling: align to a combined human+mouse reference, then classify reads.
#   Options: a combined Cell Ranger reference (human+mouse) + post-hoc split,
#            or Xenome/XenoCell for read classification.
# ---------------------------------------------------------------------
set -euo pipefail
DATA="${SCRNA_DATA:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"
ID=${1:?usage: 08_pdx_human_mouse.sh <SAMPLE_ID>}

# Strategy: build/download a COMBINED human(GRCh38)+mouse(GRCm39) reference,
# run cellranger count against it, then keep barcodes assigned to human genes.
# This dataset (GSE174609) is NOT PDX — this script is provided so the pipeline
# covers Part 8 for when you move to xenograft data.
echo "PDX split for $ID is a template — wire in your combined reference / Xenome step."
echo "Output target: $DATA/data/fastq/${ID}_human/  -> then 02_align_cellranger.sh"
