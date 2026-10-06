#!/bin/bash
#SBATCH --job-name=scRNA_star_index
#SBATCH --output=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq/results/logs/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=08:00:00

# Prepares the three STARsolo prerequisites that 00_setup.slurm did not cover:
#   1. 10x barcode whitelist  (3M-february-2018.txt)
#   2. Ensembl GRCh38 r110 GTF
#   3. STAR genome index
#
# Re-runnable: each step is skipped if the output already exists.
# Run once, before the Snakemake pipeline:
#   sbatch code/02b_build_star_index.slurm

set -euo pipefail

echo "========== Job started at: $(date) =========="

export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq
THREADS=${SLURM_NTASKS}

echo "Code root : $SCRNA_PROJECT"
echo "Data root : $SCRNA_DATA"
echo "Threads   : $THREADS"

# ── 1. 10x barcode whitelist ──────────────────────────────────────────────────
WL_DEST="${SCRNA_DATA}/software/3M-february-2018.txt"
if [[ -f "${WL_DEST}" ]]; then
    echo "Whitelist already present — skipping."
else
    # Cell Ranger 10.x stores it gzipped inside the installation
    WL_GZ="${SCRNA_DATA}/software/cellranger-10.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
    if [[ -f "${WL_GZ}" ]]; then
        echo "Extracting whitelist from Cell Ranger installation..."
        zcat "${WL_GZ}" > "${WL_DEST}"
    else
        echo "Not found in Cell Ranger — downloading from 10x Genomics..."
        curl -fL \
            "https://cf.10xgenomics.com/supp/cell-exp/3M-february-2018.txt.gz" \
          | zcat > "${WL_DEST}"
    fi
    echo "Whitelist ready: ${WL_DEST}"
fi

# ── 2. Ensembl GRCh38 release-110 GTF ────────────────────────────────────────
GTF_DIR="${SCRNA_DATA}/references/GRCh38"
GTF="${GTF_DIR}/Homo_sapiens.GRCh38.110.gtf"
mkdir -p "${GTF_DIR}"

if [[ -f "${GTF}" ]]; then
    echo "GTF already present — skipping download."
else
    echo "Downloading Ensembl GRCh38 release-110 GTF..."
    curl -fL \
        "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz" \
      | zcat > "${GTF}"
    echo "GTF ready: ${GTF}"
fi

# ── 3. STAR genome index ──────────────────────────────────────────────────────
IDX_DIR="${SCRNA_DATA}/references/star_index_GRCh38"
# Reuse the genome FASTA already in the Cell Ranger reference (saves ~4 GB download)
GENOME_FA="${SCRNA_DATA}/references/cellranger/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

if [[ -d "${IDX_DIR}" && -f "${IDX_DIR}/SA" ]]; then
    echo "STAR index already present — skipping."
else
    echo ""
    echo "Building STAR genome index (typically 2–4 h on 16 cores)..."
    mkdir -p "${IDX_DIR}"

    # Activate the scrna conda env (contains STAR)
    set +eu
    source ~/.bashrc
    mamba-haining
    micromamba activate scrna
    set -eu

    STAR \
        --runMode genomeGenerate \
        --genomeDir  "${IDX_DIR}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF}" \
        --runThreadN ${THREADS}

    echo "STAR index built: ${IDX_DIR}"
fi

echo ""
echo "========== STAR setup complete at: $(date) =========="
echo ""
echo "Next step:"
echo "  sbatch code/snakemake.slurm"
