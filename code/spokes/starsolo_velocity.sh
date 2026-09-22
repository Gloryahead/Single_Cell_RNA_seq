#!/usr/bin/env bash
# ---------------------------------------------------------------------
# starsolo_velocity.sh — STARsolo alignment producing spliced/unspliced
#                        counts needed by scVelo (Part 13).
# Run: bash code/spokes/starsolo_velocity.sh
# Env: scrna conda env (STAR is installed there)
# Outputs: results/velocity/<sample>/Solo.out/
# ---------------------------------------------------------------------
set -euo pipefail
DATA="${SCRNA_DATA:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"

THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-16}}"
STAR_REF="$DATA/references/star_index_GRCh38"
GTF="$DATA/references/GRCh38/Homo_sapiens.GRCh38.110.gtf"
FASTQ_ROOT="$DATA/data/fastq"
VEL_OUT="$DATA/results/velocity"
SW="$DATA/software"

mkdir -p "$VEL_OUT" "$DATA/references/GRCh38"

# ── Step 1: Build STAR genome index (once) ────────────────────────────────────
if [ ! -d "$STAR_REF" ]; then
  echo "Building STAR genome index (takes ~30 min)..."
  GENOME_FA="$DATA/references/GRCh38/GRCh38.primary_assembly.genome.fa"
  if [ ! -f "$GENOME_FA" ]; then
    wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
         -O "${GENOME_FA}.gz" && gunzip "${GENOME_FA}.gz"
    wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz \
         -O "${GTF}.gz" && gunzip "${GTF}.gz"
  fi
  mkdir -p "$STAR_REF"
  STAR --runMode genomeGenerate \
    --genomeDir       "$STAR_REF" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile     "$GTF" \
    --runThreadN      "$THREADS" \
    --genomeSAindexNbases 14
fi

# ── Step 2: Download 10x Chromium v3 barcode whitelist (once) ────────────────
WHITELIST="$SW/3M-february-2018.txt"
if [ ! -f "$WHITELIST" ]; then
  mkdir -p "$SW"
  wget https://github.com/10XGenomics/cellranger/raw/main/lib/python/cellranger/barcodes/3M-february-2018.txt.gz \
       -O "${WHITELIST}.gz" && gunzip "${WHITELIST}.gz"
fi

# ── Step 3: Run STARsolo per sample (Healthy + Post only for velocity) ────────
SAMPLES=(Healthy_1 Healthy_2 Healthy_3 Healthy_4
         Post_Patient_1 Post_Patient_2 Post_Patient_3 Post_Patient_4)

for SAMPLE in "${SAMPLES[@]}"; do
  OUT="$VEL_OUT/$SAMPLE"
  if [ -d "$OUT/Solo.out" ]; then
    echo "$SAMPLE already done, skipping."
    continue
  fi
  R1="$FASTQ_ROOT/${SAMPLE}/${SAMPLE}_S1_L001_R1_001.fastq.gz"
  R2="$FASTQ_ROOT/${SAMPLE}/${SAMPLE}_S1_L001_R2_001.fastq.gz"
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "FASTQ not found for $SAMPLE — run 01_download_fastq.sh first"; continue
  fi
  echo "Running STARsolo: $SAMPLE ..."
  mkdir -p "$OUT"
  STAR \
    --soloType            CB_UMI_Simple \
    --soloCBwhitelist     "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures        Gene Velocyto \
    --soloOutDir          "$OUT/Solo.out" \
    --soloOutFormatFeaturesGeneField3 "Gene Expression" \
    --genomeDir           "$STAR_REF" \
    --sjdbGTFfile         "$GTF" \
    --readFilesIn         "$R2" "$R1" \
    --readFilesCommand    zcat \
    --outSAMtype          BAM SortedByCoordinate \
    --outSAMattributes    NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outFileNamePrefix   "$OUT/" \
    --runThreadN          "$THREADS" \
    --limitBAMsortRAM     60000000000
  echo "STARsolo done: $SAMPLE -> $OUT/Solo.out/"
done

echo "Part 13 shell step complete. Run next: python code/spokes/velocity_scvelo.py"
