#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 01_download_fastq.sh — fetch GSE174609 runs from SRA, rename for Cell Ranger.
#                                                          (Tutorial Part 1a)
# Run : bash code/01_download_fastq.sh              # all 12
#       bash code/01_download_fastq.sh Healthy_1     # just one
#       THREADS=12 bash code/01_download_fastq.sh # override thread count
# Env : micromamba activate scrna   (sra-tools; pigz optional but faster)
# ---------------------------------------------------------------------
set -euo pipefail
# Activate haining group micromamba environments
set +eu; source ~/.bashrc; set -eu
mamba-haining
cd "${SCRNA_PROJECT:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"
mkdir -p results/logs data/fastq

# --- performance knobs ---
# node-local scratch for fasterq-dump temp (NEVER /xdisk). Falls back sensibly.
SCRATCH="${SLURM_TMPDIR:-${TMPDIR:-/tmp/$USER}}"
# threads = your SLURM allocation if set, else 8. Don't exceed reserved cores.
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
# parallel gzip if available, else plain gzip
if command -v pigz >/dev/null 2>&1; then ZIP="pigz -p $THREADS"; else ZIP="gzip"; fi
echo "scratch=$SCRATCH  threads=$THREADS  zip=${ZIP%% *}"

PAIRS=(
  Healthy_1:SRR14575500 Healthy_2:SRR14575501 Healthy_3:SRR14575502 Healthy_4:SRR14575503
  Pre_Patient_1:SRR14575504 Pre_Patient_2:SRR14575505 Pre_Patient_3:SRR14575506 Pre_Patient_4:SRR14575507
  Post_Patient_1:SRR14575508 Post_Patient_2:SRR14575509 Post_Patient_3:SRR14575510 Post_Patient_4:SRR14575511
)

download_one () {
  local ID="$1" SRR="$2"
  echo "[$ID / $SRR] downloading"
  local tmp="$SCRATCH/fqd.$SRR"; mkdir -p "$tmp"

  # 1) prefetch: download .sra ONCE to local scratch (resumable, no conversion)
  prefetch "$SRR" -O "$tmp"

  # 2) convert the LOCAL .sra with temp ALSO on local scratch -> writes FASTQ to $tmp
  fasterq-dump "$tmp/$SRR/$SRR.sra" \
    --split-files --include-technical \
    --threads "$THREADS" --temp "$tmp" \
    -O "$tmp"

  # 3) identify R1 (~28bp) vs R2 (cDNA) by READ LENGTH
  local R1="" R2="" f L
  for f in "$tmp/${SRR}"_*.fastq; do
    L=$(head -2 "$f" | tail -1 | tr -d '\n' | wc -c)
    if   [ "$L" -le 30 ]; then R1="$f"
    elif [ "$L" -ge 50 ]; then R2="$f"; fi
  done
  if [ -z "$R1" ] || [ -z "$R2" ]; then
    echo "[$ID] Could not ID R1/R2 by length — inspect manually"; rm -rf "$tmp"; return 1
  fi

  # 4) compress straight into the final /xdisk location (one write to the shared FS)
  mkdir -p "data/fastq/${ID}"
  $ZIP -c "$R1" > "data/fastq/${ID}/${ID}_S1_L001_R1_001.fastq.gz"
  $ZIP -c "$R2" > "data/fastq/${ID}/${ID}_S1_L001_R2_001.fastq.gz"

  rm -rf "$tmp"                        # wipe node-local scratch
  echo "[$ID] ready -> data/fastq/${ID}/"
}

for pair in "${PAIRS[@]}"; do
  id=${pair%%:*}; srr=${pair##*:}
  if [ "$#" -eq 0 ] || [ "$1" = "$id" ]; then
    download_one "$id" "$srr"
  fi
done
echo "All requested downloads complete."
