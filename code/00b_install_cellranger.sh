#!/usr/bin/env bash
# ---------------------------------------------------------------------
# 00b_install_cellranger.sh — install Cell Ranger into the project's software/
# and expose it on PATH. Runs from ANY directory (anchors to $SCRNA_PROJECT).
#
# Cell Ranger is NOT a conda package (10x license forbids redistribution), so
# it lives beside the scrna env and is prepended to PATH via a sourced snippet.
#
# ONE-TIME SETUP (interactive session on a compute node — NOT the login node):
#   export SCRNA_PROJECT=/home/u11/maarowosegbe/Single_Cell_RNA_seq
#   export SCRNA_DATA=/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq
#   bash $SCRNA_PROJECT/code/00b_install_cellranger.sh
# NOTE: the signed CR_URL below EXPIRES. If you see "Access Denied", regenerate
#       it from https://www.10xgenomics.com/support/software/cell-ranger/downloads
#       and paste the new one into CR_URL.
# ---------------------------------------------------------------------
set -euo pipefail

PROJECT="${SCRNA_PROJECT:-/home/u11/maarowosegbe/Single_Cell_RNA_seq}"
DATA="${SCRNA_DATA:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"

VERSION="10.1.0"
SW="$DATA/software"
DEST="$SW/cellranger-${VERSION}"
mkdir -p "$SW"

# signed URL from the 10x download page (expires — refresh if it fails):
CR_URL="https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.1.0.tar.gz?Expires=1790164870&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=DpPd91mIuiNfKqsav9~TKCFxjAH~bNeIA-Ojg6XvMsuhpZlO99TUTm~OK-Ze~BvnmhlY4rWCbYySR72NAOpBnSNRuFAxNUWdV2txoX8pUkYzF7XDHoLBn3kdb4XNTXfs6t0iUCAxzM2uxyZe34z3Lb~PPMn3zMCHrpBDKp4116ErKlodw-xY34Z~8je5x8dlng9L~Kkb9pJsJpx4Zh~7gWGz4OKRPx7mqvwbWt04mm83xOgqGp3AEPDN78Te-1xegusrDVUYradejngev-vafRn3KAgtLTR~DLDcGMe9WT5qH10Vr8iDhXSbW750aT06XfWbaVU1F~rEqkAsQbhOdQ__"

if [[ -x "$DEST/cellranger" ]]; then
  echo "Cell Ranger already installed at $DEST"
else
  [[ -z "$CR_URL" ]] && { echo "ERROR: CR_URL empty — paste a fresh signed URL from 10x."; exit 1; }
  echo "Downloading Cell Ranger ${VERSION} ..."
  curl -fL "$CR_URL" -o "$SW/cellranger-${VERSION}.tar.gz" || {
    echo "Download failed. If HTTP 403/AccessDenied, the signed URL EXPIRED — get a fresh one."; exit 1; }
  echo "Unpacking ..."
  tar -xzvf "$SW/cellranger-${VERSION}.tar.gz" -C "$SW"
  rm -f "$SW/cellranger-${VERSION}.tar.gz"
fi

# sanity check
"$DEST/cellranger" --version

# activation snippet to source from the SLURM job
cat > "$SW/activate_cellranger.sh" <<ACT
# source this to put Cell Ranger on PATH:
export PATH="${DEST}:\$PATH"
ACT
echo
echo "Installed. In your SLURM job, AFTER 'micromamba activate scrna', add:"
echo "    source \"\${SCRNA_DATA:-$DATA}/software/activate_cellranger.sh\""
