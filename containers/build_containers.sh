#!/usr/bin/env bash
# ---------------------------------------------------------------------
# build_containers.sh — build Apptainer SIF images from .def files.
#
# Run on a compute node (NOT the login node — needs several GB RAM):
#   bash containers/build_containers.sh
#
# The resulting .sif files are referenced by the Snakemake workflow
# via the `singularity:` directive in each rule.
# Flags:
#   --only-r        build only the R container
#   --only-python   build only the Python container
# ---------------------------------------------------------------------
set -euo pipefail

CODE="${SCRNA_PROJECT:-/home/u11/maarowosegbe/Single_Cell_RNA_seq}"
DATA="${SCRNA_DATA:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"
mkdir -p "$DATA/containers"

BUILD_R=true
BUILD_PY=true
case "${1:-}" in
  --only-r)      BUILD_PY=false ;;
  --only-python) BUILD_R=false  ;;
esac

if $BUILD_R; then
  echo "Building R container -> $DATA/containers/scrna_r.sif ..."
  apptainer build --fakeroot "$DATA/containers/scrna_r.sif" "$CODE/containers/scrna_r.def"
  echo "Done: $DATA/containers/scrna_r.sif"
fi

if $BUILD_PY; then
  echo "Building Python container -> $DATA/containers/scrna_python.sif ..."
  apptainer build --fakeroot "$DATA/containers/scrna_python.sif" "$CODE/containers/scrna_python.def"
  echo "Done: $DATA/containers/scrna_python.sif"
fi

echo "All containers built. Run Snakemake with:"
echo "  snakemake --profile profile/slurm --use-singularity -j 8 all"
