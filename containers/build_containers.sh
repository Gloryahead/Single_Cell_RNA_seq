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

PROJECT="${SCRNA_PROJECT:-/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq}"
cd "$PROJECT"

BUILD_R=true
BUILD_PY=true
case "${1:-}" in
  --only-r)      BUILD_PY=false ;;
  --only-python) BUILD_R=false  ;;
esac

if $BUILD_R; then
  echo "Building R container (scrna_r.sif) ..."
  apptainer build --fakeroot containers/scrna_r.sif containers/scrna_r.def
  echo "Done: containers/scrna_r.sif"
fi

if $BUILD_PY; then
  echo "Building Python container (scrna_python.sif) ..."
  apptainer build --fakeroot containers/scrna_python.sif containers/scrna_python.def
  echo "Done: containers/scrna_python.sif"
fi

echo "All containers built. Run Snakemake with:"
echo "  snakemake --profile profile/slurm --use-singularity -j 8 all"
