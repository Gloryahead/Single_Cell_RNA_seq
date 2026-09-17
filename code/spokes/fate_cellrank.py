#!/usr/bin/env python
# ---------------------------------------------------------------------
# fate_cellrank.py — cell-fate probabilities (CellRank).          (Part 14)
# Python spoke: continues from the velocity output (spoke_velocity.h5ad).
# Env: micromamba activate scrna-velocity
# Run: python code/spokes/fate_cellrank.py
# HONESTY NOTE: CellRank needs a real directional process. The NGS101 series
#   found it WASN'T appropriate for static PBMCs — use this on developmental
#   data (or skip for PBMC). Kept here for completeness of the 17-part set.
# ---------------------------------------------------------------------
import os, yaml, scanpy as sc, cellrank as cr

PROJECT = os.environ.get("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
os.chdir(PROJECT)
cfg = yaml.safe_load(open("config.yaml"))
out = cfg["project"]["outdir"]

adata = sc.read_h5ad(f"{out}/objects/spoke_velocity.h5ad")

vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
combined = 0.8 * vk + 0.2 * ck

g = cr.estimators.GPCCA(combined)
g.fit(cluster_key="celltype", n_states=None)
g.compute_fate_probabilities()
g.plot_fate_probabilities(save=f"{out}/plots/spoke_cellrank_fate.png")
adata.write(f"{out}/objects/spoke_cellrank.h5ad")
print("cellrank done ->", f"{out}/objects/spoke_cellrank.h5ad")
