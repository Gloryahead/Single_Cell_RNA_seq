#!/usr/bin/env python
# ---------------------------------------------------------------------
# velocity_scvelo.py — RNA velocity (scVelo).                     (Part 13)
# Python spoke: reads the .h5ad the hub exported (05_annotated.h5ad).
# Env: micromamba activate scrna-velocity
# PREREQ: spliced/unspliced layers. Cell Ranger alone does NOT give these —
#         run velocyto on the BAM, or align with STARsolo --soloFeatures Velocyto,
#         then merge that loom into this AnnData.
# Run: python code/spokes/velocity_scvelo.py
# ---------------------------------------------------------------------
import os, glob, yaml
import numpy as np
import scanpy as sc
import scvelo as scv

PROJECT = os.environ.get("SCRNA_PROJECT", "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")
os.chdir(PROJECT)
cfg = yaml.safe_load(open("config.yaml"))
out = cfg["project"]["outdir"]
os.makedirs(f"{out}/objects", exist_ok=True)
os.makedirs(f"{out}/plots", exist_ok=True)
scv.settings.set_figure_params("scvelo")

# Load the hub AnnData (exported by 05_annotate.R)
adata = sc.read_h5ad(f"{out}/objects/05_annotated.h5ad")

# Load and merge STARsolo Velocyto output (spliced/unspliced/ambiguous)
samples = cfg.get("samples", [])
vel_dir = "results/velocity"
adatas = []
for s in samples:
    sid = s["id"]
    sp_path = f"{vel_dir}/{sid}/Solo.out/Velocyto/filtered"
    if not os.path.isdir(sp_path):
        print(f"  STARsolo output missing for {sid}: {sp_path} — skipping")
        continue
    a = sc.read_mtx(f"{sp_path}/spliced.mtx").T
    genes = open(f"{sp_path}/features.tsv").read().splitlines()
    barcodes = [f"{sid}_{b}" for b in open(f"{sp_path}/barcodes.tsv").read().splitlines()]
    a.var_names = genes
    a.obs_names = barcodes
    u = sc.read_mtx(f"{sp_path}/unspliced.mtx").T
    u.var_names = genes; u.obs_names = barcodes
    am = sc.read_mtx(f"{sp_path}/ambiguous.mtx").T
    am.var_names = genes; am.obs_names = barcodes
    a.layers["spliced"] = a.X
    a.layers["unspliced"] = u.X
    a.layers["ambiguous"] = am.X
    adatas.append(a)

if not adatas:
    raise RuntimeError("No STARsolo outputs found. Run starsolo_velocity.sh first.")

ldata = adatas[0].concatenate(adatas[1:]) if len(adatas) > 1 else adatas[0]

# Keep only barcodes present in the hub Seurat object
common = adata.obs_names.intersection(ldata.obs_names)
adata = adata[common].copy()
ldata = ldata[common].copy()
adata.layers["spliced"] = ldata[common].layers["spliced"]
adata.layers["unspliced"] = ldata[common].layers["unspliced"]
adata.layers["ambiguous"] = ldata[common].layers["ambiguous"]

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="X_umap", color="celltype",
                                 save=f"{out}/plots/spoke_velocity.png")
adata.write(f"{out}/objects/spoke_velocity.h5ad")
print("velocity done ->", f"{out}/objects/spoke_velocity.h5ad")
