#!/usr/bin/env python
# ---------------------------------------------------------------------
# velocity_scvelo.py — RNA velocity (scVelo).                     (Part 13)
# Tutorial approach (ngs101.com Part 13): R exports metadata + UMAP as
# CSVs; this script builds AnnData from STARsolo Velocyto matrices and
# joins those CSVs. No h5ad from R needed — eliminates the sceasy OOM.
# Env: micromamba activate scrna-velocity
# PREREQ: STARsolo --soloFeatures Velocyto outputs
# Run: python code/spokes/velocity_scvelo.py
# ---------------------------------------------------------------------
import os, yaml
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvelo as scv
from scipy.sparse import csr_matrix
from scipy.io import mmread

PROJECT = os.environ.get("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
os.chdir(PROJECT)
cfg = yaml.safe_load(open("config.yaml"))
out = cfg["project"]["outdir"]
os.makedirs(f"{out}/objects", exist_ok=True)
os.makedirs(f"{out}/plots", exist_ok=True)
scv.settings.set_figure_params("scvelo")

# Load Seurat metadata and UMAP exported by 05_annotate.R
meta_df = pd.read_csv(f"{out}/objects/05_annotated_meta.csv", index_col="barcode")
umap_df = pd.read_csv(f"{out}/objects/05_annotated_umap.csv", index_col="barcode")

# Load STARsolo Velocyto outputs (spliced / unspliced / ambiguous)
samples = cfg.get("samples", [])
VEL_COND = {"Healthy", "Post"}
vel_dir = f"{out}/velocity"
ldata_list = []
for s in samples:
    if s.get("condition") not in VEL_COND:
        continue
    sid = s["id"]
    sp_path = f"{vel_dir}/{sid}/Solo.out/Velocyto/filtered"
    if not os.path.isdir(sp_path):
        print(f"  STARsolo Velocyto output missing for {sid}: {sp_path} — skipping")
        continue
    genes    = open(f"{sp_path}/features.tsv").read().splitlines()
    barcodes = [f"{sid}_{b}" for b in open(f"{sp_path}/barcodes.tsv").read().splitlines()]
    sp = csr_matrix(mmread(f"{sp_path}/spliced.mtx").T)
    un = csr_matrix(mmread(f"{sp_path}/unspliced.mtx").T)
    am = csr_matrix(mmread(f"{sp_path}/ambiguous.mtx").T)
    a = ad.AnnData(X=sp)
    a.var_names = genes
    a.obs_names = barcodes
    a.layers["spliced"]   = sp
    a.layers["unspliced"] = un
    a.layers["ambiguous"] = am
    ldata_list.append(a)

if not ldata_list:
    raise RuntimeError("No STARsolo Velocyto outputs found. Run starsolo first.")

ldata = ad.concat(ldata_list, join="inner") if len(ldata_list) > 1 else ldata_list[0]

# Keep only barcodes present in the Seurat object (survived QC + annotation)
shared = meta_df.index.intersection(ldata.obs_names)
if len(shared) == 0:
    raise RuntimeError("No shared barcodes between STARsolo output and Seurat metadata.")
print(f"Shared barcodes: {len(shared)}")

adata = ldata[shared].copy()
adata.obs = meta_df.loc[shared]
adata.obsm["X_umap"] = umap_df.loc[shared, ["UMAP_1", "UMAP_2"]].values

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="X_umap", color="celltype",
                                 save=f"{out}/plots/spoke_velocity.png")
adata.write(f"{out}/objects/spoke_velocity.h5ad")
print("velocity done ->", f"{out}/objects/spoke_velocity.h5ad")
