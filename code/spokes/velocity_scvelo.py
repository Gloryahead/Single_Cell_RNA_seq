#!/usr/bin/env python
# ---------------------------------------------------------------------
# velocity_scvelo.py — RNA velocity (scVelo).                     (Part 13)
#
# Strategy:
#   R exports metadata + UMAP as CSVs (05_annotate.R).
#   AnnData built from STARsolo Velocyto matrices; CSVs joined to obs.
#   Stochastic mode: recover_dynamics() OOMs at 40 GB on 32K PBMC cells.
#   Stochastic estimates velocity from 2nd-order moments — valid, widely used.
#   Dynamical-mode-only analyses (latent_time, phase portraits by fit_likelihood,
#   differential_kinetic_test) are omitted.
#
# In : results/velocity/{sample}/Solo.out/Velocyto/filtered/
#      results/objects/05_annotated_meta.csv
#      results/objects/05_annotated_umap.csv
# Out: results/objects/spoke_velocity.h5ad
#      results/plots/spoke_velocity_*.png
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
n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
os.makedirs(f"{out}/objects", exist_ok=True)
os.makedirs(f"{out}/plots", exist_ok=True)

scv.settings.set_figure_params("scvelo")
scv.settings.figdir     = f"{out}/plots"
scv.settings.n_jobs     = n_cores
sc.settings.n_jobs      = n_cores

# ---- Load Seurat metadata and UMAP ----------------------------------------
meta_df = pd.read_csv(f"{out}/objects/05_annotated_meta.csv", index_col="barcode")
umap_df = pd.read_csv(f"{out}/objects/05_annotated_umap.csv", index_col="barcode")

# ---- Load STARsolo Velocyto outputs ---------------------------------------
samples = cfg.get("samples", [])
VEL_COND = {"Healthy", "Post"}
vel_dir  = f"{out}/velocity"

ldata_list = []
for s in samples:
    if s.get("condition") not in VEL_COND:
        continue
    sid    = s["id"]
    sp_path = f"{vel_dir}/{sid}/Solo.out/Velocyto/filtered"
    if not os.path.isdir(sp_path):
        print(f"  STARsolo Velocyto output missing for {sid}: {sp_path} — skipping")
        continue

    # features.tsv: tab-delimited GeneID, GeneName[, Type]
    feat_df  = pd.read_csv(f"{sp_path}/features.tsv", sep="\t", header=None)
    gene_ids = feat_df[0].values
    gene_sym = feat_df[1].values if feat_df.shape[1] > 1 else gene_ids

    # STARsolo barcodes have no -1 suffix; Seurat merged barcodes = "{sid}_{bc}-1"
    barcodes = [f"{sid}_{b}-1"
                for b in open(f"{sp_path}/barcodes.tsv").read().splitlines()]

    sp = csr_matrix(mmread(f"{sp_path}/spliced.mtx").T)
    un = csr_matrix(mmread(f"{sp_path}/unspliced.mtx").T)
    am = csr_matrix(mmread(f"{sp_path}/ambiguous.mtx").T)

    a = ad.AnnData(X=sp)
    a.var_names          = pd.Index(gene_ids)
    a.var["gene_symbols"] = gene_sym
    a.obs_names          = barcodes
    a.layers["spliced"]   = sp
    a.layers["unspliced"] = un
    a.layers["ambiguous"] = am
    ldata_list.append(a)

if not ldata_list:
    raise RuntimeError("No STARsolo Velocyto outputs found. Run starsolo first.")

# Concatenate across samples; inner join keeps only genes present in all samples
ldata = ad.concat(ldata_list, join="inner") if len(ldata_list) > 1 else ldata_list[0]

# Switch var_names to gene symbols (Ensembl IDs saved in var["gene_ids"])
ldata.var["gene_ids"] = ldata.var_names.values
ldata.var_names       = pd.Index(ldata.var["gene_symbols"].values)
ldata.var_names_make_unique()

# ---- Intersect with Seurat metadata ---------------------------------------
shared = meta_df.index.intersection(ldata.obs_names)
if len(shared) == 0:
    raise RuntimeError("No shared barcodes between STARsolo output and Seurat metadata.")
print(f"Shared barcodes: {len(shared)}")

adata = ldata[shared].copy()
adata.obs = meta_df.loc[shared].copy()
adata.obsm["X_umap"] = umap_df.loc[shared, ["UMAP_1", "UMAP_2"]].values.astype(float)

# ---- Preprocessing --------------------------------------------------------
# Step 1: filter genes with insufficient shared (spliced+unspliced) counts
scv.pp.filter_genes(adata, min_shared_counts=20)

# Steps 2-4: normalize, log-transform, select top HVGs
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=True)
print(f"After HVG selection: {adata.n_vars} genes, {adata.n_obs} cells")

# Step 5: compute moments (also calls PCA + neighbors internally)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# ---- RNA velocity (stochastic mode) ---------------------------------------
# Stochastic: velocity estimated from 2nd-order moments of splicing kinetics.
# Dynamical mode calls recover_dynamics() which OOMed at 40 GB on 32K PBMC cells.
velocity_ok = False
try:
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata, n_jobs=n_cores)
    velocity_ok = True
    adata.write(f"{out}/objects/spoke_velocity.h5ad")
    print("Checkpoint saved -> spoke_velocity.h5ad")
except Exception as e:
    print(f"WARNING: velocity computation failed ({type(e).__name__}: {e})")
    print("Saving AnnData without velocity layers.")
    adata.write(f"{out}/objects/spoke_velocity.h5ad")

if not velocity_ok:
    print("velocity done (no plots) ->", f"{out}/objects/spoke_velocity.h5ad")
    raise SystemExit(0)

# ---- Velocity confidence --------------------------------------------------
try:
    scv.tl.velocity_confidence(adata)
    scv.pl.scatter(adata, c=["velocity_length", "velocity_confidence"],
                   cmap="coolwarm", perc=[5, 95],
                   save=f"{out}/plots/spoke_velocity_confidence.png",
                   show=False)
    print("Plot: velocity confidence")
except Exception as e:
    print(f"[skip] velocity confidence: {e}")

# ---- Combined stream plot (all conditions) --------------------------------
try:
    scv.pl.velocity_embedding_stream(
        adata, basis="X_umap", color="celltype",
        title="RNA velocity — all conditions",
        save=f"{out}/plots/spoke_velocity_stream.png",
        show=False)
    print("Plot: velocity stream (all conditions)")
except Exception as e:
    print(f"[skip] stream plot: {e}")

# ---- Arrow embedding plot -------------------------------------------------
try:
    scv.pl.velocity_embedding(
        adata, basis="X_umap", color="celltype",
        arrow_length=3, arrow_size=2,
        save=f"{out}/plots/spoke_velocity_arrows.png",
        show=False)
    print("Plot: velocity arrows")
except Exception as e:
    print(f"[skip] arrow plot: {e}")

# ---- Per-condition stream plots -------------------------------------------
for cond in sorted(adata.obs["condition"].unique()):
    try:
        sub = adata[adata.obs["condition"] == cond].copy()
        scv.pl.velocity_embedding_stream(
            sub, basis="X_umap", color="celltype",
            title=f"RNA velocity — {cond}",
            save=f"{out}/plots/spoke_velocity_stream_{cond}.png",
            show=False)
        print(f"Plot: velocity stream ({cond})")
    except Exception as e:
        print(f"[skip] stream plot {cond}: {e}")

# ---- Phase portrait scatter (top velocity genes) --------------------------
try:
    # Rank by R² of the steady-state fit; fall back to any velocity gene
    vg = adata.var[adata.var.get("velocity_genes", pd.Series(True, index=adata.var_names))]
    if "velocity_r2" in vg.columns:
        top_genes = vg.nlargest(10, "velocity_r2").index.tolist()
    else:
        top_genes = vg.index[:10].tolist()
    if len(top_genes) >= 4:
        scv.pl.scatter(adata, var_names=top_genes, ncols=5,
                       save=f"{out}/plots/spoke_velocity_phase_portraits.png",
                       show=False)
        print(f"Plot: phase portraits ({len(top_genes)} genes)")
except Exception as e:
    print(f"[skip] phase portrait: {e}")

# ---- Top velocity genes ranked list (CSV) ---------------------------------
try:
    scv.tl.rank_velocity_genes(adata, groupby="celltype", min_corr=0.3)
    vg_df = scv.get_df(adata, "rank_velocity_genes/names")
    vg_df.to_csv(f"{out}/tables/spoke_velocity_top_genes.csv")
    print(f"Table: top velocity genes per cell type")
except Exception as e:
    print(f"[skip] rank_velocity_genes: {e}")

# ---- Final save -----------------------------------------------------------
adata.write(f"{out}/objects/spoke_velocity.h5ad")
print("velocity done ->", f"{out}/objects/spoke_velocity.h5ad")
