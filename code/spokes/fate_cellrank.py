#!/usr/bin/env python
# ---------------------------------------------------------------------
# fate_cellrank.py — cell-fate probabilities (CellRank 2).         (Part 14)
#
# BIOLOGICAL NOTE: CellRank computes absorption probabilities on a velocity-
# derived Markov chain. This is meaningful ONLY for developing/differentiating
# tissues where velocity arrows encode a genuine directional process. The
# GSE174609 PBMC dataset contains already-differentiated mature cells. The
# NGS101 tutorial (Part 14) explicitly found that CellRank assigned 85-92%
# B-cell fate to every cell type in this dataset. Results here are kept for
# completeness of the 17-part set. Interpret with caution.
#
# HPC NOTE: The default GPCCA Schur decomposition uses SLEPc (PETSc), which
# calls MPI init. On Puma, MPICH OFI MPI init Abort()s inside Singularity
# because no InfiniBand fabric is present. Abort() bypasses Python except.
# Fix: conda install petsc slepc with --without-mpi, OR set env vars below,
# OR request an MPI-capable node via --constraint=ib. Rule is commented out
# of rule all until the container issue is resolved.
#
# In : results/objects/spoke_velocity.h5ad
# Out: results/objects/spoke_cellrank.h5ad
#      results/tables/spoke_cellrank_drivers.csv
#      results/plots/spoke_cellrank_*.png
# ---------------------------------------------------------------------
import os, sys, warnings, logging, yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv
import cellrank as cr

# Attempt to avoid MPICH OFI Abort in Singularity (no InfiniBand fabric)
os.environ.setdefault("PETSC_OPTIONS",      "-mpi_fake")
os.environ.setdefault("OMPI_MCA_btl",       "^openib")
os.environ.setdefault("FI_PROVIDER",         "sockets")

warnings.filterwarnings("ignore", category=FutureWarning)
logging.getLogger("cellrank").setLevel(logging.ERROR)
sc.settings.verbosity  = 1
scv.settings.verbosity = 1

PROJECT = os.environ.get("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
os.chdir(PROJECT)
cfg = yaml.safe_load(open("config.yaml"))
out = cfg["project"]["outdir"]
os.makedirs(f"{out}/objects", exist_ok=True)
os.makedirs(f"{out}/plots",   exist_ok=True)
os.makedirs(f"{out}/tables",  exist_ok=True)

scv.settings.figdir = f"{out}/plots"
sc.settings.figdir  = f"{out}/plots"

# ---- Step 1: Load velocity AnnData ----------------------------------------
adata = sc.read_h5ad(f"{out}/objects/spoke_velocity.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# ---- Step 2: Visualize velocity field before handing to CellRank ----------
if "velocity_graph" in adata.uns:
    try:
        # velocity_graph already computed by scVelo; just re-plot
        scv.pl.velocity_embedding_stream(
            adata, basis="X_umap", color="celltype",
            title="RNA Velocity Stream (pre-CellRank check)",
            save=f"{out}/plots/spoke_cellrank_velocity_check.png", show=False)
        print("Plot: velocity stream check")
    except Exception as e:
        print(f"[skip] velocity stream check: {e}")

# ---- Step 3: Build the combined transition matrix -------------------------
# VelocityKernel reads from adata.uns["velocity_graph"] (set by scv.tl.velocity_graph)
# ConnectivityKernel reads from adata.obsp["connectivities"] (set by sc.pp.neighbors)
has_velocity = "velocity_graph" in adata.uns and "velocity" in adata.layers

if has_velocity:
    print("Building 80% VelocityKernel + 20% ConnectivityKernel")
    vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    combined = 0.8 * vk + 0.2 * ck
else:
    print("WARNING: No velocity data found — using ConnectivityKernel only.")
    print("Fate probabilities will reflect kNN topology, not developmental direction.")
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    combined = ck

print(f"Transition matrix shape: {combined.transition_matrix.shape}")

# ---- Step 4: GPCCA macrostates --------------------------------------------
g = cr.estimators.GPCCA(combined)

n_celltypes = adata.obs["celltype"].nunique()
n_states    = min(15, max(5, n_celltypes * 2))
print(f"\nCell types detected: {n_celltypes}  →  n_states = {n_states}")
print("Cell type counts:")
print(adata.obs["celltype"].value_counts().to_string())

try:
    # Step 4.1-4.2: Schur decomposition
    # n_components > n_states (extra components aid stability)
    n_components = min(n_states + 5, adata.n_obs - 1)
    g.compute_schur(n_components=n_components)
    print("Schur decomposition complete.")

    # Step 4.3: compute macrostates
    g.compute_macrostates(n_states=n_states, cluster_key="celltype")

    # Step 4.4: coarse-T stability check
    stability = pd.Series(
        np.diag(g.coarse_T.values),
        index=g.coarse_T.index
    ).sort_values(ascending=False)
    print("\nMacrostate stability (diagonal of coarse_T):")
    print(stability.round(3).to_string())

    off_diag = g.coarse_T.values - np.diag(np.diag(g.coarse_T.values))
    print(f"\nLargest negative off-diagonal: {off_diag.min():.4f}")

    # Step 4.5: macrostate UMAP
    try:
        g.plot_macrostates(which="all", basis="X_umap",
                           legend_loc="right margin",
                           title=f"Macrostates (k={n_states})",
                           save=f"{out}/plots/spoke_cellrank_macrostates.png",
                           show=False)
        print("Plot: macrostates UMAP")
    except Exception as e:
        print(f"[skip] macrostates plot: {e}")

    # Verify macrostate → cell type mapping
    if "macrostates_fwd" in adata.obs.columns:
        ms_df = adata.obs[["celltype", "macrostates_fwd"]].dropna()
        print("\nMacrostate — cell type crosstab (top 20 rows):")
        print(pd.crosstab(ms_df["macrostates_fwd"], ms_df["celltype"]).head(20))

except Exception as e:
    print(f"FATAL: Schur/macrostate computation failed: {e}")
    print("This may be the Singularity MPICH/OFI Abort(). See header note.")
    # Save minimal output and exit gracefully
    adata.write(f"{out}/objects/spoke_cellrank.h5ad")
    sys.exit(0)

# ---- Step 5: Terminal and initial states -----------------------------------
try:
    # Step 5.1: automatic stability filter
    g.predict_terminal_states()
    if g.terminal_states is not None:
        print("\nTerminal states (automatic filter):")
        print(g.terminal_states.value_counts().to_string())
    else:
        print("WARNING: No macrostates cleared the default stability threshold.")
        # Fallback: use top-3 most stable macrostates
        top3 = stability.index[:3].tolist()
        print(f"Fallback: manually setting top-3 stable macrostates as terminal: {top3}")
        g.set_terminal_states(states=top3)

    # Step 5.3: validation crosstab
    if "term_states_fwd" in adata.obs.columns:
        ts_df = adata.obs[["celltype", "term_states_fwd"]].dropna()
        print("\nTerminal state — cell type crosstab:")
        ct = pd.crosstab(ts_df["term_states_fwd"], ts_df["celltype"])
        print(ct)

    # Step 5.4: terminal states UMAP
    try:
        g.plot_macrostates(which="terminal", basis="X_umap",
                           legend_loc="right margin",
                           title="Terminal States",
                           save=f"{out}/plots/spoke_cellrank_terminal.png",
                           show=False)
        print("Plot: terminal states UMAP")
    except Exception as e:
        print(f"[skip] terminal states plot: {e}")

    # Step 5.5: initial state
    try:
        g.predict_initial_states(n_states=1, allow_overlap=True)
        if g.initial_states is not None:
            print("\nInitial state:")
            print(g.initial_states.value_counts().to_string())
    except Exception as e:
        print(f"[skip] initial state: {e}")

except Exception as e:
    print(f"[skip] terminal/initial state computation: {e}")

# ---- Step 6: Fate probabilities -------------------------------------------
try:
    g.compute_fate_probabilities()
    print(f"\nFate probabilities computed.")
    fate_names = list(g.fate_probabilities.names)
    print(f"Terminal fates: {fate_names}")
    print(f"Fate probability matrix shape: {g.fate_probabilities.X.shape}")
except Exception as e:
    print(f"FATAL: compute_fate_probabilities failed: {e}")
    adata.write(f"{out}/objects/spoke_cellrank.h5ad")
    sys.exit(0)

# ---- Step 7: Visualize fate probabilities ----------------------------------

# Step 7.1: one UMAP per fate
try:
    fate_probs_arr = np.array(g.fate_probabilities)
    for i, fate in enumerate(fate_names):
        adata.obs[f"fate_{fate}"] = fate_probs_arr[:, i]

    n_fates = len(fate_names)
    fig, axes = plt.subplots(1, n_fates, figsize=(5 * n_fates, 5))
    if n_fates == 1:
        axes = [axes]
    for ax, fate in zip(axes, fate_names):
        sc.pl.embedding(adata, basis="X_umap", color=f"fate_{fate}",
                        color_map="RdYlBu_r", vmin=0, vmax=1,
                        title=f"Fate: {fate}", ax=ax, show=False)
    plt.tight_layout()
    plt.savefig(f"{out}/plots/spoke_cellrank_fate_probabilities.png", dpi=100)
    plt.close()
    print("Plot: fate probability UMAPs")
except Exception as e:
    print(f"[skip] fate probability UMAPs: {e}")

# Step 7.2: aggregate fate probabilities by cell type (validation heatmap)
try:
    fate_cols = [f"fate_{f}" for f in fate_names if f"fate_{f}" in adata.obs.columns]
    if fate_cols:
        fate_by_ct = adata.obs.groupby("celltype")[fate_cols].mean()
        fate_by_ct.columns = fate_names

        # PBMC diagnostic: flag B-cell domination (>80% of any cell type)
        print("\nMean fate probabilities by cell type:")
        print(fate_by_ct.round(3).to_string())
        dom = (fate_by_ct > 0.80).any(axis=0)
        if dom.any():
            dominated = dom[dom].index.tolist()
            print(f"\nWARNING: Fate(s) {dominated} dominate >80% of at least one cell type.")
            print("For a static PBMC dataset this reflects kNN topology, not developmental fate.")

        # Heatmap
        import seaborn as sns
        fig, ax = plt.subplots(figsize=(max(6, n_fates * 1.5),
                                        max(4, fate_by_ct.shape[0] * 0.5)))
        sns.heatmap(fate_by_ct, annot=True, fmt=".2f", cmap="YlOrRd",
                    vmin=0, vmax=1, ax=ax, linewidths=0.5)
        ax.set_title("Mean Fate Probability by Cell Type")
        ax.set_xlabel("Terminal Fate"); ax.set_ylabel("Cell Type")
        plt.tight_layout()
        plt.savefig(f"{out}/plots/spoke_cellrank_fate_heatmap.png", dpi=100)
        plt.close()
        print("Plot: fate probability heatmap")

        # Save table
        fate_by_ct.to_csv(f"{out}/tables/spoke_cellrank_fate_by_celltype.csv")
except Exception as e:
    print(f"[skip] fate probability heatmap: {e}")

# ---- Step 8: Lineage driver genes -----------------------------------------
try:
    drivers_list = []
    for fate in fate_names:
        fate_col = f"fate_{fate}"
        if fate_col not in adata.obs.columns:
            continue
        try:
            dr = g.compute_lineage_drivers(
                lineages=[fate],
                cluster_key="celltype",
                return_drivers=True,
                use_raw=False
            )
            if dr is not None:
                dr_df = dr.copy()
                dr_df["fate"] = fate
                drivers_list.append(dr_df.head(50))
                print(f"  Lineage drivers ({fate}): top 5 = "
                      f"{dr_df.index[:5].tolist()}")
        except Exception as e2:
            print(f"  [skip] lineage drivers ({fate}): {e2}")
            continue

    if drivers_list:
        all_drivers = pd.concat(drivers_list)
        all_drivers.to_csv(f"{out}/tables/spoke_cellrank_drivers.csv")
        print(f"Table: lineage drivers saved ({len(all_drivers)} rows)")

        # Step 8.3: top driver genes UMAP (top 2 per fate)
        driver_genes = []
        for fate in fate_names:
            sub = all_drivers[all_drivers["fate"] == fate]
            driver_genes.extend(sub.index[:2].tolist())
        driver_genes = [g for g in dict.fromkeys(driver_genes)
                        if g in adata.var_names][:8]
        if driver_genes:
            try:
                fig, axes = plt.subplots(2, 4, figsize=(20, 8))
                axes = axes.flatten()
                for i, gene in enumerate(driver_genes[:8]):
                    sc.pl.embedding(adata, basis="X_umap", color=gene,
                                    color_map="RdYlBu_r", ax=axes[i], show=False)
                for j in range(len(driver_genes), 8):
                    axes[j].axis("off")
                plt.suptitle("Top Lineage Driver Genes", y=1.02)
                plt.tight_layout()
                plt.savefig(f"{out}/plots/spoke_cellrank_driver_genes.png", dpi=100)
                plt.close()
                print(f"Plot: driver gene UMAPs ({', '.join(driver_genes)})")
            except Exception as e:
                print(f"[skip] driver gene UMAPs: {e}")
    else:
        # Create empty CSV so Snakemake output exists
        pd.DataFrame(columns=["gene", "fate"]).to_csv(
            f"{out}/tables/spoke_cellrank_drivers.csv")

except Exception as e:
    print(f"[skip] lineage drivers: {e}")
    pd.DataFrame(columns=["gene", "fate"]).to_csv(
        f"{out}/tables/spoke_cellrank_drivers.csv")

# ---- Save -----------------------------------------------------------------
adata.write(f"{out}/objects/spoke_cellrank.h5ad")
print("\ncellrank done ->", f"{out}/objects/spoke_cellrank.h5ad")
print("NOTE: For a static PBMC dataset, fate probabilities reflect kNN topology.")
print("      Interpret these results as a methods demonstration only.")
