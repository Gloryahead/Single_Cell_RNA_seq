"""
Snakemake workflow for GSE174609 scRNA-seq pipeline (Parts 1-17).

Run (conda envs, no containers):
  snakemake --profile profile/slurm -j 8 all

Run (Apptainer containers — build them first with containers/build_containers.sh):
  snakemake --profile profile/slurm --use-singularity -j 8 all

Run a single target:
  snakemake --profile profile/slurm results/objects/04_clustered.rds
"""

import yaml, os

# ── Config & paths ────────────────────────────────────────────────────────────
cfg  = yaml.safe_load(open("config.yaml"))
OUT  = cfg["project"]["outdir"]
# Code (scripts, config, envs) lives in home; all data/results live on xdisk.
CODE = os.environ.get("SCRNA_PROJECT", "/home/u11/maarowosegbe/Single_Cell_RNA_seq")
DATA = os.environ.get("SCRNA_DATA",    "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")

SAMPLES  = [s["id"] for s in cfg["samples"]]
HEALTHY  = [s["id"] for s in cfg["samples"] if s["condition"] == "Healthy"]
POST     = [s["id"] for s in cfg["samples"] if s["condition"] == "Post"]
VEL_SAMP = HEALTHY + POST      # samples to run STARsolo on (for velocity)

# Container SIF files are large — store on xdisk, not in home
R_SIF  = f"{DATA}/containers/scrna_r.sif"
PY_SIF = f"{DATA}/containers/scrna_python.sif"

# ── Helpers ───────────────────────────────────────────────────────────────────
def rscript(script):
    return f"SCRNA_PROJECT={CODE} SCRNA_DATA={DATA} Rscript {script}"

def pyscript(script):
    return f"SCRNA_PROJECT={CODE} SCRNA_DATA={DATA} python {script}"

# ── Master target ─────────────────────────────────────────────────────────────
localrules: qc_report, setup_networks

rule all:
    input:
        # QC report (generated from Cell Ranger metrics_summary.csv for every sample)
        f"{OUT}/qc_report.html",
        # Core pipeline
        f"{OUT}/objects/03_filtered.rds",
        f"{OUT}/objects/04_clustered.rds",
        f"{OUT}/objects/05_annotated.rds",
        f"{OUT}/tables/06_de_post_vs_healthy.csv",
        # Spokes
        # spoke_monocle3.rds    — monocle3 not installed in scrna_r.sif; rebuild container
        f"{OUT}/objects/spoke_slingshot.rds",
        # spoke_cellchat.rds    — CellChat not installed in scrna_r.sif; rebuild container
        # spoke_nichenet_ligands.csv — nichenetr not installed in scrna_r.sif; rebuild container
        # spoke_copykat.rds     — copykat not installed in scrna_r.sif; rebuild container
        f"{OUT}/objects/spoke_hdwgcna.rds",
        f"{OUT}/objects/spoke_velocity.h5ad",
        f"{OUT}/tables/spoke_velocity_top_genes.csv",
        # spoke_cellrank.h5ad — MPICH/OFI MPI init fails in Singularity on Puma (no IB fabric);
        #   C-level Abort() is uncatchable; also CellRank is inappropriate for static PBMCs
        # spoke_scplotter_02_umap_celltypes.png — scplotter not installed in scrna_r.sif; rebuild container
        f"{OUT}/tables/spoke_tf_activity.csv",
        f"{OUT}/tables/spoke_pathway_activity.csv",

# ── QC report (runs locally after Cell Ranger; no SLURM job needed) ──────────
rule qc_report:
    input:
        expand(f"{OUT}/counts/{{sample}}/outs/metrics_summary.csv", sample=SAMPLES)
    output:
        html = f"{OUT}/qc_report.html",
        csv  = f"{OUT}/qc_report.csv",
    shell:
        pyscript("code/qc_report.py")

# ── Part 1: Cell Ranger alignment ─────────────────────────────────────────────
rule cellranger_count:
    input:
        r1 = f"{DATA}/data/fastq/{{sample}}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2 = f"{DATA}/data/fastq/{{sample}}/{{sample}}_S1_L001_R2_001.fastq.gz",
    output:
        directory(f"{OUT}/counts/{{sample}}/outs/filtered_feature_bc_matrix"),
    params:
        ref = f"{DATA}/references/cellranger/refdata-gex-GRCh38-2024-A",
        out = f"{OUT}/counts",
    resources:
        cpus_per_task = 22,
        mem_mb        = 112640,  # 22 × 5120 MB = 110 GB (matches --localmem=105)
        runtime       = 1440,    # 24 h
    # Cell Ranger is NOT containerised (10x license restriction)
    shell:
        """
        source {DATA}/software/activate_cellranger.sh
        cd {params.out}
        rm -rf {wildcards.sample}
        cellranger count \
          --id={wildcards.sample} \
          --transcriptome={params.ref} \
          --fastqs={DATA}/data/fastq/{wildcards.sample} \
          --sample={wildcards.sample} \
          --localcores={resources.cpus_per_task} \
          --localmem=105 \
          --create-bam=false
        """

# ── Part 2: QC filtering ──────────────────────────────────────────────────────
rule qc_filter:
    input:
        filtered = expand(f"{OUT}/counts/{{sample}}/outs/filtered_feature_bc_matrix",
                          sample=SAMPLES),
        raw      = expand(f"{OUT}/counts/{{sample}}/outs/raw_feature_bc_matrix",
                          sample=SAMPLES),
    output:
        f"{OUT}/objects/03_filtered.rds",
    resources:
        cpus_per_task = 8,
        mem_mb        = 51200,   # 8 × 5120 MB → 50 GB (SoupX loads raw + filtered per sample)
        runtime       = 240,     # ScaleData + RunUMAP on 12-sample merged object needs ~3 h
    singularity: R_SIF
    shell:
        rscript("code/03_qc_filter.R")

# ── Part 3: Integration & clustering ─────────────────────────────────────────
rule integrate_cluster:
    input:  f"{OUT}/objects/03_filtered.rds"
    output: f"{OUT}/objects/04_clustered.rds"
    resources:
        cpus_per_task = 8,
        mem_mb        = 40960,   # 8 × 5120 MB = 40 GB
        runtime       = 120,
    singularity: R_SIF
    shell:
        rscript("code/04_cluster_integrate.R")

# ── Part 4: Cell type annotation ─────────────────────────────────────────────
rule annotate:
    input:  f"{OUT}/objects/04_clustered.rds"
    output:
        rds      = f"{OUT}/objects/05_annotated.rds",
        meta_csv = f"{OUT}/objects/05_annotated_meta.csv",
        umap_csv = f"{OUT}/objects/05_annotated_umap.csv",
    resources:
        cpus_per_task = 20,
        mem_mb        = 102400,  # 20 × 5120 MB = 100 GB (3× SingleR + scType + FindAllMarkers)
        runtime       = 360,     # 6 h: Monaco + DICE add ~2× SingleR runtime
    singularity: R_SIF
    shell:
        rscript("code/05_annotate.R")

# ── Part 5: DE + pathway ──────────────────────────────────────────────────────
rule dge_pathway:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/tables/06_de_post_vs_healthy.csv"
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB
        runtime       = 360,     # 6 h: FindMarkers per cell type × 3 contrasts adds ~90 min
    singularity: R_SIF
    shell:
        rscript("code/06_dge_pathway.R")

# ── Part 7: Trajectory ────────────────────────────────────────────────────────
rule monocle3:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_monocle3.rds"
    resources:
        cpus_per_task = 8,
        mem_mb        = 40960,   # 8 × 5120 MB = 40 GB
        runtime       = 120,
    singularity: R_SIF
    shell:
        rscript("code/spokes/trajectory_monocle3.R")

rule slingshot:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_slingshot.rds"
    resources:
        cpus_per_task = 8,
        mem_mb        = 40960,   # 8 × 5120 MB = 40 GB
        runtime       = 120,     # gene-pseudotime Spearman correlations add ~30-60 min
    singularity: R_SIF
    shell:
        rscript("code/spokes/trajectory_slingshot.R")

# ── Part 9: CellChat ──────────────────────────────────────────────────────────
rule cellchat:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_cellchat.rds"
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB (CellChat is heavy with 106K cells)
        runtime       = 240,
    singularity: R_SIF
    shell:
        rscript("code/spokes/ccc_cellchat.R")

# ── Part 10: NicheNet ─────────────────────────────────────────────────────────
rule nichenet:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/tables/spoke_nichenet_ligands.csv"
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB
        runtime       = 240,
    singularity: R_SIF
    shell:
        rscript("code/spokes/ccc_nichenet.R")

# ── Part 11: CopyKAT ──────────────────────────────────────────────────────────
rule copykat:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_copykat.rds"
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB
        runtime       = 360,
    singularity: R_SIF
    shell:
        rscript("code/spokes/cnv_copykat.R")

# ── Part 12: hdWGCNA ──────────────────────────────────────────────────────────
rule hdwgcna:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_hdwgcna.rds"
    resources:
        cpus_per_task = 18,
        mem_mb        = 92160,   # 18 × 5120 MB = 90 GB (Seurat obj ~40 GB + WGCNA TOM + ModuleEigengenes OOMed at 60 GB)
        runtime       = 240,
    singularity: R_SIF
    shell:
        rscript("code/spokes/coexpression_hdwgcna.R")

# ── Part 13: STARsolo (per sample) ────────────────────────────────────────────
rule starsolo:
    input:
        r1 = f"{DATA}/data/fastq/{{sample}}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2 = f"{DATA}/data/fastq/{{sample}}/{{sample}}_S1_L001_R2_001.fastq.gz",
    output:
        directory(f"{OUT}/velocity/{{sample}}/Solo.out"),
    params:
        ref      = f"{DATA}/references/star_index_GRCh38",
        gtf      = f"{DATA}/references/GRCh38/Homo_sapiens.GRCh38.110.gtf",
        wl       = f"{DATA}/software/3M-february-2018.txt",
        out_pref = f"{OUT}/velocity/{{sample}}/",
    resources:
        cpus_per_task = 16,
        mem_mb        = 81920,   # 16 × 5120 MB = 80 GB
        runtime       = 480,
    singularity: R_SIF   # STAR is in the R env
    shell:
        """
        mkdir -p {params.out_pref}
        STAR \
          --soloType CB_UMI_Simple \
          --soloCBwhitelist {params.wl} \
          --soloCBstart 1 --soloCBlen 16 \
          --soloUMIstart 17 --soloUMIlen 12 \
          --soloFeatures Gene Velocyto \
          --soloOutFormatFeaturesGeneField3 "Gene Expression" \
          --genomeDir {params.ref} \
          --sjdbGTFfile {params.gtf} \
          --readFilesIn {input.r2} {input.r1} \
          --readFilesCommand zcat \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
          --outFileNamePrefix {params.out_pref} \
          --runThreadN {resources.cpus_per_task} \
          --limitBAMsortRAM 60000000000
        """

# ── Part 13: scVelo ───────────────────────────────────────────────────────────
rule scvelo:
    input:
        meta_csv = f"{OUT}/objects/05_annotated_meta.csv",
        umap_csv = f"{OUT}/objects/05_annotated_umap.csv",
        starsolo = expand(f"{OUT}/velocity/{{sample}}/Solo.out",
                          sample=VEL_SAMP),
    output:
        h5ad = f"{OUT}/objects/spoke_velocity.h5ad",
        csv  = f"{OUT}/tables/spoke_velocity_top_genes.csv",
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB (stochastic mode + moments on 32K cells)
        runtime       = 480,     # extra 2 h: per-condition plots + HVG selection + top-genes rank
    singularity: PY_SIF
    shell:
        pyscript("code/spokes/velocity_scvelo.py")

# ── Part 14: CellRank ─────────────────────────────────────────────────────────
rule cellrank:
    input:  f"{OUT}/objects/spoke_velocity.h5ad"
    output:
        h5ad    = f"{OUT}/objects/spoke_cellrank.h5ad",
        drivers = f"{OUT}/tables/spoke_cellrank_drivers.csv",
    resources:
        cpus_per_task = 8,
        mem_mb        = 40960,   # 8 × 5120 MB = 40 GB
        runtime       = 240,
    singularity: PY_SIF
    shell:
        pyscript("code/spokes/fate_cellrank.py")

# ── Part 15: scplotter ────────────────────────────────────────────────────────
rule scplotter:
    input:  f"{OUT}/objects/05_annotated.rds"
    output:
        umap     = f"{OUT}/plots/spoke_scplotter_02_umap_celltypes.png",
        comp     = f"{OUT}/plots/spoke_scplotter_13_composition_by_condition.png",
    resources:
        cpus_per_task = 8,
        mem_mb        = 40960,   # 8 × 5120 MB = 40 GB
        runtime       = 120,     # FindAllMarkers + condition DE add ~45-60 min
    singularity: R_SIF
    shell:
        rscript("code/spokes/viz_scplotter.R")

# ── Network pre-download (localrule — runs on orchestrator node, has internet) ─
# Downloads CollecTRI and PROGENy CSVs before Parts 16/17 compute jobs start.
# Being a localrule means it runs on the Snakemake orchestrator node (the job
# submitted via snakemake.slurm), never dispatched to a child SLURM job.
# If you prefer to run it manually first: sbatch code/setup_networks.slurm
rule setup_networks:
    output:
        # Only PROGENy is mandatory — no fallback exists if it fails.
        # CollecTRI is optional: if collectri_human.csv is absent, grn_decoupler.R
        # falls back to GTRD (MSigDB C3). Declaring it here would cause
        # MissingOutputException whenever the OmnipathR static-table bug fires.
        progeny = f"{OUT}/references/progeny_human.csv",
    singularity: R_SIF
    shell:
        rscript("code/setup_networks.R")

# ── Part 16: decoupleR / CollecTRI ────────────────────────────────────────────
rule decoupler:
    input:
        rds     = f"{OUT}/objects/05_annotated.rds",
        progeny = f"{OUT}/references/progeny_human.csv",  # ensures setup_networks runs first
    output:
        pb  = f"{OUT}/tables/spoke_tf_activity.csv",
        ct  = f"{OUT}/tables/spoke_tf_activity_by_celltype.csv",
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB (pseudobulk + 500-cell downsample ULM)
        runtime       = 120,     # 2 ULM runs + 4 plots (no download — CSV pre-loaded)
    singularity: R_SIF
    shell:
        rscript("code/spokes/grn_decoupler.R")

# ── Part 17: PROGENy ─────────────────────────────────────────────────────────
rule progeny:
    input:
        rds     = f"{OUT}/objects/05_annotated.rds",
        progeny = f"{OUT}/references/progeny_human.csv",
    output:
        pb  = f"{OUT}/tables/spoke_pathway_activity.csv",
        ct  = f"{OUT}/tables/spoke_pathway_activity_by_celltype.csv",
    resources:
        cpus_per_task = 12,
        mem_mb        = 61440,   # 12 × 5120 MB = 60 GB (pseudobulk + 500-cell downsample MLM)
        runtime       = 90,      # 2 MLM runs + 3 plots (no download — CSV pre-loaded)
    singularity: R_SIF
    shell:
        rscript("code/spokes/pathway_progeny.R")
