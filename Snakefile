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
PROJ = os.environ.get("SCRNA_PROJECT",
                      "/xdisk/haining/maarowosegbe/Single_Cell_RNA_seq")

SAMPLES  = [s["id"] for s in cfg["samples"]]
HEALTHY  = [s["id"] for s in cfg["samples"] if s["condition"] == "Healthy"]
POST     = [s["id"] for s in cfg["samples"] if s["condition"] == "Post"]
VEL_SAMP = HEALTHY + POST      # samples to run STARsolo on (for velocity)

# Container paths (only used with --use-singularity)
R_SIF  = f"{PROJ}/containers/scrna_r.sif"
PY_SIF = f"{PROJ}/containers/scrna_python.sif"

# ── Helpers ───────────────────────────────────────────────────────────────────
def rscript(script):
    return f"SCRNA_PROJECT={PROJ} Rscript {script}"

def pyscript(script):
    return f"SCRNA_PROJECT={PROJ} python {script}"

# ── Master target ─────────────────────────────────────────────────────────────
rule all:
    input:
        # Core pipeline
        f"{OUT}/objects/03_filtered.rds",
        f"{OUT}/objects/04_clustered.rds",
        f"{OUT}/objects/05_annotated.rds",
        f"{OUT}/tables/06_de_post_vs_healthy.csv",
        # Spokes
        f"{OUT}/objects/spoke_monocle3.rds",
        f"{OUT}/objects/spoke_slingshot.rds",
        f"{OUT}/objects/spoke_cellchat.rds",
        f"{OUT}/tables/spoke_nichenet_ligands.csv",
        f"{OUT}/objects/spoke_copykat.rds",
        f"{OUT}/objects/spoke_hdwgcna.rds",
        f"{OUT}/objects/spoke_velocity.h5ad",
        f"{OUT}/objects/spoke_cellrank.h5ad",
        f"{OUT}/plots/spoke_scplotter_umap.png",
        f"{OUT}/tables/spoke_tf_activity.csv",
        f"{OUT}/tables/spoke_pathway_activity.csv",

# ── Part 1: Cell Ranger alignment ─────────────────────────────────────────────
rule cellranger_count:
    input:
        r1 = f"{PROJ}/data/fastq/{{sample}}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2 = f"{PROJ}/data/fastq/{{sample}}/{{sample}}_S1_L001_R2_001.fastq.gz",
    output:
        directory(f"{OUT}/counts/{{sample}}/outs/filtered_feature_bc_matrix"),
    params:
        ref = f"{PROJ}/references/cellranger/refdata-gex-GRCh38-2024-A",
        out = f"{OUT}/counts",
    resources:
        mem_mb   = 110000,
        runtime  = 1440,    # 24 h
        cpus_per_task = 16,
    # Cell Ranger is NOT containerised (10x license restriction)
    shell:
        """
        source {PROJ}/software/activate_cellranger.sh
        cd {params.out}
        rm -rf {wildcards.sample}
        cellranger count \
          --id={wildcards.sample} \
          --transcriptome={params.ref} \
          --fastqs={PROJ}/data/fastq/{wildcards.sample} \
          --sample={wildcards.sample} \
          --localcores={resources.cpus_per_task} \
          --localmem=105 \
          --create-bam=false
        """

# ── Part 2: QC filtering ──────────────────────────────────────────────────────
rule qc_filter:
    input:
        expand(f"{OUT}/counts/{{sample}}/outs/filtered_feature_bc_matrix",
               sample=SAMPLES),
    output:
        f"{OUT}/objects/03_filtered.rds",
    resources:
        mem_mb   = 65536,
        runtime  = 120,
        cpus_per_task = 8,
    singularity: R_SIF
    shell:
        rscript("code/03_qc_filter.R")

# ── Part 3: Integration & clustering ─────────────────────────────────────────
rule integrate_cluster:
    input:  f"{OUT}/objects/03_filtered.rds"
    output: f"{OUT}/objects/04_clustered.rds"
    resources:
        mem_mb   = 65536,
        runtime  = 120,
        cpus_per_task = 8,
    singularity: R_SIF
    shell:
        rscript("code/04_cluster_integrate.R")

# ── Part 4: Cell type annotation ─────────────────────────────────────────────
rule annotate:
    input:  f"{OUT}/objects/04_clustered.rds"
    output:
        rds  = f"{OUT}/objects/05_annotated.rds",
        h5ad = f"{OUT}/objects/05_annotated.h5ad",
    resources:
        mem_mb   = 65536,
        runtime  = 120,
        cpus_per_task = 8,
    singularity: R_SIF
    shell:
        rscript("code/05_annotate.R")

# ── Part 5: DE + pathway ──────────────────────────────────────────────────────
rule dge_pathway:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/tables/06_de_post_vs_healthy.csv"
    resources:
        mem_mb   = 65536,
        runtime  = 180,
        cpus_per_task = 8,
    singularity: R_SIF
    shell:
        rscript("code/06_dge_pathway.R")

# ── Part 7: Trajectory ────────────────────────────────────────────────────────
rule monocle3:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_monocle3.rds"
    resources:
        mem_mb   = 65536,
        runtime  = 120,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/trajectory_monocle3.R")

rule slingshot:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_slingshot.rds"
    resources:
        mem_mb   = 32768,
        runtime  = 60,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/trajectory_slingshot.R")

# ── Part 9: CellChat ──────────────────────────────────────────────────────────
rule cellchat:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_cellchat.rds"
    resources:
        mem_mb   = 65536,
        runtime  = 240,
        cpus_per_task = 8,
    singularity: R_SIF
    shell:
        rscript("code/spokes/ccc_cellchat.R")

# ── Part 10: NicheNet ─────────────────────────────────────────────────────────
rule nichenet:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/tables/spoke_nichenet_ligands.csv"
    resources:
        mem_mb   = 65536,
        runtime  = 240,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/ccc_nichenet.R")

# ── Part 11: CopyKAT ──────────────────────────────────────────────────────────
rule copykat:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_copykat.rds"
    resources:
        mem_mb   = 65536,
        runtime  = 360,
        cpus_per_task = 8,
    singularity: R_SIF
    shell:
        rscript("code/spokes/cnv_copykat.R")

# ── Part 12: hdWGCNA ──────────────────────────────────────────────────────────
rule hdwgcna:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/objects/spoke_hdwgcna.rds"
    resources:
        mem_mb   = 65536,
        runtime  = 240,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/coexpression_hdwgcna.R")

# ── Part 13: STARsolo (per sample) ────────────────────────────────────────────
rule starsolo:
    input:
        r1 = f"{PROJ}/data/fastq/{{sample}}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2 = f"{PROJ}/data/fastq/{{sample}}/{{sample}}_S1_L001_R2_001.fastq.gz",
    output:
        directory(f"{OUT}/velocity/{{sample}}/Solo.out"),
    params:
        ref      = f"{PROJ}/references/star_index_GRCh38",
        gtf      = f"{PROJ}/references/GRCh38/Homo_sapiens.GRCh38.110.gtf",
        wl       = f"{PROJ}/software/3M-february-2018.txt",
        out_pref = f"{OUT}/velocity/{{sample}}/",
    resources:
        mem_mb   = 80000,
        runtime  = 480,
        cpus_per_task = 16,
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
          --soloOutDir {params.out_pref}/Solo.out \
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
        h5ad    = f"{OUT}/objects/05_annotated.h5ad",
        starsolo = expand(f"{OUT}/velocity/{{sample}}/Solo.out",
                          sample=VEL_SAMP),
    output: f"{OUT}/objects/spoke_velocity.h5ad"
    resources:
        mem_mb   = 65536,
        runtime  = 360,
        cpus_per_task = 8,
    singularity: PY_SIF
    shell:
        pyscript("code/spokes/velocity_scvelo.py")

# ── Part 14: CellRank ─────────────────────────────────────────────────────────
rule cellrank:
    input:  f"{OUT}/objects/spoke_velocity.h5ad"
    output: f"{OUT}/objects/spoke_cellrank.h5ad"
    resources:
        mem_mb   = 65536,
        runtime  = 240,
        cpus_per_task = 8,
    singularity: PY_SIF
    shell:
        pyscript("code/spokes/fate_cellrank.py")

# ── Part 15: scplotter ────────────────────────────────────────────────────────
rule scplotter:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/plots/spoke_scplotter_umap.png"
    resources:
        mem_mb   = 32768,
        runtime  = 60,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/viz_scplotter.R")

# ── Part 16: decoupleR / CollecTRI ────────────────────────────────────────────
rule decoupler:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/tables/spoke_tf_activity.csv"
    resources:
        mem_mb   = 32768,
        runtime  = 60,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/grn_decoupler.R")

# ── Part 17: PROGENy ─────────────────────────────────────────────────────────
rule progeny:
    input:  f"{OUT}/objects/05_annotated.rds"
    output: f"{OUT}/tables/spoke_pathway_activity.csv"
    resources:
        mem_mb   = 32768,
        runtime  = 60,
        cpus_per_task = 4,
    singularity: R_SIF
    shell:
        rscript("code/spokes/pathway_progeny.R")
