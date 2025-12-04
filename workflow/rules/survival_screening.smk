configfile: "config.yaml"
import os

def survival_output(wildcards):
    with open(config["TCGA_cohorts"]) as f:
        TCGA_PROJECTS = [line.strip() for line in f if line.strip()]
    # Automatically add SKCM_prim and SKCM_met
    if "TCGA-SKCM" in TCGA_PROJECTS:
        TCGA_PROJECTS += ["TCGA-SKCM_prim", "TCGA-SKCM_met"]
    
    return  expand("<results>/screening/survival/{project}/survival_pval_filtered.tsv",
               project=TCGA_PROJECTS)

rule all_survival:
    input:
        survival_output

rule survival_screening:
    input:
        signatures_file = os.path.abspath(config["signatures_file"]),
        cohorts = config["TCGA_cohorts"],
        deseq2_file = "<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds"
    params:
        THRESHOLD = config["THRESHOLD"],
        DPI = config["DPI"]
    output:
        "<results>/screening/survival/{project}/survival_pval_filtered.tsv"
    threads:
        config["resources"]["survival_screening"]["threads"]
    resources:
        mem = config["resources"]["survival_screening"]["mem"],
        time = config["resources"]["survival_screening"]["time"]
    conda:
        "../../envs/tcga_smk.yml"
    log:
        "<logs>/survival_screening_{project}.log"
    shell:
        """
        Rscript scripts/survival_screening.R {params.DPI} {params.THRESHOLD} {input.signatures_file} {wildcards.project} {output} > {log} 2>&1
        """
