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

rule all_merge:
    input:
        "<results>/screening/survival/merged_per_signature.xlsx"

rule merge_survival_results:
    input:
        cohorts = os.path.abspath(config["TCGA_cohorts"]),
        tsv_files = lambda wildcards: survival_output(wildcards)
    # params:
    #     out_path_survival = "<results>/screening/survival/"
    output:
        "<results>/screening/survival/merged_per_signature.xlsx"
    threads:
        config["resources"]["merge_survival_results"]["threads"]
    resources:
        mem = config["resources"]["merge_survival_results"]["mem"],
        time = config["resources"]["merge_survival_results"]["time"]
    conda:
        "../../envs/merge_smk.yml"
    log:
        "<logs>/merge_survival_results.log"
    shell:
        """
        Rscript scripts/merge_survival_results.R {input.cohorts} {output} > {log} 2>&1
        """
