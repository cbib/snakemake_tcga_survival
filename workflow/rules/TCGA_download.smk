configfile: "config.yaml"

def load_tcga_projects(wildcards):
    with open(config["TCGA_cohorts"]) as f:
        TCGA_PROJECTS = [line.strip() for line in f if line.strip()]
    return  expand("<results>/rds/{project}_STAR_Counts.rds",
               project=TCGA_PROJECTS)

rule all_tcga:
    input:
        load_tcga_projects

rule TCGA_download:
    input:
        config["TCGA_cohorts"]
    output:
        "<results>/rds/{project}_STAR_Counts.rds"
    threads:
        config["resources"]["TCGA_download"]["threads"]
    resources:
        mem = config["resources"]["TCGA_download"]["mem"],
        time = config["resources"]["TCGA_download"]["time"]
    conda:
        "../../envs/tcga_smk.yml"
    log:
        "<logs>/TCGA_download_{project}.log"
    shell:
        """
        Rscript scripts/TCGA_download.R {wildcards.project} {output} > {log} 2>&1
        """
