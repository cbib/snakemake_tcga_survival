configfile: "config.yaml"

def load_tcga_cohorts(wildcards):
    with open(config["TCGA_cohorts"]) as f:
        TCGA_cohortS = [line.strip() for line in f if line.strip()]
    return  expand("<results>/rds/{cohort}_STAR_Counts.rds",
               cohort=TCGA_cohortS)

rule all_tcga:
    input:
        load_tcga_cohorts

rule TCGA_download:
    input:
        config["TCGA_cohorts"]
    output:
        "<results>/rds/{cohort}_STAR_Counts.rds"
    threads:
        config["resources"]["TCGA_download"]["threads"]
    resources:
        mem = config["resources"]["TCGA_download"]["mem"],
        time = config["resources"]["TCGA_download"]["time"]
    conda:
        "../../envs/tcga_smk.yml"
    log:
        "<logs>/TCGA_download_{cohort}.log"
    retries:
        3
    shell:
        """
        Rscript scripts/TCGA_download.R {wildcards.cohort} {output} > {log} 2>&1
        """

rule rule_separate_cohorts:
    input:
        "<results>/rds/TCGA-SKCM_STAR_Counts.rds"
    output:
        outfile_prim = "<results>/rds/TCGA-SKCM_met_STAR_Counts.rds",
        outfile_met = "<results>/rds/TCGA-SKCM_prim_STAR_Counts.rds"
    threads:
        config["resources"]["TCGA_download"]["threads"]
    conda:
        "../../envs/tcga_smk.yml"
    resources:
        mem = config["resources"]["TCGA_download"]["mem"],
        time = config["resources"]["TCGA_download"]["time"]
    log:
        "<logs>/TCGA_split_SKCM.log"
    shell:
        """
        Rscript scripts/TCGA_split_SKCM.R {input} {output.outfile_prim} {output.outfile_met} > {log} 2>&1
        """
