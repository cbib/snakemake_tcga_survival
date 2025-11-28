configfile: "config.yaml"

def deseq2_output(wildcards):
    with open(config["TCGA_cohorts"]) as f:
        TCGA_PROJECTS = [line.strip() for line in f if line.strip()]
    return  expand("<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds",
               project=TCGA_PROJECTS)

rule all_deseq2:
    input:
        deseq2_output

rule DESeq2_normalization:
    input:
        cohorts = config["TCGA_cohorts"],
        rds_file = "<results>/rds/{project}_STAR_Counts.rds",
        biomart_file = "<results>/biomart/biomart_protein_coding_genes.csv"
    output:
        "<results>/DESeq2_normalized/{project}_STAR_Counts_DESeq2.rds"
    threads:
        config["resources"]["DESeq2_normalization"]["threads"]
    resources:
        mem = config["resources"]["DESeq2_normalization"]["mem"],
        time = config["resources"]["DESeq2_normalization"]["time"]
    conda:
        "../../envs/tcga_smk.yml"
    log:
        "<logs>/DESeq2_normalization_{project}.log"
    shell:
        """
        Rscript scripts/DESeq2_normalization.R {output} {input.rds_file} {input.biomart_file} > {log} 2>&1
        """
