configfile: "config.yaml"

rule biomaRt_download:
    output:
        "<results>/biomart/biomart_protein_coding_genes.csv"
    threads:
        config["resources"]["biomaRt_download"]["threads"]
    resources:
        mem = config["resources"]["biomaRt_download"]["mem"],
        time = config["resources"]["biomaRt_download"]["time"]
    conda:
        "../../envs/tcga_smk.yml"
    log:
        "<logs>/biomaRt_download.log"
    retries:
        3
    shell:
        """
        Rscript scripts/biomaRt_download.R {output} > {log} 2>&1
        """
