from snakemake.utils import validate
import pandas as pd
import sys
# configuration
configfile: "config.yaml"
# validate(config, schema="schemas/config.schema.yaml")

# include rule bundles
include: "workflow/rules/biomaRt_download.smk"
include: "workflow/rules/TCGA_download.smk"
include: "workflow/rules/DESeq2_normalization.smk"
include: "workflow/rules/survival_screening.smk"
include: "workflow/rules/merge_survival_results.smk"

# helper functions
# def survival_output(wildcards):
#     with open(config["TCGA_cohorts"]) as f:
#         TCGA_PROJECTS = [line.strip() for line in f if line.strip()]
#     return  expand("<results>/screening/survival/{project}/survival_pval_filtered.tsv",
#                project=TCGA_PROJECTS)

rule all:
    input:
        "<results>/screening/survival/merged_per_signature.xlsx"