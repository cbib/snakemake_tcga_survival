library(DESeq2)
library(TCGAbiolinks)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
rds_file <- args[2]
biomart_file <- args[3]
outdir <- dirname(outfile)
GDC_dir <- dirname(outdir)

# create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
print("directory created")
# set working directory to output directory
setwd(GDC_dir)

# load data
data = readRDS(rds_file)
# remove "Solid Tissue Normal" samples
samples_to_keep <- colData(data)$sample_type != "Solid Tissue Normal"
data <- data[, samples_to_keep]
# read biomaRt file
genes = read.csv(biomart_file)
# create dataframe with data rownames, and then remove everything after "."
data_rownames <- data.frame(ensembl_gene_id = rownames(data))
data_rownames$ensembl_gene_id <- sub("\\..*$", "", data_rownames$ensembl_gene_id)
# merge with genes to get hgnc_symbol
data_merged <- merge(data_rownames, genes, by = "ensembl_gene_id")
# remove duplicates based on hgnc_symbol, keeping the first occurrence
data_merged <- data_merged[!duplicated(data_merged$hgnc_symbol), ]
# filter matrix to keep only rows in data_merged
data_filtered <- data[match(data_merged$ensembl_gene_id, sub("\\..*$", "", rownames(data))), ]
# convert rownames to hgnc_symbol
rownames(data_filtered) <- data_merged$hgnc_symbol
# normalize using DESeq2
dds <- DESeqDataSetFromMatrix(countData = assay(data_filtered),
                              colData = colData(data),
                              design = ~1)
# save normalized object
saveRDS(dds, file = outfile)