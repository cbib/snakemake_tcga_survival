library(biomaRt)
library(curl)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
outdir <- dirname(outfile)

# create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Use the Ensembl Genes dataset
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Retrieve mapping for protein-coding genes
genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters    = "biotype",
  values     = "protein_coding",
  mart       = mart
)

# save file
write.csv(genes, file = outfile, row.names = FALSE)
