suppressPackageStartupMessages({
    library(TCGAbiolinks)
})

args <- commandArgs(trailingOnly = TRUE)
project <- args[1]
outfile <- args[2]
outdir <- dirname(outfile)
GDC_dir <- dirname(outdir)

# create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# set working directory to output directory
setwd(GDC_dir)

# create query for project
query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

# download the project with retry logic
GDCdownload(query = query,
            method = "api",
            files.per.chunk = NULL)
data <- GDCprepare(query = query)

# save the data as .rds file
cat("Saving data for project", project, "to", outfile, "\n")
saveRDS(data, file = outfile)
