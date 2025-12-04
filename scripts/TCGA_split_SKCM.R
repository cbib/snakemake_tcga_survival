suppressPackageStartupMessages({
    library(TCGAbiolinks)
})

args <- commandArgs(trailingOnly = TRUE)
inputfile <- args[1]
outfile_prim <- args[2]
outfile_met <- args[3]

# read SKCM cohort
data <- readRDS(inputfile)
# save object for the 103 primary melanomas
data_prim = data[, data$definition == "Primary solid Tumor"]
cat("Saving SKCM_prim data to", outfile_prim, "\n")
saveRDS(data_prim, file = outfile_prim)

# save object for the 369 metastatic melanomas
data_met = data[, data$definition %in% c("Metastatic", "Additional Metastatic")]
cat("Saving SKCM_prim data to", outfile_met, "\n")
saveRDS(data_met, file = outfile_met)
