library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
tcga_file <- args[1]
dir_list <- as.vector(read.table(tcga_file, header = FALSE)[,1])
dir_list <- c(dir_list, "TCGA-SKCM_prim", "TCGA-SKCM_met")
output <- args[2]
out_path_survival <- dirname(output)
# set working directory to survival
setwd(out_path_survival)

# helper function to merge survival tables
merge_survival_pvalues <- function(input_filename = "survival_pval.tsv",
                                   output_filename = "survival_pval_merged.xlsx") {
  
  # create blank excel workbook
  OUT <- createWorkbook()
  
  # loop through dirs to summarize the information
  for (dir in dir_list){
    # add one sheet per cohort
    addWorksheet(OUT, dir)
    # read file to save
    file <- file.path(dir, input_filename)
    data <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    # write the survival data
    writeData(OUT, sheet = dir, x = data, rowNames = TRUE, colNames = TRUE)
    if (dir == dir_list[length(dir_list)]){
      saveWorkbook(OUT, output_filename)
      print("Merged file saved")
    }
  }
}

# merge survival pval tables
merge_survival_pvalues(input_filename = "survival_pval.tsv",
                       output_filename = "survival_pval_merged.xlsx")
# merge filtered tables
merge_survival_pvalues(input_filename = "survival_pval_filtered.tsv",
                       output_filename = "survival_pval_filtered_merged.xlsx")

# save also a reorganized excel with the summary per gene/signature
# empty list to store excel each sheet data
data_list <- list()
# loop through tsvs
for (i in seq_along(dir_list)) {
  # get file path 
  file_path <- file.path(dir_list[i], "survival_pval.tsv")
  # read tsv
  dat <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  # save tsv in list
  data_list[[dir_list[i]]] <- dat
}

# get gene/signature names
all_signatures <- rownames(data_list[[1]])

# create empty excel book
OUT <- createWorkbook()

# build a table where rows = dirs and columns = percentiles for each signature
for (signature in all_signatures) {
  # container for one feature across all directories
  df_row <- data.frame(matrix(nrow = length(data_list), ncol = ncol(data_list[[1]])))
  colnames(df_row) <- colnames(data_list[[1]])
  rownames(df_row) <- names(data_list)
  # fill from each dataset
  for (dir in names(data_list)) {
    df_row[dir, ] <- data_list[[dir]][signature, ]
  }
  # add worksheet
  addWorksheet(OUT, signature)
  writeData(OUT, sheet = signature, x = df_row,
            rowNames = TRUE, colNames = TRUE)
}

# save workbook
saveWorkbook(OUT, "merged_per_signature.xlsx")
message("Saved merged_per_signature.xlsx")
