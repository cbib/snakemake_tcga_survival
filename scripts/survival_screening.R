library(DESeq2)
library(TCGAbiolinks)
library(ggplot2)
library(dplyr)
library(GSVA)
library(survminer)
library(survival)

# define parameters
args <- commandArgs(trailingOnly = TRUE)
DPI <- args[1] %>% as.numeric()
THRESHOLD <- args[2] %>% as.numeric()
signatures_file <- args[3]
project <- args[4]
survival_table <- args[5]
GDC_dir <- survival_table %>% dirname() %>% dirname() %>% dirname() %>% dirname() 
# correct aspect ratio of plots
FACTOR = DPI / 100

# set working directory to output directory
setwd(GDC_dir)
# create PCA/project and survival/project directories
dir.create(paste0("./screening/PCA/", project), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0("./screening/survival/", project), showWarnings = FALSE,  recursive = FALSE)


#################################################################
###################### SURVIVAL SCREENING  ######################
#################################################################

# load normalized data
dds = readRDS(paste0("./DESeq2_normalized/", project,"_STAR_Counts_DESeq2.rds"))
# read gene_signatures.txt (tab separated)
gene_signatures <- read.delim(signatures_file, header = TRUE, stringsAsFactors = FALSE)
# load vector with percentiles
percentiles <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.5)

# classify signatures in < 2 genes and >= 2 genes (column values)
single_gene_signatures <- gene_signatures[sapply(gene_signatures, function(x) length(unlist(strsplit(x, ","))) < 2)]
multiple_gene_signatures <- gene_signatures[sapply(gene_signatures, function(x) length(unlist(strsplit(x, ","))) >= 2)]
# identify if any of the multiple gene signatures end with _UP or _DOWN and they have the same character string before that suffix
names_multiple <- names(multiple_gene_signatures)
base_names <- sub("(_UP|_DOWN)$", "", names_multiple)
duplicated_base_names <- base_names[duplicated(base_names)]
# extract normlalized counts for ssGSEA
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
# perform ssgsea for multiple gene signatures
scores <- gsva(norm_counts, as.list(multiple_gene_signatures), method = "ssgsea", ssgsea.norm = TRUE)
scores <- as.data.frame(t(scores))

# combine _UP and _DOWN signatures with the same base name
for (base_name in duplicated_base_names) {
  up_name <- paste0(base_name, "_UP")
  down_name <- paste0(base_name, "_DOWN")
  if (up_name %in% colnames(scores) & down_name %in% colnames(scores)) {
    combined_scores <- scores[, up_name] - scores[, down_name]
    scores <- cbind(scores, combined_scores)
    colnames(scores)[ncol(scores)] <- paste0(base_name, "_COMBINED")
  }
}

# stratify using percentiles and save categorical variables in colData
for (signature in colnames(scores)) {
  for (percentile in percentiles) {
    threshold_low <- quantile(scores[, signature], probs = percentile)
    threshold_high <- quantile(scores[, signature], probs = 1 - percentile)
    categorical_var <- ifelse(scores[, signature] <= threshold_low, "Low",
                              ifelse(scores[, signature] >= threshold_high, "High", "Intermediate"))
    col_name <- paste0(signature, "_", percentile * 100, "pct")
    colData(dds)[, col_name] <- categorical_var
  }
}

# add the single_gene_signatures to the scores
for (signature in names(single_gene_signatures)) {
  scores[, signature] <- norm_counts[signature, ]
  for (percentile in percentiles) {
    threshold_low <- quantile(scores[, signature], probs = percentile)
    threshold_high <- quantile(scores[, signature], probs = 1 - percentile)
    categorical_var <- ifelse(scores[, signature] <= threshold_low, "Low",
                              ifelse(scores[, signature] >= threshold_high, "High", "Intermediate"))
    col_name <- paste0(signature, "_", percentile * 100, "pct")
    colData(dds)[, col_name] <- categorical_var
  }
}

# create empty p-value table
surv_pval_mat <- matrix(
  NA,
  nrow = length(colnames(scores)),
  ncol = length(percentiles),
  dimnames = list(colnames(scores), paste0(percentiles * 100, "pct"))
)

# loop through signatures and percentiles to plot PCA and survival
for (signature in colnames(scores)) {
  for (percentile in percentiles) {
    # prevent the loop from stopping if an error occurs
    tryCatch({
      col_name <- paste0(signature, "_", percentile * 100, "pct")
      dds_filtered <- colData(dds)
      # obtain the indices of non-Intermediate patients
      pat_to_remove <- which(dds_filtered[, col_name] != "Intermediate")
      # remove intermediate patients
      dds_filtered <- dds_filtered[pat_to_remove, ]
      # replace NA time to death with last follow-up time for censored cases
      notDead <- is.na(dds_filtered$days_to_death)
      dds_filtered$days_to_death[notDead] <- dds_filtered$days_to_last_follow_up[notDead]
      # create event column (s = TRUE if dead, FALSE if alive)
      dds_filtered$s <- grepl("dead|deceased", dds_filtered$vital_status, ignore.case = TRUE)
      # create grouping factor
      dds_filtered$type <- factor(dds_filtered[[col_name]], levels = c( "Low", "High"))
      # keep only required columns
      dds_filtered_surv <- dds_filtered[, c("days_to_death", "s", "type")]
      # convert days to months to ease readability
      dds_filtered_surv$days_to_death <- dds_filtered_surv$days_to_death / 30.5
      # survival model
      surv_obj <- Surv(time = as.numeric(dds_filtered_surv$days_to_death), event = dds_filtered_surv$s)
      fit <- survfit(surv_obj ~ type, data = dds_filtered_surv)
      # extract p-value
      # create the Kaplan-Meier plot only if p-value < THRESHOLD
      # create the PCA only under the same condition
      pval <- surv_pvalue(fit)$pval
      # save the p-value in the matrix
      surv_pval_mat[signature, paste0(percentile * 100, "pct")] <- pval
      if (pval < THRESHOLD) {
        # survival legend title
        if (pval < 0.00001) {
          surv_title <- paste0("Expression groups (p-val < 0.0001)")
        } else {
          surv_title <- paste0("Expression groups (p-val = ", round(pval, 5), ")")
        }
        # custom legend labels with sample sizes
        label.add.n <- function(x) {
          n <- sum(dds_filtered_surv$type == x)
          paste0(x, " (n=", n, ")")
        }
        legend_labels <- sapply(levels(dds_filtered_surv$type), label.add.n)

        # plot survival
        p <- ggsurvplot(
          fit,
          data = dds_filtered_surv,
          theme = theme_survminer(),
          risk.table = FALSE,
          pval = FALSE,
          conf.int = FALSE,
          fontsize = 50,
          xlab = "Time (months)",
          ylab = "Survival probability",
          legend.title = surv_title,
          legend.labs = legend_labels,
          palette = c("#255BA8", "#ED412B"),
          break.time.by = 12
        )
        ggsave(filename = paste0("./screening/survival/", project, "/Survival_", signature, "_", percentile * 100, "pct.png"), plot = p$plot, dpi = DPI, width = 6, height = 4)
        # filter dds object for PCA
        dds_filtered_pca <- dds[, rownames(dds_filtered)]
        # calculate PC1 and PC2 for PCA plot
        pca_data <- plotPCA(vst(dds_filtered_pca), intgroup = col_name, returnData = TRUE)
        percentVar <- round(100 * attr(pca_data, "percentVar"), 2)
        # plot PCA
        ggplot(pca_data, aes(PC1, PC2, color = get(col_name))) +
        geom_point(size = 3) +
        scale_color_manual(values = c("Low" = "#255BA8", "High" = "#ED412B")) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(paste0("PCA - ", signature, " - ", percentile * 100, "pct")) +
        theme_minimal() +
        labs(color = "Expression groups")
        # save PCA as png
        ggsave(filename = paste0("./screening/PCA/", project, "/PCA_", signature, "_", percentile * 100, "pct.png"), dpi = DPI, width = 6 * FACTOR, height = 4 * FACTOR)
    }
    }, error = function(e) {

      # Print a warning but continue the loop
      message("⚠️ Error with signature '", signature,
              "', percentile ", percentile,
              ": ", conditionMessage(e))
      # store NA upon error
      surv_pval_mat[signature, paste0(percentile * 100, "pct")] <- NA
      # Continue to next iteration
      return(NULL)
    })
  }
}

# save scores matrix
write.table(
  scores,
  file = paste0("./screening/survival/", project,"/patient_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# save patient stratification on colData
tmp <- as.data.frame(colData(dds))
logical <- endsWith(colnames(tmp), "pct")
write.table(
  tmp[, logical],
  file = paste0("./screening/survival/", project,"/patient_scores_categorical.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# save survival p-value matrix
write.table(
  surv_pval_mat,
  file = paste0("./screening/survival/", project,"/survival_pval.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# create a filtered version of the p-value matrix
surv_pval_filtered <- surv_pval_mat
surv_pval_filtered[surv_pval_filtered >= THRESHOLD] <- ""

# save filtered p-value matrix
write.table(
  surv_pval_filtered,
  file = survival_table,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)