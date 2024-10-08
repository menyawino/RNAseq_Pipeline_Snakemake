#!/usr/bin/env Rscript

library(optparse)
library(sleuth)
library(ggplot2)
library(dplyr)

# Option parsing
option_list <- list(
  make_option(c("--output"), type = "character", help = "Output file path for differential expression results"),
  make_option(c("--report"), type = "character", help = "Output file path for the PDF report"),
  make_option(c("--log"), type = "character", help = "Log file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# get sample_id from the kallisto output directory
sample_id <- dir(file.path("analysis/006_count/kallisto/"))

# Load metadata
sample_metadata <- read.csv("samples_r.csv")

# Replace sample with sample_id in metadata to match sample and lane names after checking that the order of sample_id is the same as sample in the metadata by checking the first chunk before the _ in the sample_id
sample_ids <- gsub("_.*", "", sample_id)

# sample_ids <- unique(sample_ids)

# match sample_id with sample in metadata
sample_metadata$sample <- sample_metadata$sample[match(sample_id, sample_metadata$sample)]

kal_dirs <- file.path("analysis/006_count/kallisto", sample_id)

# Add path to metadata
sample_metadata <- dplyr::mutate(sample_metadata, path = kal_dirs)

sample_metadata


#Prepare sleuth object
so <- sleuth_prep(sample_metadata, ~condition)


# Fit models for statistical testing
so <- sleuth_fit(so, ~condition, "full")


# Perform Wald test for the 'condition' factor (specify the level you want to test)
so <- sleuth_wt(so, "conditionHCM")

# Get Wald test results for the specified condition
sleuth_results <- sleuth_results(so, "conditionHCM", show_all = TRUE)  # The second argument refers to the specific coefficient being tested

# Calculate fold change (FC) from the beta values (log2 fold change)
sleuth_results$fold_change <- 2^sleuth_results$b

# Write results to file (manual)
# write.table(sleuth_results, file = "sleuth_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Write results to file
write.table(sleuth_results, file = paste0(opt$output, "differentially_expressed_genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Get significant results
sleuth_significant <- sleuth_results %>%
  filter((fold_change >= 2 | fold_change <= 0.5) & pval < 0.05)

# Add log2 fold change to the significant results
sleuth_significant$log2_fold_change <- log2(sleuth_significant$fold_change)

# Remove things that are not genes from target_id pattern
sleuth_significant$target_id <- gsub("\\..*", "", sleuth_significant$target_id)

# Sort significant results by log2 fold change
sleuth_significant <- arrange(sleuth_significant, desc(log2_fold_change))

# Print significant results
head(sleuth_significant, 10)
tail(sleuth_significant, 10)


# write significant results to file
write.table(sleuth_significant, file = paste0(opt$output, "significant_genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


# Start PDF
pdf(paste0(opt$report, "report.pdf"))

# plot p-value distribution
hist(sleuth_results$pval, breaks = 50, main = "P-Value Distribution")

# Plot log2 fold change distribution
hist(sleuth_results$log2_fold_change, breaks = 50, main = "Log2 Fold Change Distribution")

# Plot MA plot
plot(so, "conditionHCM")


# Plot PCA with group the points with a circle around them
pca_data <- plot_pca(so, color_by = "condition")

# Print PCA data
print(pca_data)

# save the plot to a file
ggsave(paste0(opt$output, "pca_plot.pdf"), pca_data)

# End PDF
dev.off()
