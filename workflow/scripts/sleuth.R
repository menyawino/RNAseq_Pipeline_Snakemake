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
sample_metadata <- read.csv("samples_r_metadata.csv")

# Replace sample with sample_id in metadata to match sample and lane names after checking that the order of sample_id is the same as sample in the metadata by checking the first chunk before the _ in the sample_id
sample_ids <- gsub("_.*", "", sample_id)
sample_metadata$sample <- sample_metadata$sample[match(sample_ids, sample_metadata$sample)]

kal_dirs <- file.path("analysis/006_count/kallisto", sample_id)

# Add path to metadata
sample_metadata



#Prepare sleuth object
sample_metadata <- dplyr::mutate(sample_metadata, path = kal_dirs)


so <- sleuth_prep(sample_metadata, ~condition)

# Fit models
so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, "reduced")

# Perform likelihood ratio test
so <- sleuth_lrt(so, "reduced", "full")

# Get results
sleuth_results <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)

# Write results to file
write.table(sleuth_results, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

# Get significant results
sleuth_significant <- dplyr::filter(sleuth_results, pval <= 0.05)

# Remove things that are not genes from target_id pattern ENST00000700062.1| "remove anything after the first ."
sleuth_significant$target_id <- gsub("\\..*", "", sleuth_significant$target_id)
head(sleuth_significant, 10)


# write significant results to file
write.table(sleuth_significant, file = "significant_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Generate a PDF report with visualizations
pdf(opt$report, width=11, height=8.5)

# number of significant results
nrow(sleuth_significant)

# plot p-value distribution
hist(sleuth_results$pval, breaks = 50, main = "P-Value Distribution")

# Plot PCA with group the points with a circle around them
pca_data <- plot_pca(so, color_by = "condition")
print(pca_data)

# End PDF
dev.off()

# Reset output
sink()
