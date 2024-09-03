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

# Redirect output to log file
sink(opt$log, append=TRUE, split=TRUE)




sample_id <- dir(file.path("analysis/006_count/kallisto/"))

# drop the first and last elements of the sample_id vector
# sample_id <- sample_id[-c(1, length(sample_id))]

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
write.table(sleuth_results, file = "dongol.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Generate a PDF report with visualizations
pdf("dongol.pdf", width=11, height=8.5)


# Plot PCA
plot_pca <- plot_pca(so, color_by = "condition")
print(plot_pca)


# # Get significant results
# sleuth_significant <- dplyr::filter(sleuth_results, qval <= 0.05)

# # plot bootstrap of the top 10 most significant genes
# plot_bootstrap(so, sleuth_significant$target_id[1:10], units = "est_counts", color_by = "condition")


# Plot for specific genes (add any genes of interest here)
genes_of_interest <- c("GENE1", "GENE2") # replace with actual gene names
for (gene in genes_of_interest) {
    plot_gene <- plot_bootstrap(so, gene_id = gene, units = "est_counts")
    print(plot_gene)
}

# End PDF
dev.off()

# Reset output
sink()
