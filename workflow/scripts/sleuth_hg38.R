#!/usr/bin/env Rscript

library(optparse)
library(dplyr)  # Ensure dplyr is loaded
source("/mnt/omar/pipelines/rnaseq_snakemake/workflow/scripts/common_functions.R")

# Option parsing
option_list <- list(
  make_option(c("--output"), type = "character", help = "Output directory for results"),
  make_option(c("--log"), type = "character", help = "Log file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Main script
build <- 'Hg38'
level <- 'Gene'

# get sample_id from the kallisto output directory
sample_id <- dir(file.path("/mnt/omar/pipelines/rnaseq_snakemake/analysis/006_count/kallisto/"))

# Load metadata
sample_metadata <- read.csv("/mnt/omar/pipelines/rnaseq_snakemake/samples_r.csv")

# Replace sample with sample_id in metadata to match sample and lane names
sample_ids <- gsub("_.*", "", sample_id)
sample_metadata$sample <- sample_metadata$sample[match(sample_id, sample_metadata$sample)]
kal_dirs <- file.path("/mnt/omar/pipelines/rnaseq_snakemake/analysis/006_count/kallisto", sample_id)
sample_metadata <- dplyr::mutate(sample_metadata, path = kal_dirs)

# Function to read abundance tsv and return a dataframe in the specified format
read_abundance_tsv <- function(file_path) {
  df <- read.table(file_path, header = TRUE, sep = "\t")
  df <- df %>%
    mutate(
      ens_transcript = sapply(strsplit(target_id, "\\|"), `[`, 1),
      ens_gene = sapply(strsplit(target_id, "\\|"), `[`, 2),
      ext_gene = sapply(strsplit(target_id, "\\|"), `[`, 6)
    ) %>%
    dplyr::select(target_id, ens_transcript, ens_gene, ext_gene)
  return(df)
}

# Transcript-level DEA
so <- sleuth_prep(sample_metadata, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
soLRT <- sleuth_lrt(so, 'reduced', 'full')

# HG38 Gene-level DEA

# Generate t2g_hg38 from Kallisto input CSV
kallisto_input_csv <- paste0("/mnt/omar/pipelines/rnaseq_snakemake/analysis/006_count/kallisto/", sample_metadata$sample[1], "/abundance.tsv")
t2g_hg38 <- read_abundance_tsv(kallisto_input_csv)
so_hg38 <- sleuth_prep(sample_metadata, extra_bootstrap_summary = TRUE, target_mapping = t2g_hg38, aggregation_column = 'ens_gene', gene_mode = TRUE)
so_hg38 <- sleuth_fit(so_hg38, ~condition, 'full')
so_hg38 <- sleuth_fit(so_hg38, ~1, 'reduced')
soLRT_hg38 <- sleuth_lrt(so_hg38, 'reduced', 'full')

# Output results
countsPerCondition <- prepOutputDF(soLRT)
sigResultsList <- rownames(countsPerCondition)

# For hg38 transcripts
countsPerCondition$gene_name <- sapply(strsplit(as.character(sigResultsList), "\\|"), `[`, 6)
write.table(countsPerCondition, paste0(opt$output, 'LRS_kallisto_', build, '_sig', level, '_0.05_foldchange.csv'), sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)

# For hg38 genes
resultsTable_hg38 <- sleuth_results(soLRT_hg38, 'reduced:full', 'lrt', show_all = FALSE)
sigResultsTable_hg38 <- dplyr::filter(resultsTable_hg38, qval <= 0.05)
sigResultsTable_hg38 <- sigResultsTable_hg38 %>% distinct()
sigResultsList_hg38 <- sigResultsTable_hg38$target_id
countsPerCondition <- fetchConditionTPMmatrix(soLRT_hg38)
countsPerCondition <- countsPerCondition[sigResultsList_hg38, ]
countsPerCondition <- countsPerCondition %>% mutate(logFoldChange = log2(HCM / Control))
countsPerCondition$qval <- resultsTable_hg38$qval[match(rownames(countsPerCondition), resultsTable_hg38$target_id)]
t2g_hg38 <- t2g_hg38[, c('ens_gene', 'ext_gene')]
t2g_hg38 <- t2g_hg38[t2g_hg38$ens_gene %in% sigResultsList_hg38, ]
t2g_hg38 <- unique(t2g_hg38)
t2g_hg38 <- t2g_hg38 %>% rename('target_id' = 'ens_gene')
countsPerCondition$target_id <- rownames(countsPerCondition)
countsPerCondition <- merge(countsPerCondition, t2g_hg38[, c("target_id", "ext_gene")], by = "target_id", all.x = TRUE)
countsPerCondition <- countsPerCondition %>% rename('gene_name' = 'ext_gene')
rownames(countsPerCondition) <- countsPerCondition$target_id
countsPerCondition <- countsPerCondition[, -1]
write.table(countsPerCondition, paste0(opt$output, 'LRS_kallisto_', build, '_sig', level, '_0.05_foldchange.csv'), sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)
