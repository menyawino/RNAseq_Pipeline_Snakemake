#!/usr/bin/env Rscript

# install required packages
install.packages("optparse")
install.packages("sleuth")
install.packages("tidyverse")
install.packages("gplots")
install.packages("biomaRt")

# install them using mamba command line (biomard it from biocmanager)



library(optparse)
library(sleuth)
library(tidyverse)
library(gplots)
library(biomaRt)

# Option parsing
option_list <- list(
  make_option(c("--output"), type = "character", help = "Output directory for results"),
  make_option(c("--log"), type = "character", help = "Log file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Helper functions
plotVolcano <- function(soLRT, pvalue_cutoff = 0.05, fold_cutoff = 0, include_inf = TRUE, level = 'Transcript', build = 'Hg38') {
  countsPerReplicate <- kallisto_table(soLRT, use_filtered = TRUE, normalized = TRUE)
  countsPerCondition <- countsPerReplicate %>%
    group_by(target_id, condition) %>%
    summarize(tpm = mean(tpm)) %>%
    pivot_wider(names_from = 'condition', values_from = 'tpm') %>%
    mutate(logFoldChange = log2(HCM / Control))
  
  resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
  countsPerCondition$qval <- resultsTable$qval[match(countsPerCondition$target_id, resultsTable$target_id)]
  
  if (!include_inf) {
    countsPerCondition <- subset(countsPerCondition, logFoldChange != 'Inf' & logFoldChange != '-Inf')
  }
  
  countsPerCondition$log_qvalues <- -log10(countsPerCondition$qval)
  countsPerCondition$direction <- ifelse(countsPerCondition$qval <= pvalue_cutoff & countsPerCondition$logFoldChange >= fold_cutoff, "upregulated",
                                         ifelse(countsPerCondition$qval <= pvalue_cutoff & countsPerCondition$logFoldChange <= -fold_cutoff, "downregulated", "not significant"))
  
  ggplot(countsPerCondition, aes(x = logFoldChange, y = log_qvalues)) +
    geom_point(aes(color = direction), size = 2) +
    scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "not significant" = "grey")) +
    theme_minimal() +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = paste(level, "Volcano Plot at Pval", pvalue_cutoff, "& LFC", fold_cutoff, "for", build)) +
    geom_vline(xintercept = c(-fold_cutoff, fold_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plotPCA <- function(countsPerReplicate, x_ax = 'PC1', y_ax = 'PC2') {
  tpm_centered <- t(countsPerReplicate - rowMeans(countsPerReplicate))
  tpm_prcomp <- prcomp(tpm_centered)
  
  sample_names <- colnames(countsPerReplicate)
  conditions <- c('HCM', 'HCM', 'HCM', 'HCM', 'HCM', 'HCM', 'HCM', 'HCM', 'HCM', 'Control', 'Control', 'Control', 'HCM', 'HCM', 'HCM')
  plot_df <- data.frame(PC1 = tpm_prcomp$x[, 1], PC2 = tpm_prcomp$x[, 2], PC3 = tpm_prcomp$x[, 3], Condition = conditions)
  
  ggplot(plot_df, aes(x = plot_df[, x_ax], y = plot_df[, y_ax], col = Condition)) +
    geom_point() +
    theme_minimal() +
    labs(title = "PCA Plot", x = x_ax, y = y_ax)
}

fetchTPMmatrix <- function(soLRT, pvalue_cutoff = 0.05, fold_cutoff = 0) {
  countsPerReplicate <- kallisto_table(soLRT, use_filtered = TRUE, normalized = TRUE)
  countsPerCondition <- countsPerReplicate %>%
    group_by(target_id, condition) %>%
    summarize(tpm = mean(tpm)) %>%
    pivot_wider(names_from = 'condition', values_from = 'tpm') %>%
    mutate(logFoldChange = log2(HCM / Control))
  
  resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
  sigResultsTable <- dplyr::filter(resultsTable, qval <= pvalue_cutoff)
  sigResultsList <- sigResultsTable$target_id
  
  countsPerCondition <- countsPerCondition[countsPerCondition$target_id %in% sigResultsList, ]
  countsPerCondition <- subset(countsPerCondition, abs(logFoldChange) >= fold_cutoff)
  
  sigResultsList <- countsPerCondition$target_id
  countsPerReplicate <- countsPerReplicate[countsPerReplicate$target_id %in% sigResultsList, ]
  
  countsPerReplicate <- countsPerReplicate %>%
    pivot_wider(names_from = 'sample', values_from = 'tpm') %>%
    as.data.frame()
  
  rownames(countsPerReplicate) <- countsPerReplicate$target_id
  countsPerReplicate <- countsPerReplicate[, -1]
  
  return(countsPerReplicate)
}

convertToLogTPMmatrix <- function(countsMatrix, include_zeros = TRUE) {
  logCounts <- log10(countsMatrix + 0.0001)
  
  if (!include_zeros) {
    logCounts[logCounts == -4] <- NA
  }
  
  return(logCounts)
}

prepOutputDF <- function(soLRT) {
  resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
  sigResultsTable <- dplyr::filter(resultsTable, qval <= 0.05)
  sigResultsList <- sigResultsTable$target_id
  
  countsPerCondition <- fetchConditionTPMmatrix(soLRT)
  countsPerCondition <- countsPerCondition[sigResultsList, ]
  countsPerCondition <- countsPerCondition %>%
    mutate(logFoldChange = log2(HCM / Control))
  
  countsPerCondition$qval <- resultsTable$qval[match(rownames(countsPerCondition), resultsTable$target_id)]
  return(countsPerCondition)
}

fetchConditionTPMmatrix <- function(soLRT) {
  countsPerReplicate <- kallisto_table(soLRT, use_filtered = TRUE, normalized = TRUE)
  countsPerCondition <- countsPerReplicate %>%
    group_by(target_id, condition) %>%
    summarize(tpm = mean(tpm)) %>%
    pivot_wider(names_from = 'condition', values_from = 'tpm') %>%
    as.data.frame()
  
  rownames(countsPerCondition) <- countsPerCondition$target_id
  countsPerCondition <- countsPerCondition[, -1]
  
  return(countsPerCondition)
}

# Main script
build <- 'Hg38'
level <- 'Gene'

# get sample_id from the kallisto output directory
sample_id <- dir(file.path("/mnt/omar/pipelines/rnaseq_snakemake/analysis/006_count/kallisto/"))

# Load metadata
sample_metadata <- read.csv("/mnt/omar/pipelines/rnaseq_snakemake/samples_r.csv")

# Replace sample with sample_id in metadata to match sample and lane names after checking that the order of sample_id is the same as sample in the metadata by checking the first chunk before the _ in the sample_id
sample_ids <- gsub("_.*", "", sample_id)

# sample_ids <- unique(sample_ids)

# match sample_id with sample in metadata
sample_metadata$sample <- sample_metadata$sample[match(sample_id, sample_metadata$sample)]

kal_dirs <- file.path("/mnt/omar/pipelines/rnaseq_snakemake/analysis/006_count/kallisto", sample_id)

# Add path to metadata
sample_metadata <- dplyr::mutate(sample_metadata, path = kal_dirs)


# Transcript-level DEA
so <- sleuth_prep(sample_metadata, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
soLRT <- sleuth_lrt(so, 'reduced', 'full')

# HG37 Gene-level DEA
mart_hg37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://grch37.ensembl.org')
t2g_hg37 <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart_hg37)
t2g_hg37 <- dplyr::rename(t2g_hg37, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so_hg37 <- sleuth_prep(sample_metadata, extra_bootstrap_summary = TRUE, target_mapping = t2g_hg37, aggregation_column = 'ens_gene', gene_mode = TRUE)
so_hg37 <- sleuth_fit(so_hg37, ~condition, 'full')
so_hg37 <- sleuth_fit(so_hg37, ~1, 'reduced')
soLRT_hg37 <- sleuth_lrt(so_hg37, 'reduced', 'full')

# HG38 Gene-level DEA
t2g_hg38 <- read.csv("kallisto_sleuth/LRS_kallisto_hg38_t2g2.csv", header = TRUE, sep = ",")
so_hg38 <- sleuth_prep(sample_metadata, extra_bootstrap_summary = TRUE, target_mapping = t2g_hg38, aggregation_column = 'ens_gene', gene_mode = TRUE)
so_hg38 <- sleuth_fit(so_hg38, ~condition, 'full')
so_hg38 <- sleuth_fit(so_hg38, ~1, 'reduced')
soLRT_hg38 <- sleuth_lrt(so_hg38, 'reduced', 'full')

# Output results
countsPerCondition <- prepOutputDF(soLRT)
sigResultsList <- rownames(countsPerCondition)

# For hg37 transcripts
mart_hg37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://grch37.ensembl.org')
t2g_hg37 <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart_hg37)
t2g_hg37 <- dplyr::rename(t2g_hg37, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2g_hg37 <- t2g_hg37[t2g_hg37$target_id %in% sigResultsList, ]
countsPerCondition$target_id <- rownames(countsPerCondition)
countsPerCondition <- merge(countsPerCondition, t2g_hg37[, c("target_id", "ext_gene")], by = "target_id", all.x = TRUE)
countsPerCondition <- countsPerCondition %>% rename('gene_name' = 'ext_gene')
rownames(countsPerCondition) <- countsPerCondition$target_id
countsPerCondition <- countsPerCondition[, -1]
write.table(countsPerCondition, paste0(opt$output, 'LRS_kallisto_', build, '_sig', level, '_0.05_foldchange.csv'), sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)

# For hg38 transcripts
countsPerCondition$gene_name <- sapply(strsplit(as.character(sigResultsList), "\\|"), `[`, 6)
write.table(countsPerCondition, paste0(opt$output, 'LRS_kallisto_', build, '_sig', level, '_0.05_foldchange.csv'), sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)

# For hg37 genes
t2g_hg37 <- t2g_hg37[, c('ens_gene', 'ext_gene')]
t2g_hg37 <- t2g_hg37[t2g_hg37$ens_gene %in% sigResultsList, ]
t2g_hg37 <- unique(t2g_hg37)
t2g_hg37 <- t2g_hg37 %>% rename('target_id' = 'ens_gene')
countsPerCondition$target_id <- rownames(countsPerCondition)
countsPerCondition <- merge(countsPerCondition, t2g_hg37[, c("target_id", "ext_gene")], by = "target_id", all.x = TRUE)
countsPerCondition <- countsPerCondition %>% rename('gene_name' = 'ext_gene')
rownames(countsPerCondition) <- countsPerCondition$target_id
countsPerCondition <- countsPerCondition[, -1]
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







