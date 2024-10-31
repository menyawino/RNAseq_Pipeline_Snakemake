library(sleuth)
library(tidyverse)
library(gplots)
library(biomaRt)

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
