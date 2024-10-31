# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install("devtools")    # only if devtools not yet installed
# BiocManager::install("pachterlab/sleuth")
library(sleuth)

# install.packages("tidyverse")
library(tidyverse)

# install.packages('gplots')
library(gplots)

##########################################################################################################################################
# Helper functions

plotVolcano <- function(pvalue_cutoff = 0.05, fold_cutoff = 0, include_inf = TRUE,level='Transcript',build='Hg38') {
  countsPerReplicate <- kallisto_table(soLRT, use_filtered= T, normalized = T)
  countsPerCondition <- countsPerReplicate %>%
    group_by(target_id, condition) %>%
    summarize(tpm = mean(tpm))
  countsPerCondition <- countsPerCondition %>%
    pivot_wider(names_from = 'condition', values_from = 'tpm')
  countsPerCondition <- countsPerCondition %>%
    mutate(logFoldChange = log2(HCM / Control))
  
  resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
  countsPerCondition$qval <- resultsTable$qval[match(countsPerCondition$target_id, resultsTable$target_id)]
  
  if(include_inf == FALSE){
    countsPerCondition <- subset(countsPerCondition, logFoldChange != 'Inf' & logFoldChange != '-Inf')
  }
  
  countsPerCondition$log_qvalues<--log(countsPerCondition$qval,10)
  countsPerCondition$direction <- ifelse(countsPerCondition$qval <= pvalue_cutoff & countsPerCondition$logFoldChange >= fold_cutoff, "upregulated",
                                         ifelse(countsPerCondition$qval <= pvalue_cutoff & countsPerCondition$logFoldChange <= -fold_cutoff, "downregulated", "not significant"))
  upregulated<-subset(countsPerCondition,direction=='upregulated')
  print(paste('Upregulated=', nrow(upregulated)))
  downregulated<-subset(countsPerCondition,direction=='downregulated')
  print(paste('Downregulated=', nrow(downregulated)))
  
  plotTitle<-paste(level,"Volcano Plot at Pval",pvalue_cutoff, "& LFC", fold_cutoff, "for ",build)
  
  ggplot(countsPerCondition, aes(x = logFoldChange, y = log_qvalues)) +
    geom_point(aes(color = direction), 
               size = 2) +  # Dynamic color based on conditions
    scale_color_manual(values = c("upregulated" = "red", 
                                  "downregulated" = "blue", 
                                  "not significant" = "grey")) +
    theme_minimal() +
    labs(x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-value", 
         title = plotTitle) +
    geom_vline(xintercept = c(-fold_cutoff, fold_cutoff), linetype = "dashed", color = "black") +  # Fold change cutoff
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "black") +  # p-value cutoff
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

plotPCA <- function(countsPerReplicate, x_ax = 'PC1', y_ax = 'PC2') {
  tpm_centered <- t(countsPerReplicate-rowMeans(countsPerReplicate))
  tpm_prcomp <- prcomp(tpm_centered)
  
  sample_names <- colnames(countsPerReplicate)
  conditions <- c('HCM','HCM','HCM','HCM','HCM','HCM','HCM','HCM','HCM','Control','Control','Control','HCM','HCM','HCM')
  plot_df <- data.frame(PC1 = tpm_prcomp$x[,1], PC2 = tpm_prcomp$x[,2], PC3 = tpm_prcomp$x[,3], Condition = conditions)
  
  ggplot(plot_df, aes(x = plot_df[,x_ax], y = plot_df[,y_ax], col = Condition)) +
    geom_point() +
    theme_minimal() +
    labs(title = "PCA Plot", x = x_ax, y = y_ax)
}


fetchTPMmatrix <- function(pvalue_cutoff = 0.05, fold_cutoff = 0) {
  countsPerReplicate <- kallisto_table(soLRT, use_filtered= T, normalized = T)
  countsPerCondition <- countsPerReplicate %>%
    group_by(target_id, condition) %>%
    summarize(tpm = mean(tpm))
  
  countsPerCondition <- countsPerCondition %>%
    pivot_wider(names_from = 'condition', values_from = 'tpm')
  
  countsPerCondition <- countsPerCondition %>%
    mutate(logFoldChange = log2(HCM / Control))
  
  resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
  sigResultsTable <- dplyr::filter(resultsTable, qval <= pvalue_cutoff)
  sigResultsList <- sigResultsTable$target_id
  
  countsPerCondition <- countsPerCondition[countsPerCondition$target_id %in% sigResultsList, ]
  countsPerCondition <- subset(countsPerCondition, abs(logFoldChange) >= fold_cutoff)
  
  sigResultsList <- countsPerCondition$target_id
  
  countsPerReplicate <- countsPerReplicate[countsPerReplicate$target_id %in% sigResultsList, ]
  
  countsPerReplicate <- countsPerReplicate[,c('target_id','sample','tpm')]
  countsPerReplicate <- countsPerReplicate %>%
    pivot_wider(names_from = 'sample', values_from = 'tpm')
  
  countsPerReplicate <- as.data.frame(countsPerReplicate)
  rownames(countsPerReplicate) <- countsPerReplicate$target_id
  countsPerReplicate <- countsPerReplicate[,-1]
  
  return(countsPerReplicate)
  
}

fetchConditionTPMmatrix <- function() {
  countsPerReplicate <- kallisto_table(soLRT, use_filtered= T, normalized = T)
  countsPerCondition <- countsPerReplicate %>%
    group_by(target_id, condition) %>%
    summarize(tpm = mean(tpm))
  countsPerCondition <- countsPerCondition %>%
    pivot_wider(names_from = 'condition', values_from = 'tpm')
  
  countsPerCondition <- as.data.frame(countsPerCondition)
  rownames(countsPerCondition) <- countsPerCondition$target_id
  countsPerCondition <- countsPerCondition[,-1]
  
  return(countsPerCondition)
  
}

convertToLogTPMmatrix <- function(countsMatrix, include_zeros = TRUE) {
  logCounts <- log10(countsMatrix+0.0001)
  
  if (include_zeros == FALSE){
    logCounts[logCounts == -4] <- NA
  }
  
  return(logCounts)
  
}

prepOutputDF <- function(){
  resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
  sigResultsTable <- dplyr::filter(resultsTable, qval <= 0.05)
  sigResultsList <- sigResultsTable$target_id
  
  countsPerCondition <- fetchConditionTPMmatrix()
  countsPerCondition <- countsPerCondition[sigResultsList, ]
  countsPerCondition <- countsPerCondition %>%
    mutate(logFoldChange = log2(HCM / Control))
  
  countsPerCondition$qval <- resultsTable$qval[match(rownames(countsPerCondition), resultsTable$target_id)]
  return(countsPerCondition)
}

##########################################################################################################################################
# Preparing the dataframe of samples, conditions and input file paths
build='Hg38'
level='Gene'

kallistoOutputs <- dir(file.path("kallisto_sleuth",paste0("kallisto_",build)))

kallistoDirs <- file.path("kallisto_sleuth",paste0("kallisto_",build), kallistoOutputs)

s2c <- read.table(file.path("kallisto_sleuth", "metadata.csv"), sep=",", header = TRUE, stringsAsFactors=FALSE)

s2c <- dplyr::mutate(s2c, path = kallistoDirs)

##########################################################################################################################################
# Transcript-level DEA

# Running LRT
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

soLRT <- sleuth_lrt(so, 'reduced', 'full')

##########################################################################################################################################
# HG37 Gene-level DEA

# Load the biomaRt package to access the BioMart database
mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'https://grch37.ensembl.org'  # Use the GRCh37 archive
)


# Retrieve transcript-to-gene mappings
t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id","ensembl_gene_id", "external_gene_name"),
  mart = mart
)

# Rename columns for clarity and compatibility with downstream analyses
t2g <- dplyr::rename(t2g,
                     target_id = ensembl_transcript_id,  # Rename 'ensembl_transcript_id' to 'target_id'
                     ens_gene = ensembl_gene_id,         # Rename 'ensembl_gene_id' to 'ens_gene'
                     ext_gene = external_gene_name       # Rename 'external_gene_name' to 'ext_gene'
)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = t2g, aggregation_column = 'ens_gene',gene_mode=TRUE)

so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

soLRT <- sleuth_lrt(so, 'reduced', 'full')

##########################################################################################################################################
# HG38 Gene-level DEA

t2g <- read.csv("kallisto_sleuth/LRS_kallisto_hg38_t2g2.csv", header = TRUE, sep = ",")

# Running LRT
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = t2g, aggregation_column = 'ens_gene',gene_mode=TRUE)

so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

soLRT <- sleuth_lrt(so, 'reduced', 'full')


##########################################################################################################################################
# First look at the data

resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
sigResultsTable <- dplyr::filter(resultsTable, qval <= 0.05)
sigResultsList <- sigResultsTable$target_id

countsPerCondition <- fetchConditionTPMmatrix()
countsPerCondition <- countsPerCondition[sigResultsList, ]
countsPerCondition <- countsPerCondition %>%
  mutate(logFoldChange = log2(HCM / Control))
countsPerCondition_noInf <- subset(countsPerCondition, logFoldChange != 'Inf' & logFoldChange != '-Inf')

# Histogram
breaks <- seq(floor(min(countsPerCondition_noInf$logFoldChange)), ceiling(max(countsPerCondition_noInf$logFoldChange)), by = 1)
xlim_min <- floor(min(countsPerCondition_noInf$logFoldChange))
xlim_max <- ceiling(max(countsPerCondition_noInf$logFoldChange))
hist(countsPerCondition_noInf$logFoldChange,xlab='LFC distribution', col = "lightblue",main = paste('Histogram of',level,'LFC Distributions For',build),breaks=breaks)

##########################################################################################################################################
# Plots

# Scatter plot
countsPerCondition <- fetchConditionTPMmatrix()
logCountsPerCondition <- convertToLogTPMmatrix(countsPerCondition, include_zeros=TRUE)

Control <- logCountsPerCondition[,c('Control')]
HCM <- logCountsPerCondition[,c('HCM')]

plot(HCM ~ Control, 
     xlim = c(floor(min(Control,HCM)),ceiling(max(Control))), 
     ylim = c(floor(min(Control,HCM)),ceiling(max(HCM))), main=paste(level,'Scatter Plot for', build))
abline(0, 1, col = "red")

# Boxplot per replicate
countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 1, fold_cutoff= 0)
logCountsPerReplicate <- convertToLogTPMmatrix(countsPerReplicate, include_zeros=TRUE)

controls <- c('Myo182', 'Myo96', 'Myoc169')
conditions <- c('23011','48977','50940MF','53311','56091Myo','56253Myo','57050Myo','58862Myo','60773MF','myo49610','myo59049','myo63988')
samples <- c(controls, conditions)
colors = c(rep("grey", 3), rep("lightblue", 12))
boxplot(logCountsPerReplicate,col = colors,names=samples,las=2,main=paste(level, "Boxplots Including TPM = 0 for", build))

logCountsPerReplicate <- convertToLogTPMmatrix(countsPerReplicate, include_zeros=FALSE)

boxplot(logCountsPerReplicate,col = colors,names=samples,las=2,main=paste(level, "Boxplots Excluding TPM = 0 for", build))

# Boxplot per condition
countsPerCondition <- fetchConditionTPMmatrix()
logCountsPerCondition <- convertToLogTPMmatrix(countsPerCondition, include_zeros=TRUE)

boxplot(logCountsPerCondition,col = c('gray','lightblue'),names=c('Control','HCM'),las=2,main=paste(level, "Boxplots Including TPM = 0 for", build))

logCountsPerCondition <- convertToLogTPMmatrix(countsPerCondition, include_zeros=FALSE)

boxplot(logCountsPerCondition,col = c('gray','lightblue'),names=c('Control','HCM'),las=2,main=paste(level, "Boxplots Excluding TPM = 0 for", build))


# Volcano plot
pvalue_cutoff=0.05
fold_cutoff=2
plotVolcano(pvalue_cutoff = pvalue_cutoff, fold_cutoff = fold_cutoff,include_inf=TRUE,level=level,build=build)
plotVolcano(pvalue_cutoff = pvalue_cutoff, fold_cutoff = fold_cutoff,include_inf=FALSE,level=level,build=build)

# PCA
# https://support.bioconductor.org/p/102053/

countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 1, fold_cutoff = 0)
plotPCA(countsPerReplicate,'PC1','PC2')
plotPCA(countsPerReplicate,'PC1','PC3')
plotPCA(countsPerReplicate,'PC2','PC3')

countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 0.05, fold_cutoff = 0)
plotPCA(countsPerReplicate,'PC1','PC2')
plotPCA(countsPerReplicate,'PC1','PC3')
plotPCA(countsPerReplicate,'PC2','PC3')

countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 0.05, fold_cutoff = 2)
plotPCA(countsPerReplicate,'PC1','PC2')
plotPCA(countsPerReplicate,'PC1','PC3')
plotPCA(countsPerReplicate,'PC2','PC3')

# Dendrogram before filtration
countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 1, fold_cutoff = 0)
plot(hclust(as.dist(1 - cor(countsPerReplicate))), main=paste(level,'Dendrogram For', build))


# Dendrogram after filtration - not taking fold change into account
countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 0.05, fold_cutoff = 0)
plot(hclust(as.dist(1 - cor(countsPerReplicate))), main=paste(level,'Dendrogram at Pval', pvalue_cutoff ,'For', build))

# Corresponding heatmap
colv = as.dendrogram(hclust(as.dist(1-cor(countsPerReplicate))))
countsPerReplicate <- as.matrix(countsPerReplicate)
heatmap.2(x = countsPerReplicate,Colv = colv,
          col = colorRampPalette(c("darkblue","white","darkred"))(100),
          trace = "none",         # Remove row/column labels
          scale = 'row',
          main = paste(level,'Heatmap at Pval', pvalue_cutoff ,'For', build),
          key = TRUE,
          symkey = FALSE, 
          keysize = 1.5,
          density.info = 'none')


# Dendrogram after filtration - taking fold change into account
countsPerReplicate <- fetchTPMmatrix(pvalue_cutoff = 0.05, fold_cutoff = 2)
plot(hclust(as.dist(1 - cor(countsPerReplicate))), main=paste(level,'Dendrogram at Pval', pvalue_cutoff , '& LFC', fold_cutoff,'For', build))

# Corresponding heatmap
colv = as.dendrogram(hclust(as.dist(1-cor(countsPerReplicate))))
countsPerReplicate <- as.matrix(countsPerReplicate)
heatmap.2(x = countsPerReplicate,Colv = colv,
          col = colorRampPalette(c("darkblue","white","darkred"))(100),  # Color scale          scale = 'row',
          trace = "none",         # Remove row/column labels
          scale = 'row',
          main =paste(level,'Dendrogram at Pval', pvalue_cutoff , '& LFC', fold_cutoff,'For', build),
          cexRow = 0.8,           # Row label size
          cexCol = 0.8,           # Column label size
          key = TRUE,             # Show color scale legend
          symkey = FALSE,          # Use asymmetric color scale legend
          keysize = 1.5,
          density.info = 'none')

##########################################################################################################################################
# Producing output files


# For hg37 transcripts
countsPerCondition <- prepOutputDF()
sigResultsList <- rownames(countsPerCondition)

mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'https://grch37.ensembl.org'  # Use the GRCh37 archive
)


t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id","ensembl_gene_id", "external_gene_name"),
  mart = mart
)

t2g <- dplyr::rename(t2g,
                     target_id = ensembl_transcript_id,  # Rename 'ensembl_transcript_id' to 'target_id'
                     ens_gene = ensembl_gene_id,         # Rename 'ensembl_gene_id' to 'ens_gene'
                     ext_gene = external_gene_name       # Rename 'external_gene_name' to 'ext_gene'
)

t2g<-t2g[t2g$target_id %in% sigResultsList, ]
countsPerCondition$target_id <- rownames(countsPerCondition)

countsPerCondition <- merge(countsPerCondition, t2g[, c("target_id", "ext_gene")], 
                            by = "target_id", all.x = TRUE)

countsPerCondition <- countsPerCondition %>%
  rename('gene_name' = 'ext_gene')

rownames(countsPerCondition)<-countsPerCondition$target_id
countsPerCondition<-countsPerCondition[,-1]

write.table(countsPerCondition, paste0('kallisto_sleuth/LRS_kallisto_',build,'_sig',level,'_0.05_foldchange.csv'), sep = ',', row.names = T, col.names=NA, quote = F)


# For hg38 transcripts
countsPerCondition <- prepOutputDF()
sigResultsList <- rownames(countsPerCondition)

countsPerCondition$gene_name <- sapply(strsplit(as.character(sigResultsList), "\\|"), `[`, 6)
write.table(countsPerCondition, paste0('kallisto_sleuth/LRS_kallisto_',build,'_sig',level,'_0.05_foldchange.csv'), sep = ',', row.names = T, col.names=NA, quote = F)


# For hg37 genes
countsPerCondition <- prepOutputDF()
sigResultsList <- rownames(countsPerCondition)

t2g<-t2g[,c('ens_gene','ext_gene')]
t2g<-t2g[t2g$ens_gene %in% sigResultsList, ]
t2g<-unique(t2g)
t2g <- t2g %>%
  rename('target_id' = 'ens_gene')
countsPerCondition$target_id <- rownames(countsPerCondition)

countsPerCondition <- merge(countsPerCondition, t2g[, c("target_id", "ext_gene")], 
                            by = "target_id", all.x = TRUE)

countsPerCondition <- countsPerCondition %>%
  rename('gene_name' = 'ext_gene')

rownames(countsPerCondition)<-countsPerCondition$target_id
countsPerCondition<-countsPerCondition[,-1]

write.table(countsPerCondition, paste0('kallisto_sleuth/LRS_kallisto_',build,'_sig',level,'_0.05_foldchange.csv'), sep = ',', row.names = T, col.names=NA, quote = F)

# For hg38 genes
resultsTable <- sleuth_results(soLRT, 'reduced:full', 'lrt', show_all = FALSE)
sigResultsTable <- dplyr::filter(resultsTable, qval <= 0.05)
sigResultsTable<-sigResultsTable[,-1]
sigResultsTable <- sigResultsTable %>% distinct()
sigResultsList <- sigResultsTable$target_id

countsPerCondition <- fetchConditionTPMmatrix()
countsPerCondition <- countsPerCondition[sigResultsList, ]
countsPerCondition <- countsPerCondition %>%
  mutate(logFoldChange = log2(HCM / Control))

countsPerCondition$qval <- resultsTable$qval[match(rownames(countsPerCondition), resultsTable$target_id)]

t2g<-t2g[,c('ens_gene','ext_gene')]
t2g<-t2g[t2g$ens_gene %in% sigResultsList, ]
t2g<-unique(t2g)
t2g <- t2g %>%
  rename('target_id' = 'ens_gene')
countsPerCondition$target_id <- rownames(countsPerCondition)

countsPerCondition <- merge(countsPerCondition, t2g[, c("target_id", "ext_gene")], 
                            by = "target_id", all.x = TRUE)

countsPerCondition <- countsPerCondition %>%
  rename('gene_name' = 'ext_gene')

rownames(countsPerCondition)<-countsPerCondition$target_id
countsPerCondition<-countsPerCondition[,-1]


write.table(countsPerCondition,  paste0('kallisto_sleuth/LRS_kallisto_',build,'_sig',level,'_0.05_foldchange.csv'), sep = ',', row.names = T,col.names=NA, quote = F)

#####################################################################################################