
# # install biocmanager
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# # install ballgown
# BiocManager::install("ballgown")
# BiocManager::install("biomaRt")
# BiocManager::install("genefilter")

# install.packages("RSkittleBrewer")
# install.packages("devtools")
# install.packages("gplots")
# install.packages("ggfortify")
# install.packages("tidyverse")
# install.packages("ggrepel")



library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(gplots)
library(ggfortify)
library(biomaRt)
library(tidyverse)
library(ggrepel)

# mamba install r-base r-ballgown r-dplyr r-devtools r-ggplot2 r-gplots r-ggfortify r-biomart r-tidyverse r-ggrepel

path = getwd()

#Read the data
pheno_data = read.csv(paste0(path, "/sampless.csv"))

#Create a ballgown object
bg = ballgown(dataDir = paste0(path, "/analysis/006_count"),
    samplePattern = "23011_S",
    pData=pheno_data)

bg_expr = texpr(bg, 'FPKM')

#Filter bg to include only transcripts with variance > 1
bg_filt = subset(bg,"genefilter::rowVars(texpr(bg)) > 0.1",
    genomesubset=TRUE)

#Differential expression analysis for the transcripts
results_transcripts = stattest(bg_filt,
    feature="transcript",
    covariate="group",
    adjustvars = c("sex", "age"),
    getFC=TRUE, 
    meas="FPKM")

#Differential expression analysis for the genes
results_genes = stattest(bg_filt,
    feature="gene",
    covariate="group",
    adjustvars = c("sex", "age"),
    getFC=TRUE,
    meas="FPKM")


#Create a data frame with transcript results
results_transcripts = data.frame(
    geneNames=ballgown::geneNames(bg_filt),
    geneIDs=ballgown::geneIDs(bg_filt),
    transcriptName=ballgown::transcriptNames(bg_filt),
    results_transcripts)

#Use biomaRt package to get the gene names and IDs
mart <- useMart(biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl")

transcripts <- results_transcripts$transcriptName

res <- getBM(attributes = c('ensembl_transcript_id', 
                        'ensembl_gene_id', 
                        'external_transcript_name',
                        'external_gene_name'),
    filters = 'ensembl_transcript_id', 
    values = transcripts,
    mart = mart)



#Add gene names and gene IDs to the transcripts results
results_transcripts <- inner_join(
    results_transcripts,
    res,
    by = c("transcriptName" = "ensembl_transcript_id"))

#Add gene names and gene IDs to the genes results
idsnnames <- results_transcripts %>% 
    distinct(ensembl_gene_id, 
    external_gene_name, 
    geneIDs)

results_genes <- inner_join(
    results_genes,
    idsnnames,
    by = c("id" = "geneIDs"))

#Sort the transcripts results by p-value
results_transcripts = arrange(results_transcripts, pval)

#Sort the genes results by p-value
results_genes = arrange(results_genes, pval)

#Add a column to the results data frame with the log2 fold change
results_genes[,"de"] = log2(results_genes[,"fc"])

write.csv(results_transcripts, 
    paste0(path, "transcript_results.csv"), 
    row.names=FALSE)

write.csv(results_genes, 
    paste0(path, "gene_results.csv"),
    row.names=FALSE)



##################################################
################## Visualization #################
##################################################

#Transfom the FPKM values to log2 scale and reorder samples
fpkm = texpr(bg_filt, meas="FPKM")
fpkm = log2(fpkm+1)

#Create a palette
tropical= c('darkorange', 
    'dodgerblue', "hotpink", 
    "limegreen", "yellow")
palette(tropical)


# Open a PDF file where we will save some plots
outfile = paste0(path, "de_plots.pdf")
pdf(file=outfile)


#Create a boxplot summary for the FPKM values for each sample
par(mar=c(5, 4, 5, 8), xpd=TRUE)

boxplot(fpkm,
    col=as.numeric(as.factor(pheno_data$population))+1,
    las=2,
    ylab='log2(FPKM+1)')

legend("topright",
    inset=c(-0.35, 0.1),
    c("HCM","Healthy"),
    fill=c(2,3),
    horiz=FALSE,
    cex=1)


which(ballgown::geneNames(bg_filt) == 'NPPA')

gene_order=156

#Compare expression of a gene for both conditions
boxplot(fpkm[gene_order,] ~ pheno_data$population,
    border=c(2,3),
    main=paste(ballgown::geneNames(bg_filt)[gene_order],' : ',
    ballgown::transcriptNames(bg_filt)[gene_order]),
    pch=19,
    xlab="Type",
    ylab='log2(FPKM+1)')

points(fpkm[gene_order,] ~ jitter(as.numeric(pheno_data$population)),
    col=as.numeric(pheno_data$population)+1,
    pch=19)


# #Create a density plot of the FPKM values
# for(i in 1:ncol(fpkm)){
#     dens = density(fpkm[,i])
#     lines(dens$x, dens$y, col=i)
# }


# #Create a scatter plot of the FPKM values
# plot(fpkm[1,] ~ fpkm[2,],
#     pch=19,
#     xlab='log2(FPKM+1) sample 1',
#     ylab='log2(FPKM+1) sample 3')



sample_ids = c('SRR24302632', 
    'SRR24302634', 'SRR24302635', 
    'SRR24302636', 'SRR24302637', 
    'SRR24302638', 'SRR24302640', 
    'SRR24302643', 'SRR24302644', 
    'SRR24302645', 'SRR24302646')

samples_hcm = c('SRR24302632', 'SRR24302634', 
    'SRR24302635', 'SRR24302636', 'SRR24302640', 
    'SRR24302643')

samples_healthy = c('SRR24302637', 'SRR24302638',
    'SRR24302644','SRR24302645', 'SRR24302646')

sample_id = 'SRR24302634'
gene_pos = 3221

#Plot the transcripts for a gene in a single sample
plotTranscripts(ballgown::geneIDs(bg)[gene_pos],
    bg,
    main=c(paste0('Gene ', 
    ballgown::geneNames(bg)[gene_pos], 
    ' in sample '),
    sample=sample_id))

#Plot the transcripts for a gene in all samples
plotTranscripts(ballgown::geneIDs(bg)[gene_pos], 
    bg,
    main=c(paste0('Gene ', 
    ballgown::geneNames(bg)[gene_pos])),
    samples=sample_ids, 
    meas='FPKM', 
    colorby='transcript')

#Plot the mean FPKM values for each gene
plotMeans(ballgown::geneIDs(bg)[gene_pos], 
    bg,
    groupvar='population', 
    meas='FPKM', 
    colorby=c('transcript'),
    labelTranscripts=TRUE)


dev.off()



##################################################
############# Advanced Visualization #############
##################################################

# Open a PDF file where we will save some plots
outfile = paste0(path, "advanced_plots.pdf")
pdf(file=outfile)

#Read the results
results_transcripts = read.csv(paste0(path,
    "transcript_results.csv"))

results_genes = read.csv(paste0(path,
    "gene_results.csv"))

#Subset the results to include only significant results
# results_transcripts = subset(results_transcripts,
#     results_transcripts$pval<0.05)

# results_genes = subset(results_genes,
#     results_genes$pval<0.05)


gene_expression = as.data.frame(gexpr(bg_filt))

new_colnames <- c('HCM1', 'HCM2', 'HCM3', 'HCM4',
    'Healthy1', 'Healthy2', 'HCM5', 'HCM6',
    'Healthy3', 'Healthy4', 'Healthy5')

colnames(gene_expression) <- new_colnames

############## Transcripts per gene ##############

transcript_gene_table = indexes(bg_filt)$t2g
head(transcript_gene_table)

counts = table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)

hist(counts,
    breaks=50,
    col=2,
    xlab="Transcripts per gene",
    main="Distribution of transcript count per gene")

legend_text = c(paste("Genes with one transcript =", c_one),
    paste("Genes with more than one transcript =",
    c_more_than_one),
    paste("Max transcripts for single gene = ",
    c_max))

legend("topright", 
    legend_text,
    lty=NULL)


########## Transcript size distribution ##########
########## Something is wrong here ###############

# full_table <- texpr(bg_filt, 'all')
# lol <- full_table$length

# hist(lol,
#     breaks=50,
#     xlab="Transcript length (bp)",
#     main="Distribution of transcript lengths",
#     col="steelblue")
    

###########  distribution of DE values ###########

sig=which(results_genes$pval<0.05)

hist(results_genes[sig,"de"],
    breaks=50,
    col="seagreen",
    xlab="log2(Fold change) HCM vs Healthy",
    main="Distribution of differential expression values")

abline(v=-2,
    col="black",
    lwd=2,
    lty=2)

abline(v=2,
    col="black",
    lwd=2,
    lty=2)

legend("topleft",
    "Fold-change > 2",
    lwd=2,
    lty=2)

#################### Heatmap #####################

sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 1)
sig_tn_de = sigp[sigde,]
head(sig_tn_de)

mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

# Set the main title for the heatmap
main_title="Sig DE Genes in HCM vs Healthy by FPKM"
par(cex.main=0.8)
sig_genes_de=sig_tn_de[,"id"]
sig_gene_names_de=sig_tn_de[,"external_gene_name"]

# Get the gene expression data for the significant genes
data_columns=c(1:11)

# Transform the gene expression data to log2 scale
data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),
    data_columns])+1)

# Create the heatmap
heatmap.2(data,
    hclustfun=myclust,
    distfun=mydist,
    na.rm = TRUE,
    scale="none",
    dendrogram="both",
    margins=c(10,4),
    Rowv=TRUE,
    Colv=TRUE,
    symbreaks=FALSE,
    key=TRUE, 
    symkey=FALSE,
    density.info="none", 
    trace="none",
    main=main_title, 
    cexRow=0.5,
    cexCol=1, 
    labRow=sig_gene_names_de,
    col=rev(heat.colors(75)))


#################### Volcano Plot #####################

#Default all genes to "no change"
results_genes$diffexpressed <- "No"

#If log2Foldchange >2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$de > 1 & results_genes$pval < 0.05] <- "Up regulated genes in HCM"

#If log2Foldchange <-2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$de < -1 & results_genes$pval < 0.05] <- "Down regulated genes in HCM"

#Create a new column to store the gene names
results_genes$gene_label <- NA

#Write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label[results_genes$diffexpressed != "No"] <- results_genes$external_gene_name[results_genes$diffexpressed != "No"]

#Plot the volcano plot
ggplot(data=results_genes[results_genes$diffexpressed != "No",],
    aes(x=de, y=-log10(pval),
    label=gene_label,
    color = diffexpressed)) +
    xlab("log2Foldchange") +
    scale_color_manual(name = "Differentially expressed",
    values=c("blue", "red")) +
    geom_point() +
    theme_minimal() +
    geom_text_repel() +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    geom_point(data = results_genes[results_genes$diffexpressed == "No",],
    aes(x=de, y=-log10(pval)), colour = "black")


#################### PCA #####################

pca = (prcomp(t(fpkm)))

autoplot(pca,
    main="PCA of FPKM for the samples",
    data=pData(bg),
    colour = "population")

dev.off()

