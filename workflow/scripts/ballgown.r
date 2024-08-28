#!/usr/bin/env Rscript

library(ballgown)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

# Load the data
bg <- ballgown(dataDir=args[1], samplePattern=args[2], meas='all')

# Filter the data
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Run the differential expression analysis

bg_filt <- stattest(bg_filt, feature='gene', covariate=args[3])

# Write the results to a file

write.table(bg_filt, file=args[4], sep="\t", quote=FALSE, row.names=FALSE)

# Plot the results

plotData <- merge(as.data.frame(transcripts(bg_filt)), as.data.frame(genes(bg_filt)), by="gene_id")

