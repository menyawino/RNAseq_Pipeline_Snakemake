library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

path = "/mnt/e/MYF/Ballgown/"

results_transcripts = read.csv("/mnt/e/MYF/Ballgown/transcript_results.csv")
results_genes = read.csv("/mnt/e/MYF/Ballgown/gene_results.csv")
annotations_ahb = read.csv("/mnt/e/MYF/Ballgown/annotations_ahb.csv")

res_ids <- inner_join(results_genes,
    annotations_ahb,
    by=c("ensembl_gene_id"="gene_id"))    

# filter over expressed genes
res_ids_OE <- dplyr::filter(res_ids, fc > 1)

# filter under expressed genes
res_ids_UE <- dplyr::filter(res_ids, fc < 1)


################################################
######### Over-representation analysis #########
################################################

allOE_genes <- as.character(res_ids_OE$ensembl_gene_id)

sigOE <- dplyr::filter(res_ids_OE, pval < 0.05)

sigOE_genes <- as.character(sigOE$ensembl_gene_id)

sigOE_genes <- unique(sigOE_genes)

# Perform GO analysis
ego <- enrichGO(gene = sigOE_genes, 
    keyType = "ENSEMBL",
    OrgDb = org.Hs.eg.db, 
    ont = "BP", 
    pAdjustMethod = "BH", 
    qvalueCutoff = 0.05, 
    readable = TRUE)

save(ego, file = paste0(path, "ego.RData"))

load(paste0(path, "ego.RData"))

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

# write.csv(cluster_summary,
#     "/mnt/e/MYF/Ballgown/cluster_summary_genes.csv")

# PDF file open
outfile = paste0(path, "func_plotss.pdf")
pdf(file=outfile)

dotplot(ego, showCategory=15)

compare_cluster_GO_emap <- enrichplot::pairwise_termsim(ego, semData = d)

emapplot(compare_cluster_GO_emap, 
    layout = "nicely",
    showCategory = 27)


# Category netplot shows the relationships between 
# the genes associated with the top five most 
# significant GO terms and the fold changes of the 
# significant genes associated with these terms (color)

OE_foldchanges <- sigOE$fc

names(OE_foldchanges) <- sigOE$ensembl_gene_id

cnetplot(ego, 
    categorySize="pvalue", 
    showCategory = 5, 
    foldChange=OE_foldchanges, 
    vertex.label.font=3)

dev.off()

################################################
########### Functional class scoring ###########
################################################

########## GSEA ##########

# Remove any NA values and duplicates
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]

# Extract the foldchanges
foldchanges <- res_entrez$fc

# Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrezid

# Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)


# GSEA using gene sets from KEGG pathways
# ordered named vector of fold changes Entrez IDs are names
gseaKEGG <- gseKEGG(geneList = foldchanges,
    organism = "hsa",
    minGSSize = 2,
    pvalueCutoff = 0.05,
    verbose = FALSE)

# GSEA using gene sets from GO terms
gseaGO <- gseGO(geneList = foldchanges, 
    OrgDb = org.Hs.eg.db, 
    ont = 'BP',
    minGSSize = 20,
    pvalueCutoff = 0.05,
    verbose = FALSE) 


save(gseaKEGG, gseaGO, file = paste0(path, "gsea.RData"))

load(paste0(path, "gsea.RData"))

# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
gseaGO_results <- gseaGO@result


# Write GSEA results to file
View(gseaKEGG_results)
View(gseaGO_results)

write.csv(gseaKEGG_results,
    "/mnt/e/MYF/Ballgown/gseaKEGG_results.csv")
write.csv(gseaGO_results,
    "/mnt/e/MYF/Ballgown/gseaGO_results.csv")



# Load the GSEA results
gseaKEGG_results <- read.csv("/mnt/e/MYF/Ballgown/gseaKEGG_results.csv")
gseaGO_results <- read.csv("/mnt/e/MYF/Ballgown/gseaGO_results.csv")

# Plot the GSEA plot for a single enriched pathway
gseaplot(gseaKEGG, geneSetID = 'hsa01230')

# Plot the GSEA plot for all enriched pathways
# first unload dplyr to avoid conflicts
detach("package:dplyr", unload=TRUE) 

# Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
    pathway.id = "hsa05150",
    species = "hsa",
    limit = list(gene = 10, # value gives the max/min limit for foldchanges
    cpd = 1))

# Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
    pathview(gene.data = foldchanges,
        pathway.id = gseaKEGG_results$ID[x],
        species = "hsa", 
        limit = list(gene = 10,
        cpd = 1))
}

purrr::map(9:length(gseaKEGG_results$ID), get_kegg_plots)

dev.off()