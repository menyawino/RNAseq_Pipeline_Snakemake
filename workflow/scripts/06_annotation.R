library(AnnotationHub)
library(ensembldb)
library(purrr)

ah <- AnnotationHub()

human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens <- human_ens[["AH75011"]]

genes(human_ens, return.type = "data.frame") %>% View()

head(results_transcripts)

#Create annotation dataframe
annotations_ahb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, entrezid, gene_biotype) %>% 
  dplyr::filter(gene_id %in% results_transcripts$ensembl_gene_id)

# Keep the first identifier for these multiple mapping cases
annotations_ahb$entrezid <- map(annotations_ahb$entrezid
   ,1) %>%  unlist()

# Identify and keep non-duplicate genes
non_duplicates_idx <- which(duplicated(annotations_ahb$entrezid) == FALSE)
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

write.csv(annotations_ahb, 
    paste0(path, "annotations_ahb.csv"),
    row.names=FALSE)
