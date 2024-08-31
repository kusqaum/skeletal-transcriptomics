#mapping all of the IMPC ensembl IDs to human ensembl IDs
library(tidyverse)
library(biomaRt)

impcList <- read.table("list_imc_ids.txt", header = T)
impcList <- c(impcList$impcGenes2.ID)

h <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
m <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)

#also get the gene symbols as well
mapped_list <- getLDS(attributes = c("ensembl_gene_id", "mgi_symbol"),
                      filters = "ensembl_gene_id",
                      values = impcList,
                      mart = m,
                      attributesL = c("ensembl_gene_id", "hgnc_symbol"),
                      martL = h)


#remove those that do not have corresponding gene symbols as there are some duplicates that 
#have the corresponding gene symbol
filteredmapped <- mapped_list %>% filter(HGNC.symbol!="") %>%
  distinct(Gene.stable.ID, .keep_all = T) %>% distinct(Gene.stable.ID.1, .keep_all = T)


allIMPC_homology <- filteredmapped %>% select(1,3)
colnames(allIMPC_homology)[1] <- "Mouse_ensembl"
colnames(allIMPC_homology)[2] <- "gene_ensembl_id"
write.table(allIMPC_homology, "allIMPC_homology.txt", sep = "\t",
            row.names = F, quote = F)