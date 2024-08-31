# this script is mapping mouse ensembl IDs in gene expression data
# to human ensembl IDs using biomaRt
library(biomaRt)
library(tidyverse)
#read in mouse ensembl ids
mouse_ensemblIDs <- read.table("mouse_ensemblIDs.txt", header = T)
mouseIDSl <- c(mouse_ensemblIDs$mouseDf.ID)

#h <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", verbose = T, host = "https://dec2021.archive.ensembl.org")
h <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
m <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)

#we just want protein coding genes
genesV222 = getLDS(attributes =  c("ensembl_gene_id", "gene_biotype"),
                   filters = c("ensembl_gene_id","biotype"),
                   values = list(mouseIDSl, "protein_coding"),
                   mart = m, 
                   attributesL = c("ensembl_gene_id"), 
                   martL = h, 
                   uniqueRows=T)


genesV4 = getLDS(attributes =  c("ensembl_gene_id", "mgi_symbol"),
                 filters = c("ensembl_gene_id"),
                 values = list(mouseIDSl),
                 mart = m, 
                 attributesL = c("ensembl_gene_id", "hgnc_symbol"), 
                 martL = h, 
                 uniqueRows=T)

humanv4 <- genesV4 %>% distinct(Gene.stable.ID.1, .keep_all = T)
write.table(humanv4, "mouseHumanEnsembl.txt", sep = "\t", row.names = F, quote = F)
humanv41 <- humanv4 %>% filter(HGNC.symbol=="")
write.table(humanv41, "noHumanSymb.txt", sep = "\t", row.names = F, quote = F)


humanx2 <- genesV222 %>% distinct(Gene.stable.ID.1, .keep_all = T) #17978
colnames(humanx2)[1] <- "ID"
colnames(humanx2)[3] <- "HumanEnsembl"
humanx2 <- humanx2 %>% select(ID, HumanEnsembl)
write.table(humanx2, "processed/Genes_Mouse_to_human.txt", quote = F, row.names = F)