#this script extracts all the protein coding human genes from biomart
library(tidyverse)
library(biomaRt)

e <- useMart("ensembl")
ensemblHuman <- useDataset("hsapiens_gene_ensembl", mart = e)
pc_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "biotype",
                  values = c("protein_coding"),
                  mart = ensemblHuman)
length(unique(pc_genes$ensembl_gene_id))

pc_genesF <- pc_genes%>%select(ensembl_gene_id)
# colnames(pc_genes)[1] <- "ID"

humn_ens <- read.table("human_ensemblIDs.txt", header = T)
# colnames(humn_ens)[1] <- "ID"
write.table(pc_genes, "processed/human_coding_genes.txt", quote = F, row.names = F, sep = "\t")
#df <- merge(humn_ens, pc_genes, by = "ID")
