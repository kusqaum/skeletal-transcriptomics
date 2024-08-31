# this script I am getting the genes annotated to a skeletal phenotype in the MGI database
# and map the human gene symbols to ensembl identifiers

library(tidyverse)
library(biomaRt)

externalDatabase <- read.delim("https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt", header = F,
                               sep = "\t")
externalDatabase$V6 <- NULL
colnames(externalDatabase)<- c("humanMarkerSymbol", "humanEntrezGeneID", "mouseMarkerSymbol", "mgiMarkerID", "mammalianPhenotypeID")
head(externalDatabase[1:5,1:5])
dim(externalDatabase)



exDb <- externalDatabase %>% filter(grepl("MP:0005390", mammalianPhenotypeID))
head(exDb[1:4,1:4])
dim(exDb)

# i need to convert the human marker symbols to ensembl gene IDs

e <- useMart("ensembl")

ensemblHuman <- useDataset("hsapiens_gene_ensembl", mart = e)
genesEnsemblandSymbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                               filters = "biotype",
                               values = c("protein_coding"),
                               mart = ensemblHuman)

length(which(genesEnsemblandSymbol$hgnc_symbol==""))
genesEnsemblandSymbol <- genesEnsemblandSymbol %>% filter(hgnc_symbol != "")
length(which(genesEnsemblandSymbol ==""))
head(genesEnsemblandSymbol)
head(exDb[1:4,1:4])
matchGS <- match(exDb[,1], genesEnsemblandSymbol[,2])
matchingGeneSymb<- genesEnsemblandSymbol[matchGS, 1]
exDb$ensembl_gene_id <- matchingGeneSymb
head(exDb)
mgiGenes <- data.frame(ensembl_gene_id = exDb[,6])
length(which(is.na(mgiGenes)))
mgiGenes$significant <- "TRUE"
mgiGenes$significant <- as.factor(mgiGenes$significant)
# 49 genes that have failed to map. 
mgiGenes <- mgiGenes %>% filter(!is.na(ensembl_gene_id))
head(mgiGenes)
# looks good

write.table(mgiGenes, "processed/annotatedMGIgenes.txt", sep = "\t", quote = F, row.names = F)
write.table(exDb, "processed/mgiHumanPhenotype.txt", sep = "\t", quote = F, row.names = F)
