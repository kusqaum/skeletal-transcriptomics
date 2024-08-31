this script I am getting the genes annotated to a skeletal phenotype in the HPO
# and map the human gene symbols to ensembl identifiers

library(tidyverse)
library(biomaRt)


humanPhenotype <- read.delim("genes_for_HP_0000924", col.names = c("geneID", "geneSymbol"), sep = "\t", strip.white = T)
# the reason match() wasn't working was because when read the hp files in , there were leading spaces at the start of each char
# so removed them using strip.white !!!!!!!!!!

e <- useMart("ensembl")

ensemblHuman <- useDataset("hsapiens_gene_ensembl", mart = e)
genesHumanEnsemblandSymbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                    filters = "biotype",
                                    values = c("protein_coding"),
                                    mart = ensemblHuman)
length(which(genesHumanEnsemblandSymbol$hgnc_symbol==""))
genesHumanEnsemblandSymbolF <- genesHumanEnsemblandSymbol %>% filter(hgnc_symbol != "")
length(which(genesHumanEnsemblandSymbolF ==""))
head(genesHumanEnsemblandSymbolF)
head(humanPhenotype)
hSP <- match(humanPhenotype[,2], genesHumanEnsemblandSymbolF[,2])
matchingHSP <- genesHumanEnsemblandSymbolF[hSP, 1]
humanPhenotype$ensembl_gene_id <- matchingHSP
head(humanPhenotype)

hpoGenes <- data.frame(ensembl_gene_id = humanPhenotype[,3])
hpoGenes$significant <- "TRUE"
str(hpoGenes$significant)
hpoGenes$significant <- as.factor(hpoGenes$significant)
length(which(is.na(hpoGenes$ensembl_gene_id)))

hpoGenes <- hpoGenes %>% filter(!is.na(hpoGenes$ensembl_gene_id))
3294-3258                                
write.table(hpoGenes, "processed/humphenotOntGenes.txt", sep = "\t", row.names = F, quote = F)
write.table(humanPhenotype, "processed/hpoHumanPhenotype.txt", sep = "\t", row.names = F, quote = F)
