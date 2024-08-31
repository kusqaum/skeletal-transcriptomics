# this script is extracting mouse gene symbols to mouse ensembl id conversion table
library(tidyverse)
hgncAllianceHomology <- read.delim("https://www.informatics.jax.org/downloads/reports/HGNC_AllianceHomology.rpt")
hgncAllianceHomology2 <- hgncAllianceHomology %>% select(1,9)
write.table(hgncAllianceHomology2, "processed/hgncAllianceHomology.txt", quote = F, sep = "\t", 
            row.names = F)