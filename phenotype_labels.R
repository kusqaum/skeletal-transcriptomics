library(tidyverse)

data <- read.csv("statistical-results-ALL.csv")
filt <- data[-which(data$top_level_mp_term_id == ""), ]

#filtering for MP:0005390- skeleton phenotype
top_level <- filt %>% filter(grepl("MP:0005390", top_level_mp_term_id)) 

#filter for significance
true <- top_level %>% filter(significant == "true")
length(unique(true$marker_symbol))
false <- top_level %>% filter(significant == "false")
length(unique(false$marker_symbol))
#now we have 1376 genes that are associated to skeletal phenotype
#and 8479 that are not

#select columns and keep unique genes
small_true <- true %>% select(marker_accession_id, marker_symbol, significant)
small_false <- false %>% select(marker_accession_id, marker_symbol, significant)
unique_true <- small_true %>% distinct(marker_symbol, .keep_all = T) # T because wanna keep all columns 
unique_false <- small_false %>% distinct(marker_symbol, .keep_all = T)

#combine
allgenes <- bind_rows(unique_true, unique_false)
allgenes <- as.logical(allgenes$significant)
write.table(allgenes, "processed/IMPC_phenotypeAssociations.txt", sep = "\t", quote = F, row.names = F)



####
# extract all IMPC genes associated with mortality/aging:
mortality_aging <- filt %>% filter(grepl("MP:0010768", top_level_mp_term_id))
mortalityT <- mortality_aging %>% filter(significant == "true") %>% 
  distinct(marker_symbol, .keep_all = T) %>% select(marker_symbol, significant, marker_accession_id)
mortalityF <- mortality_aging %>% filter(significant=="false") %>% distinct(marker_symbol, .keep_all = T) %>%
  select(marker_symbol, significant , marker_accession_id)
allMortality <- rbind(mortalityT, mortalityF)

allMortality$significant <- as.logical(allMortality$significant)
allMortalityPhen <- allMortality[order(allMortality[,"marker_symbol"], -allMortality[,'significant']),]
allMortAging <- allMortalityPhen[!duplicated(allMortalityPhen$marker_symbol),]
length(which(allMortAging$significant=="TRUE"))
length(which(allMortAging$significant=="FALSE"))
head(allMortAging)
dim(allMortAging)
length(which(allMortAging$significant=="TRUE"))
length(which(allMortAging$significant=="FALSE"))
write.table(allMortAging, "processed/mortalityPhenotype.txt", sep = "\t", quote = F, row.names = F)

