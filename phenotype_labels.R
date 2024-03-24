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
