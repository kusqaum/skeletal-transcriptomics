#gene expression preparation
library(tidyverse)
library(plyr)


# testfile <- readRDS("raw/rnaSeq/GSE99662/GSE99662_txi.RDS")
# testfile2 <- readRDS("raw/rnaSeq/GSE162691/GSE162691_txi.RDS")


# test_abundance <- testfile$abundance
# test_abundance2 <- testfile2$abundance
# head(testfile$abundance)
# testfile[["abundance"]]


#we have a vector of all the files
files <- list.files("raw/rnaSeq", full.names = F)
files <- files[-grep("test", files, fixed = T)]
files <- files[-grep("-log", files, fixed = T)]
files <- files[-grep(".tar.gz", files, fixed = T)]


#want to pass this vector to readRDS using lapply
abundance_data <- lapply(files, function(x) {
  f <- readRDS(paste0("raw/rnaSeq/",x,"/",x, "_txi.RDS"))
  abundance <- f$abundance
  return(abundance)
  
})


#human to human every transcript(geneID) is in right order

#create an empty list to separate each species with
allHuman <- list()
allMouse <- list()
# allRat <- list()
# allPig <- list()
# allCow <- list()
# allHorse <- list()
# allZebrafish <- list()
#loop through each df in list of abundance dfs
for (d in abundance_data){
  #search for those that have human IDs
  human <- grep(pattern = "ENSG", row.names(d))
  #search for those that have mouse IDs
  mouse <- grep(pattern = "ENSMUSG", row.names(d))
  #rat <- grep(pattern = "ENSRNOG", row.names(d))
  #pig <- grep(pattern = "ENSSSC", row.names(d))
  #cow <- grep(pattern = "ENSBTA", row.names(d))
  #horse <- grep(pattern = "ENSECAG", row.names(d))
  #zebrafish <- grep(pattern = "ENSDARG", row.names(d))
  if (any(human)){
    allHuman <- append(allHuman, list(d))
  }
  else if (any(mouse)){
    allMouse <- append(allMouse, list(d))
  }
  #else if (any(rat)){
   # allRat <- append(allRat, list(d))
  #}
  #else if(any(pig)){
   # allPig <- append(allPig, list(d))
  #}
  #else if (any(cow)){
   # allCow <- append(allCow, list(d))
  #}
  #else if (any(horse)){
   # allHorse <- append(allHorse, list(d))
  #}
  #else if (any(zebrafish)){
   # allZebrafish <- append(allZebrafish, list(d))
  #}
}



#bind all human DFs- change to a dataframe cause is matrix
humanDf <- data.frame(do.call("cbind", allHuman))

#now merge all mouse Dfs
#change all lists to dfs

allMouse_df <- lapply(allMouse, function(r){
  dataset <- data.frame(list(r))
})

# the same for the other animals data
#allRat_df <- lapply(allRat, function(s){
#  dataset <- data.frame(list(s))
#})

#allPig_df <- lapply(allPig, function(t){
 # dataset <- data.frame(list(t))
#})

#allCow_df <- lapply(allCow, function(u){
 # dataset <- data.frame(list(u))
#})

#allHorse_df <- lapply(allHorse, function(v){
 # dataset <- data.frame(list(v))
#})

#allZebrafish_df <- lapply(allZebrafish, function(w){
 # dataset <- data.frame(list(w))
#})


#mouse arent in same order so gonna add an extra column of their IDs
for (m in 1:length(allMouse_df)){
  allMouse_df[[m]]$ID <- rownames(allMouse_df[[m]])
}



# addd extra column for rats as well as the other species
#for (n in 1:length(allRat_df)) {
#  allRat_df[[n]]$ID <- rownames(allRat_df[[n]])
#}

#for (p in 1:length(allPig_df)) {
#  allPig_df[[p]]$ID <- rownames(allPig_df[[p]])
#}
  
#for (c in 1:length(allCow_df)) {
#  allCow_df[[c]]$ID <- rownames(allCow_df[[c]])
#}  

#for (h in 1:length(allHorse_df)) {
#  allHorse_df[[h]]$ID <- rownames(allHorse_df[[h]])
  
#}

#for (z in 1:length(allZebrafish_df)) {
#  allZebrafish_df[[z]]$ID <- rownames(allZebrafish_df[[z]])
  
#}

#now I need to map the animal ENSEMBL IDs to human
mouseDf <- join_all(allMouse_df, by = "ID", type = "full")
#ratDf <- join_all(allRat_df, by = "ID", type = "full")
#pigDf <- join_all(allPig_df, by = "ID", type = "full")
#cowDf <- join_all(allCow_df, by = "ID", type = "full")
#do horse and zebrafish also:

#going to leave this for now, just incase we want the column names instead of rownames later?
#rownames(mouseDf) <- mouseDf$ID; mouseDf$ID <- NULL

mouse_ensemble <- data.frame(mouseDf$ID)
#rat_ensembl <- data.frame(ratDf$ID)
write.table(mouse_ensemble, "raw/mouse_ensembleIDs.txt", sep = "\t", quote = F, row.names = F)

mouse_homology <- read.table("processed/Genes_Mouse_to_Human.txt", header = T)

mouseDf2 <- merge(mouseDf, mouse_homology, by = "ID")
#make the human IDs the row names and then drop those columns because don't need anymore
rownames(mouseDf2) <- mouseDf2$HumanEnsembl; mouseDf2$ID <- NULL; mouseDf2$HumanEnsembl <- NULL
####merge mouse df with human df####
human_mouseDf <- merge(mouseDf2, humanDf, by=0) 
#idk about this TBH I am not going to merge just yet.
#humanMouse <- merge(humanDf, mouseDf2, by = 0)
rownames(human_mouseDf) <- human_mouseDf$Row.names; human_mouseDf$Row.names <- NULL

#
#read in file created outside this script to create a df of just human protein coding genes
#from biomart
human_coding <- read.table("processed/human_coding_genes.txt", header = T, sep = "\t")
human_coding <- human_coding %>% select(ensembl_gene_id)
# humanDf$ensembl_gene_id <- rownames(humanDf) # uncomment this if you want to use just human gene expression data
human_mouseDf$ensembl_gene_id <- row.names(human_mouseDf)

# then merge with large humanDf so just retaining IDs common to both. remove the cols made for 
#merging
proteinCoding <- merge(human_coding, human_mouseDf, by ="ensembl_gene_id"); human_mouseDf$ensembl_gene_id <- NULL
#so now we only have 16984 genes
#humanPC <- merge(human_coding, humanDf, by = "ensembl_gene_id"); humanDf$ensembl_gene_id <- NULL#; humanPC$ensembl_gene_id <- NULL
#rownames(humanPC)<- humanPC$ensembl_gene_id # uncomment this line and linei above if you want to use human dfs only
rownames(proteinCoding)<- proteinCoding$ensembl_gene_id
forDimensionReduction <- proteinCoding
forDimensionReduction$ensembl_gene_id <- NULL
saveRDS(forDimensionReduction, "processed/geneExpressForDimRed.rds")
#write.table(humanPC, "processed/HumanCoding_rnaSeq.txt", row.names = F, sep = "\t", quote = F) # remember that one of the
#columns is called ensembl_gene_id for later

#so just keeping the protein coding genes from the human dataframe.leaves me with a total of 
#20,430 protein coding genes
#so just keeping the protein coding genes from the merged dataframe.leaves me with a total of 
#16984 protein coding genes

#Now read in IMPC labels
##### Integrating gene expression with mouse phenotype data
impcGenes <- read.table("processed/IMPC_phenotypeAssociations.txt", header = T)



#there are some duplicated genes that are associated and also not associated
# with a phenotype. I want to keep any gene that has been associated
#with a skeletal phenotype in any case:
impcGenes2 <- impcGenes[order(impcGenes[,"marker_symbol"], -impcGenes[,"significant"]),]
impcGenes2 <- impcGenes2[!duplicated(impcGenes2$marker_symbol),]
length(which(impcGenes2$significant=="TRUE"))#1376
length(which(impcGenes2$significant=="FALSE"))#7116

#read in mouse symbol to ensembl id data file
hgncAllianceHomology <- read.table("processed/hgncAllianceHomology.txt", header = T)

matchs_allHom <- match(impcGenes2$marker_symbol, hgncAllianceHomology[,1])
ensembl_id <- hgncAllianceHomology[matchs_allHom, 2]
impcGenes2$ID <- ensembl_id

sum(is.na(impcGenes2$ID))
# nA <- subset(impcGenes3, is.na(ID))
# the_null <- subset(impcGenes3, ID == "null")
# na <- impcGenes2 %>% filter(is.na(ID)) # there are 18 na values
# null <- impcGenes2 %>% filter(ID == "null") #and 5 null

#now write to table
# write.table(nA, "processed/cpg.txt", sep = "\t", row.names = F, quote = F)
#the rest NAs are cpgs
impcGenes2 <- impcGenes2 %>% filter(!is.na(ID)) %>%
  filter(ID != "null")


#so now I need to transfer all IMPC genes
impc_mouse_ensembl <- data.frame(impcGenes2$ID)
write.table(impc_mouse_ensembl, "processed/list_impc_ids.txt", row.names = F, sep = "\t", quote = F)
#so that I can convert them to human orthologs using biomart
#and then i need to read them in....



# impcHomology <- read.table("processed/allIMPC_homology.txt", header = T)
impcHomology <- read.table("processed/impcHomologyAllLatest.txt", header = T)
impcHomology_temporary <- impcHomology
#only managed to get mappings for 8003 genes out of 8474. 472 genes lost 

matchingIMPC <- match(impcGenes2$ID, impcHomology[,1])
ensemblHuman <- impcHomology[matchingIMPC, 2]
impcGenes2$ensembl_gene_id <- ensemblHuman
sum(is.na(impcGenes2$ensembl_gene_id))
#468 
#i have checked some of these and I don't think there is a mapping for them. So I will have to leave them out :
na_2 <- impcGenes2 %>% filter(is.na(ensembl_gene_id))


final_impcGenes <- impcGenes2 %>% select(ensembl_gene_id, significant) %>% filter(!is.na(ensembl_gene_id))
sum(final_impcGenes$significant == "TRUE")#1320
sum(final_impcGenes$significant == "FALSE")#6681


#then can merge this df with humanPC by "ensemblgene_id"
# dataWithLabels <- merge(humanPC, final_impcGenes, by = "ensembl_gene_id") # uncomment if you are just working with human df
#THE GE df has still got the rownames as a column called ensembl_gene_id. so will merge using that
mouseHumanWithLabels <- merge(proteinCoding, final_impcGenes, by="ensembl_gene_id")
#now change ensembl IDs to rownames:
# rownames(dataWithLabels)<- dataWithLabels$ensembl_gene_id; dataWithLabels$ensembl_gene_id <- NULL
rownames(mouseHumanWithLabels)<- mouseHumanWithLabels$ensembl_gene_id; mouseHumanWithLabels$ensembl_gene_id <-NULL
# write.table(dataWithLabels, "processed/geneExpressionDataWithLabels.txt", row.names = T, sep = "\t", quote = F)
saveRDS(mouseHumanWithLabels, "processed/mouseHumanWithLabels.rds")

##########import new impc phenotypes######
#just looking at a mortality
impcGenesMortality <- read.table("processed/mortalityPhenotype.txt", header = T, sep = "\t")
str(impcGenesMortality$significant)
impcGenesMortality$significant <- as.logical(impcGenesMortality$significant)
str(impcGenesMortality$significant)
head(impcGenesMortality)
str(impcGenesMortality)

impcGenesNew <- impcGenesMortality[order(impcGenesMortality[,"marker_symbol"], -impcGenesMortality[,"significant"]),]
head(impcGenesNew)
impcGenesNew2 <- impcGenesNew[!duplicated(impcGenesNew$marker_symbol),]
length(which(impcGenesNew2 == "TRUE"))
length(which(impcGenesNew2=="FALSE"))
#there is an imbalance so will again have to downsample

# now match the impc gene symbols to gene symbols from the mgi database
head(hgncAllianceHomology)
mortalityMatch_symb <- match(impcGenesNew2$marker_symbol, hgncAllianceHomology[,1])
mouse_ensemb_id_matching <- hgncAllianceHomology[mortalityMatch_symb, 2]
#add on the mouse ensembl ids as extra column now
impcGenesNew2$ensembl_id <- mouse_ensemb_id_matching
####'
impcGenesNewDF <- data.frame(mouseEnsembl = impcGenesNew2$ensembl_id)
write.table(impcGenesNewDF, "processed/mortalityPhen.txt", row.names = F, sep = "\t", quote = F)
####
#now map mouse ensembl to human ensembl:
mortalityPhenotyepHomology <- read.table("processed/mouseMortalityGenesHomology.txt", header = T, sep = "\t")
match_impcMortality <- match(impcGenesNew2$ensembl_id, mortalityPhenotyepHomology$Mouse_ensembl)
ensemblMortality <- mortalityPhenotyepHomology[match_impcMortality, 2]
impcGenesNew2$ensembl_gene_id <- ensemblMortality
head(impcGenesNew2)

processedMortalityGenes <- impcGenesNew2 %>% select(ensembl_gene_id, significant)%>%
  filter(!is.na(ensembl_gene_id))
head(processedMortalityGenes)

rownames(processedMortalityGenes)<- processedMortalityGenes$ensembl_gene_id; processedMortalityGenes$ensembl_gene_id<-NULL
head(processedMortalityGenes)
processedMortalityGenes$significant <- as.factor(processedMortalityGenes$significant)
write.table(processedMortalityGenes ,"processed/processedMortalityLabels.txt", sep = "\t", quote = F)



