#gene expression preparation
library(tidyverse)
library(plyr)


testfile <- readRDS("raw/rnaSeq/GSE99662/GSE99662_txi.RDS")
testfile2 <- readRDS("raw/rnaSeq/GSE162691/GSE162691_txi.RDS")


test_abundance <- testfile$abundance
test_abundance2 <- testfile2$abundance
head(testfile$abundance)
testfile[["abundance"]]


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

mouse_homology <- read.table("processed/Genes_Mouse_to_human.txt", header = T)

mouseDf2 <- merge(mouseDf, mouse_homology, by = "ID")
#make the human IDs the row names and then drop those columns because don't need anymore
rownames(mouseDf2) <- mouseDf2$HumanEnsembl; mouseDf2$ID <- NULL; mouseDf2$HumanEnsembl <- NULL
####merge mouse df with human df####
human_mouseDf <- merge(mouseDf2, humanDf, by=0)
#idk about this TBH I am not going to merge just yet.
#humanMouse <- merge(humanDf, mouseDf2, by = 0)
rownames(human_mouseDf) <- human_mouseDf$Row.names; human_mouseDf$Row.names <- NULL

#going to write out the human genes so that i can use biomart to just get the protein coding 
#human genes - might not actually need to write out the file tbh I'll just read the pc genes in straight from 
#biomart
#human_ensembl <- data.frame(rownames(humanDf))
#write.table(human_ensembl, "raw/human_ensemblIDs.txt", sep = "\t", quote = F, row.names = F)

#read in file created outside this script to create a df of just human protein coding genes
#from biomart
human_coding <- read.table("processed/human_coding_genes.txt", header = T)
humanDf$ensembl_gene_id <- rownames(humanDf)
human_mouseDf$ensembl_gene_id <- row.names(human_mouseDf)

# then merge with large humanDf so just retaining IDs common to both. remove the cols made for 
#merging
proteinCoding <- merge(human_coding, human_mouseDf, by ="ensembl_gene_id"); human_mouseDf$ensembl_gene_id <- NULL
humanPC <- merge(human_coding, humanDf, by = "ensembl_gene_id"); humanDf$ensembl_gene_id <- NULL#; humanPC$ensembl_gene_id <- NULL
rownames(humanPC)<- humanPC$ensembl_gene_id
rownames(proteinCoding)<- proteinCoding$ensembl_gene_id
write.table(humanPC, "processed/HumanCoding_rnaSeq.txt", row.names = F, sep = "\t", quote = F) # remember that one of the
#columns is called ensembl_gene_id for later

#so just keeping the protein coding genes from the human dataframe.leaves me with a total of 
#20,430 protein coding genes

##### Integrating gene expression with mouse phenotype data
impcGenes <- read.table("processed/IMPC_phenotypeAssociations.txt", header = T)

#impcGenes_ordered <- impcGenes[order(impcGenes$marker_symbol),]


#there are some duplicated genes that are associated and also not associated
# with a phenotype. I want to keep any gene that has been associated
#with a skeletal phenotype in any case:
impcGenes2 <- impcGenes[order(impcGenes[,"marker_symbol"], -impcGenes[,"significant"]),]
impcGenes2 <- impcGenes2[!duplicated(impcGenes2$marker_symbol),]

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

#fill in the ones that we can...
impcGenes2$ID[is.na(impcGenes2$ID) & impcGenes2$marker_symbol == "4932438H23Rik"] <- "ENSMUSG00000039851"
impcGenes2$ID[is.na(impcGenes2$ID) & impcGenes2$marker_symbol == "Ankrd36"] <- "ENSMUSG00000020481"
impcGenes2$ID[is.na(impcGenes2$ID) & impcGenes2$marker_symbol == "B430306N03Rik"] <- "ENSMUSG00000043740"
impcGenes2$ID[is.na(impcGenes2$ID) & impcGenes2$marker_symbol == "Gm11639"] <- "ENSMUSG00000040838"
impcGenes2$ID[is.na(impcGenes2$ID) & impcGenes2$marker_symbol == "Gm2694"] <- "ENSMUSG00000097248"

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


#need to now map the mouse ensembl IDs to human IDs...
#in_homology <- match(impcGenes2$ID, mouse_homology[,1])
# so now I need to gather the non mapped impc genes:

impcHomology <- read.table("processed/allIMPC_homology.txt", header = T)
#only managed to get mappings for 8015 genes out of 8474. 459 genes lost :/
#mouse_homology2 <- mouse_homology
#colnames(mouse_homology2)[1] <- "Mouse_ensembl"
#colnames(mouse_homology2)[2] <- "Human_ensembl"
#mouse_homology_forIMPC <- rbind(mouse_homology2, impcHomology)

matchingIMPC <- match(impcGenes2$ID, impcHomology[,1])
ensemblHuman <- impcHomology[matchingIMPC, 2]
impcGenes2$ensembl_gene_id <- ensemblHuman
sum(is.na(impcGenes2$ensembl_gene_id))
#i have checked some of these and I don't think there is a mapping for them. So I will have to leave them out :
na_2 <- impcGenes2 %>% filter(is.na(ensembl_gene_id))


final_impcGenes <- impcGenes2 %>% select(ensembl_gene_id, significant) %>% filter(!is.na(ensembl_gene_id))
sum(final_impcGenes$significant == "TRUE")#1323
sum(final_impcGenes$significant == "FALSE")#6692


#then can merge this df with humanPC by "ensemblgene_id"
dataWithLabels <- merge(humanPC, final_impcGenes, by = "ensembl_gene_id")
mouseHumanWithLabels <- merge(proteinCoding, final_impcGenes, by="ensembl_gene_id")
#now change ensembl IDs to rownames:
rownames(dataWithLabels)<- dataWithLabels$ensembl_gene_id; dataWithLabels$ensembl_gene_id <- NULL
rownames(mouseHumanWithLabels)<- mouseHumanWithLabels$ensembl_gene_id; mouseHumanWithLabels$ensembl_gene_id <-NULL
write.table(dataWithLabels, "processed/geneExpressionDataWithLabels.txt", row.names = T, sep = "\t", quote = F)
saveRDS(mouseHumanWithLabels, "processed/mouseHumanWithLabels.rds")
#mapped <- mouse_homology[in_homology, 2]
#impcGenes2$ID_humn <- mapped
#not_mapped <- impcGenes2 %>% filter(is.na(ID_humn))

#not_mapped_ids <- data.frame(not_mapped$ID)
#write.table(not_mapped_ids, "processed/not_mapped_IMPCgenes.txt", quote = F, row.names = F, sep = "\t")

#impcGenes_mapped <- merge(impcGenes2, mouse_homology, by = "ID")
#colnames(impcGenes_mapped)[5] <- "ensembl_gene_id"



#onlyHuman <- read.table("processed/HumanCoding_rnaSeq.txt", header = T)
#now i want to merge with mouse data...


