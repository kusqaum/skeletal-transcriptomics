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
allRat <- list()
allPig <- list()
allCow <- list()
allHorse <- list()
allZebrafish <- list()
#loop through each df in list of abundance dfs
for (d in abundance_data){
  #search for those that have human IDs
  human <- grep(pattern = "ENSG", row.names(d))
  #search for those that have mouse IDs
  mouse <- grep(pattern = "ENSMUSG", row.names(d))
  rat <- grep(pattern = "ENSRNOG", row.names(d))
  pig <- grep(pattern = "ENSSSC", row.names(d))
  cow <- grep(pattern = "ENSBTA", row.names(d))
  horse <- grep(pattern = "ENSECAG", row.names(d))
  zebrafish <- grep(pattern = "ENSDARG", row.names(d))
  if (any(human)){
    allHuman <- append(allHuman, list(d))
  }
  else if (any(mouse)){
    allMouse <- append(allMouse, list(d))
  }
  else if (any(rat)){
    allRat <- append(allRat, list(d))
  }
  else if(any(pig)){
    allPig <- append(allPig, list(d))
  }
  else if (any(cow)){
    allCow <- append(allCow, list(d))
  }
  else if (any(horse)){
    allHorse <- append(allHorse, list(d))
  }
  else if (any(zebrafish)){
    allZebrafish <- append(allZebrafish, list(d))
  }
}



#bind all human DFs- change to a dataframe cause is matrix
humanDf <- data.frame(do.call("cbind", allHuman))

#now merge all mouse Dfs
#change all lists to dfs

allMouse_df <- lapply(allMouse, function(r){
  dataset <- data.frame(list(r))
})

# the same for the other animals data
allRat_df <- lapply(allRat, function(s){
  dataset <- data.frame(list(s))
})

allPig_df <- lapply(allPig, function(t){
  dataset <- data.frame(list(t))
})

allCow_df <- lapply(allCow, function(u){
  dataset <- data.frame(list(u))
})

allHorse_df <- lapply(allHorse, function(v){
  dataset <- data.frame(list(v))
})

allZebrafish_df <- lapply(allZebrafish, function(w){
  dataset <- data.frame(list(w))
})



for (m in 1:length(allMouse_df)){
  allMouse_df[[m]]$ID <- rownames(allMouse_df[[m]])
}



# addd extra column for rats as well as the other species
for (n in 1:length(allRat_df)) {
  allRat_df[[n]]$ID <- rownames(allRat_df[[n]])
}

for (p in 1:length(allPig_df)) {
  allPig_df[[p]]$ID <- rownames(allPig_df[[p]])
}
  
for (c in 1:length(allCow_df)) {
  allCow_df[[c]]$ID <- rownames(allCow_df[[c]])
}  

for (h in 1:length(allHorse_df)) {
  allHorse_df[[h]]$ID <- rownames(allHorse_df[[h]])
  
}

for (z in 1:length(allZebrafish_df)) {
  allZebrafish_df[[z]]$ID <- rownames(allZebrafish_df[[z]])
  
}

#now I need to map the animal ENSEMBL IDs to human
mouseDf <- join_all(allMouse_df, by = "ID", type = "full")
#ratDf <- join_all(allRat_df, by = "ID", type = "full")
#pigDf <- join_all(allPig_df, by = "ID", type = "full")
#cowDf <- join_all(allCow_df, by = "ID", type = "full")
#do horse and zebrafish also:

#going to leave this for now, just incase we want the column names instead of rownames later?
rownames(mouseDf) <- mouseDf$ID; mouseDf$ID <- NULL

mouse_ensemble <- data.frame(mouseDf$ID)
rat_ensembl <- data.frame(ratDf$ID)
write.table(mouse_ensemble, "raw/mouse_ensembleIDs.txt", sep = "\t", quote = F, row.names = F)
write.table(rat_ensembl, "raw/rat_ensembleIDs.txt", sep = "\t", quote = F, row.names = F)

mouse_homology <- read.table("processed/Genes_Mouse_to_human.txt", header = T)

mouseDf2 <- merge(mouseDf, mouse_homology, by = "ID")
#make the human IDs the row names and then drop those columns because don't need anymore
rownames(mouseDf2) <- mouseDf2$HumanEnsembl; mouseDf2$ID <- NULL; mouseDf2$HumanEnsembl <- NULL

#idk about this TBH I am not going to merge just yet.
#humanMouse <- merge(humanDf, mouseDf2, by = 0)
#rownames(humanMouse) <- humanMouse$Row.names; humanMouse$Row.names <- NULL

#going to write out the human genes so that i can use biomart to just get the protein coding 
#human genes - might not actually need to write out the file tbh
#human_ensembl <- data.frame(rownames(humanDf))
#write.table(human_ensembl, "raw/human_ensemblIDs.txt", sep = "\t", quote = F, row.names = F)

#read in file created outside this script to create a df of just human protein coding genes
#from biomart
human_coding <- read.table("processed/human_coding_genes.txt", header = T)
humanDf$ensembl_gene_id <- rownames(humanDf)
# then merge with large humanDf so just retaining IDs common to both. remove the cols made for 
#merging
humanPC <- merge(human_coding, humanDf, by = "ensembl_gene_id"); humanDf$ensembl_gene_id <- NULL; humanPC$ensembl_gene_id <- NULL

#so just keeping the protein coding genes from the human dataframe.leaves me with a total of 
#20,430 protein coding genes
humanPC

