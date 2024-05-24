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

#create an empty list
allHuman <- list()
allMouse <- list()
allRat <- list()
allPig <- list()
allCow <- list()
#loop through each df in list of abundance dfs
for (d in abundance_data){
  #search for those that have human IDs
  human <- grep(pattern = "ENSG", row.names(d))
  #search for those that have mouse IDs
  mouse <- grep(pattern = "ENSMUSG", row.names(d))
  rat <- grep(pattern = "ENSRNOG", row.names(d))
  pig <- grep(pattern = "ENSSSC", row.names(d))
  cow <- grep(pattern = "ENSBTA", row.names(d))
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
}



#bind all human DFs
humanDf <- data.frame(do.call("cbind", allHuman))

#now merge all mouse Dfs
#change all lists to dfs
allMouse_df <- lapply(allMouse, function(r){
  dataset <- data.frame(list(r))
  
})
for (m in 1:length(allMouse_df)){
  allMouse_df[[m]]$ID <- rownames(allMouse_df[[m]])
}

mouseDf <- join_all(allMouse_df, by = "ID", type = "full")
rownames(mouseDf) <- mouseDf$ID; mouseDf$ID <- NULL


#need some logical step for mouse IDs- what is it? how to convert to human and then merge (need to go from mouse to human)
#for (i in files) {
 # f <- readRDS(paste0("raw/rnaSeq/",i,"/",i,"_txi.RDS"))
#}
