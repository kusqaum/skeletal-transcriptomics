#integrate network data with gene expression data
library(tidyverse)
#might move this to a different script
library(tidymodels)
library(themis)
library(vip)
library(doParallel)
library(foreach)
library(NMF)


labelsForNet <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
labelsForNet$significant <- as.factor(labelsForNet$significant)


#embeddings created using pecanpy in python

#importing node2vec embedding created in python
#split into 500 dims because that is the dimensions used for pecanPy
readNetworkData <- function(embData, dimens){
  #dimens is the number of dimensions used for pecanpy
  embedding <- read.table(embData, sep = "\t", skip = 1)
  emb <- embedding %>%
    separate(colnames(embedding), c("hgnc_symbol", paste0("D", 1:dimens)), sep = " ")
  #rownames(emb) <- emb$D0; emb$D0 <- NULL
  return(emb)
}


emb <- readNetworkData("processed/networkEdgeList.emb", dimens = 500)
dim(emb)

#now need to convert the genes to ensembl IDs
netGenes <- data.frame(hgnc_symbol = emb$hgnc_symbol)
dim(netGenes)
# write out now to convert
write.table(netGenes, "processed/networkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)

##now read in the mapped file

networkGenesWithEnsembl <- read.table("raw/allNetworkGenesMapped.txt", sep = "\t", header = T)

networkData <- merge(networkGenesWithEnsembl, emb, by="hgnc_symbol")
head(networkData[1:4,1:3])
rownames(networkData)<- networkData$ensembl_gene_id; networkData$hgnc_symbol<- NULL; networkData$ensembl_gene_id<-NULL
str(networkData$D2)
dim(networkData)
head(networkData[,488:500])

#now need to merge with gene exp
networkData <- mutate_all(networkData, function(x) as.numeric(as.character(x)))
str(networkData$D1)
networkDataLabelled <- merge(networkData, labelsForNet, by=0); rownames(networkDataLabelled) <- networkDataLabelled$Row.names; networkDataLabelled$Row.names <- NULL
saveRDS(networkData, "processed/processedNetworkEmb.rds")
saveRDS(networkDataLabelled, "processed/processedNetworkEmbLabelled.rds")

sVGeneExpress <- readRDS("processed/fullGeneExpressForNMF.rds")
geneExpressWithNet <- merge(sVGeneExpress, networkData, by = 0)
head(geneExpressWithNet[1:3,1:3])

dim(geneExpressWithNet)

rownames(geneExpressWithNet)<- geneExpressWithNet$Row.names; geneExpressWithNet$Row.names <- NULL
#negative so need to do offset
minVal <- abs(min(geneExpressWithNet))
minVal
geneExpressWithNetoffSet <- geneExpressWithNet+minVal
ranks = c(5,10,50,100,150,200,500)
noofruns <- 2
set.seed(1234)

# cl <- makeCluster(25, outfile="processed/nmfSvWithNetwork.txt")
# registerDoParallel(cl)
# # integratedDataPcaRes <- prcomp(fullGPCA, scale=T)
# print("running NMF")
resInteg.multiRank <- nmf(geneExpressWithNetoffSet, rank = ranks, nrun=noofruns, seed = 123456)


saveRDS(resInteg.multiRank,"processed/resInteg.multiRank.rds")
sVwithNetworkNMF <- readRDS("processed/resInteg.multiRank.rds")
nmfNetFitCl <- list()
WmatricesNet <- list()

for (rank in 1:length(sVwithNetworkNMF$measures$rank)){
  #print(rank)
  nmfNetFitCl[[rank]] <- sVwithNetworkNMF$fit[[rank]]
  WmatricesNet[[rank]] <- nmfNetFitCl[[rank]]@fit@W
}
matricesWithNet <- lapply(WmatricesNet, function(x){
  data.frame(list(x))
})


#need to get labelled genes...
getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}


sVwithNetworkLabelled <- lapply(matricesWithNet, getLabelledGenesFctn, labelsForNet)
saveRDS(sVwithNetworkLabelled, "processed/SVwithNetworkLabelled.rds")
###
# now i want to CONCATENATE network data
#first get NMF features:
nmfData <- readRDS("processed/unlabelledNMFDfs.rds")

concatNMFData <- function(dataList, netData, knownLabels){
  concat <- merge(dataList, netData, by=0);rownames(concat)<- concat$Row.names; concat$Row.names <- NULL
  concatL <- merge(concat, knownLabels, by=0);rownames(concatL)<- concatL$Row.names; concatL$Row.names <- NULL
  return(concatL)
}

concatList <- lapply(X = nmfData, FUN = concatNMFData, networkData, labelsForNet)
saveRDS(concatList, "processed/concatNMFsvNetwk.rds")
#machine learning moved to next script..