#integrate network data with gene expression data
library(tidyverse)
#might move this to a different script
library(tidymodels)
library(themis)
library(vip)
library(doParallel)
library(foreach)
library(NMF)


#embeddings created using pecanpy in python
embeddings <- read.table("processed/networkEdgeList.emb", sep = "\t", skip = 1)
colnames(embeddings)

#read in gene expression data:
# fullGeforPCA <- readRDS("processed/fullGeneExpressForPCA.rds")
# head(fullGeforPCA[1:5,1:3])
# #and labels:
# labelsFullDf <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
# head(labelsFullDf)
#importing node2vec result now
#split into 1000 dims because that is the dimensions used for pecanPy
readNetworkData <- function(embData, dimens){
  #dimens is the number of dimensions used for pecanpy
  embedding <- read.table(embData, sep = "\t", skip = 1)
  emb <- embedding %>%
    separate(colnames(embedding), c("hgnc_symbol", paste0("D", 1:dimens)), sep = " ")
  #rownames(emb) <- emb$D0; emb$D0 <- NULL
  return(emb)
}

#just separate because all has been put into one column
# embTest <- embeddings %>%
#   separate(colnames(embeddings), c("hgnc_symbol", paste0("D", 1:500)), sep = " ")

# embTest <- mutate_all(embTest, function(x) as.numeric(as.character(x)))
# rownames(emb) <- emb$D0


emb <- readNetworkData("processed/networkEdgeList.emb", dimens = 500)
#origEmb <-readNetworkData("processed/origNetworkEdgeList.emb", 1000)
#dim(origEmb)
dim(emb)
#head(origEmb[1:4,1:6])
#head(origEmb[,999:1001])
#now need to convert the genes to ensembl IDs 
netGenes <- data.frame(hgnc_symbol = emb$hgnc_symbol)
dim(netGenes)
networkGenes <- data.frame(hgnc_symbol = emb$hgnc_symbol)
write.table(netGenes, "processed/networkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)
#write.table(networkGenes,"processed/AllnetworkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)

##now read in the mapped file 

networkGenesWithEnsembl <- read.table("raw/allNetworkGenesMapped.txt", sep = "\t", header = T)
# networkGenesWithEnsembl <- networkGenesWithEnsembl %>% 
#   distinct(hgnc_symbol, .keep_all = T)
#embTest$hgnc_symbol <- embTest$D0
#rownames(networkGenesWithEnsembl)<- networkGenesWithEnsembl$hgnc_symbol
# netData_2 <- merge(networkGenesWithEnsembl, embTest, by="hgnc_symbol")
networkData <- merge(networkGenesWithEnsembl, emb, by="hgnc_symbol")
head(networkData[1:4,1:3])
rownames(networkData)<- networkData$ensembl_gene_id; networkData$hgnc_symbol<- NULL; networkData$ensembl_gene_id<-NULL
str(networkData$D2)
dim(networkData)
head(networkData[,488:500])

#now need to merge with gene exp
networkData <- mutate_all(networkData, function(x) as.numeric(as.character(x)))

sVGeneExpress <- readRDS("processed/fullGeneExpressForNMF.rds")
geneExpressWithNet <- merge(sVGeneExpress, networkData, by = 0)
head(geneExpressWithNet[1:3,1:3])

dim(geneExpressWithNet)

rownames(geneExpressWithNet)<- geneExpressWithNet$Row.names; geneExpressWithNet$Row.names <- NULL
minVal <- abs(min(geneExpressWithNet))
minVal
geneExpressWithNetoffSet <- geneExpressWithNet+minVal
ranks = c(5,10,50,100,150,200,500)
noofruns <- 2
set.seed(1234)

cl <- makeCluster(25, outfile="processed/nmfSvWithNetwork.txt")
registerDoParallel(cl)
# integratedDataPcaRes <- prcomp(fullGPCA, scale=T)
print("running NMF")
resInteg.multiRank <- nmf(geneExpressWithNetoffSet, rank = ranks, nrun=noofruns, seed = 123456)
#

saveRDS(resInteg.multiRank,"processed/resInteg.multiRank.rds")


# integPCs <- as.data.frame(integratedDataPcaRes$x)
newpcares <- prcomp(networkData, scale. = T)
integPCs <- newpcares$x
head(integPCs[1:5,1:4])
integ_100dim <- as.data.frame(integPCs)[,1:10]
# integ_50dim <- as.data.frame(integPCs)[,1:50]

# integlabelled_50dim <- merge(integ_50dim, labelsFullDf, by=0); integlabelled_50dim$Row.names<-NULL
integlabelled_100dim <- merge(integ_100dim, labelsFullDf, by=0); integ_100dim$Row.names<-NULL
rownames(integlabelled_100dim) <- integlabelled_100dim$Row.names; integlabelled_100dim$Row.names <-NULL
# #split data and create rec
integsplit100 <- split_processingData_fctn(integlabelled_100dim, 0.8)
# integsplit50<- split_processingData_fctn(integlabelled_50dim, 0.8)
# 
# #ml
# 
plan(multisession, workers=20)
mlresultWithPPI_100dim <- justToTestRF(integsplit100, algorithm = "PCA")
# mlresultWithPPI_50dim <- justToTestRF(integsplit50, algorithm = "PCA")
# 
# ##




