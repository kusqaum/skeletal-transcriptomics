#02_exlploratoryDataVisualisation
# library(tidyverse)
library(ggplot2)
library(preprocessCore)
library(NMF)
#library(doParallel)


####read in gene expression data
fullGeneExpress <- readRDS("processed/geneExpressForDimRed.rds")

mouseHumanWithLabs <- readRDS("processed/mouseHumanWithLabels.rds")
mouseHumanWithLabs$significant <- as.factor(mouseHumanWithLabs$significant)
labelsFullDf <- data.frame("significant"=mouseHumanWithLabs$significant, row.names = rownames(mouseHumanWithLabs))
str(labelsFullDf$significant)
labelsFullDf$significant <- as.factor(labelsFullDf$significant)
write.table(labelsFullDf, "processed/labelsFullDf.txt", row.names = T, sep = "\t", quote = F)
#remove rows summing to 0
fullGeneExpress <- fullGeneExpress[rowSums(fullGeneExpress)>0,] 
#here there are no 0s but just incase this does happen, remember their rownames
filtRows <- rownames(fullGeneExpress)
####log transform and normalise####

head(fullGeneExpress[1:5,1:5])
min(fullGeneExpress)#need to to pseudocount
log2(90.19820+1)
# 
fullGeneExpressL <- log2(fullGeneExpress+1)
min(fullGeneExpressL)
head(fullGeneExpressL[1:5,1:5])

# now let's normalise quantiles
fullGeneExpressN <- normalize.quantiles(as.matrix(fullGeneExpressL))
min(fullGeneExpressN)
head(fullGeneExpressN[1:5,1:5])
row.names(fullGeneExpressN)<- filtRows
#we don't really care about colnames anyways
fullGeneExpressforNMF <- as.data.frame(fullGeneExpressN)

# lowestValue <- abs(min(fullGeneExpressforNMF)) # there are no negative values now

#these are what will be used for NMF and PCA
#fullGeneExpressforNMF[1,1]+lowestValue
#fullGeneExpressforNMF <- fullGeneExpressforNMF+lowestValue
head(fullGeneExpressforNMF[1:5,1:5])
fullGeneExpressforPCA <- as.data.frame(fullGeneExpressN)
head(fullGeneExpressforPCA[1:5,1:5])
saveRDS(fullGeneExpressforPCA, "processed/fullGeneExpressForPCA.rds")
saveRDS(fullGeneExpressforNMF, "processed/fullGeneExpressForNMF.rds")
# sum(smallGeneExpress2$significant == "TRUE")
# sum(smallGeneExpress2$significant == "FALSE")
#now we have 1300 and 6578 positive and negative genes

#nmfSeed()
####perform NMF####
#determine the sequence
ranks = c(5,10,50,100,150,200,500)
# 
#and the number of ranks
noofruns <- 2
#set seed as well
set.seed(1234)

print("running NMF")
res.multiRank <- nmf(fullGeneExpressforNMF, rank = ranks, nrun=noofruns, seed = 123456)
saveRDS(res.multiRank, "processed/res.multiRank.rds")

#look at performance measures of the factorisation
summary(res.multiRank)
# consensusmap(res.multiRank, labCol = NA, labRow = 1)

#plot(res.multiRank)

res.multiRank$measures$rank


#basis_matrix <- list()
nmfFitCl <- list()
Wmatrices <- list()

for (rank in 1:length(res.multiRank$measures$rank)){
  #print(rank)
  nmfFitCl[[rank]] <- res.multiRank$fit[[rank]]
  Wmatrices[[rank]] <- nmfFitCl[[rank]]@fit@W
  #print(nmfFitCl[[rank]])
  #basis_matrix <- append(basis_matrix,nmfFitCl[[rank]]@fit@W)
}
basis_matrices <- lapply(Wmatrices, function(w){
  data.frame(list(w))
})

saveRDS(basis_matrices, "processed/unlabelledNMFDfs.rds")
#here, from NMF reduction result I will get up to only the genes that have labels from each matrix:
# so that can be used for machine learning:
# labelsFullDf[,1]<- as.factor(labelsFullDf[,1])
getLabelledGenesFctn <- function(matrixList, knownLabels){
  knownLabels[,1] <- as.factor(knownLabels[,1])
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}
#this is for NMF machine learning -cool
labelledGenesNMFRes <- map(.x = Wmatrices, .f = getLabelledGenesFctn, labelsFullDf)
#dim to check merged correctly just in case
dim(labelledGenesNMFRes[[1]])
head(labelledGenesNMFRes[[1]])
labelledGenesNMFResTest <- labelledGenesNMFRes


####perform PCA####

#perform PCA. each of our data points needs to be a gene so not going to transpose
pca_res <- prcomp(fullGeneExpressforPCA, scale = T)
pcs <- data.frame(pca_res$x)
dim(pcs)
pcsLabelled <- merge(pcs, labelsFullDf, by=0)
dim(pcsLabelled)
row.names(pcsLabelled) <- pcsLabelled$Row.names; pcsLabelled$Row.names <- NULL
dim(pcsLabelled)

varexplained <-summary(pca_res)$importance[2,]*100
ggplot(pcsLabelled, aes(x=PC1, y=PC2))+
  geom_point(aes(col = significant))+
  xlab(paste("PC1", round(varexplained[1],1), "%"))+
  ylab(paste("PC2", round(varexplained[2], 1), "%"))+
  theme_bw(base_size = 24)

tmp <- t.test(y=as.logical(pcsLabelled$significant), x=pcsLabelled[,518])
tmp$p.value
tmp$statistic
t.test(pcsLabelled$PC2 ~ mouseHumanWithLabs$significant)$p.value

####do pca on shuffled data and then let's do machine learning
# shuffledPCA<- randomize(fullGeneExpressforPCA); rownames(shuffledPCA)<- rownames(fullGeneExpressforPCA)
# pcaRes_shuffled <- prcomp(shuffledPCA, scale. = T)
# shuffPCs <- pcaRes_shuffled$x
# #just test on dim 50 or 100
# shuffPCsdim50 <- shuffPCs[,1:50]
# #add on labels:
# shuffPCsdim50_L <- merge(shuffPCsdim50, labelsFullDf, by=0)
# rownames(shuffPCsdim50_L)<- shuffPCsdim50_L$Row.names ;shuffPCsdim50_L$Row.names<- NULL
# head(shuffPCsdim50_L)

#split and create rec- this was just temporary so I need to run it in full in the model training where 
#the functions for machine learning are
# shuff_Split <- split_processingData_fctn(shuffPCsdim50_L, proportion = 0.8)
# shuff_MLRes <- justToTestRF(shuff_Split, algorithm = "PCA")

removeLabelsFctn <- function(dataWithLabels){
  dataWithLabels[,ncol(dataWithLabels)] <- NULL
  return(dataWithLabels)
}

#also create a scree plot for pca
var_explained <- (pca_res$sdev)^2/sum((pca_res$sdev)^2)
varExpl_df <- data.frame(PC = 1:length(var_explained),
                         varianceExpl = (var_explained),
                         cumulativeVariance = cumsum(var_explained))

ggplot(varExpl_df[1:200,], aes(x = PC, y = 100*cumulativeVariance)) +
  geom_bar(stat = 'identity', col = 'black', fill="lightgrey") +
  #geom_point()+
  #geom_line()+
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 90, col = "red", lty = "dashed")+
  #scale_x_continuous(breaks = seq(1, nrow(varExpl_df[1:200,]), 1)) +
  xlab("PC index") + ylab("% explained variance by principal components")

  #xlim(0, 10)


##### changes for PCA dimension reduction here: #####


#loop through length
#here creating a list for all dimensions when applied PCA
for (i in 1:length(ranks)) {
  #create empty lists
  pcList <- list()
  #pcList2 <- list()
  #loop through each of the dimensions
  for (l in (ranks)) {
    #grab each dimension
    dime <- pcsLabelled[,1:l]
    #dime2<- pcs[,1:l]
    pcList <- append(pcList, list(dime))
    #pcList2<- append(pcList2, list(dime2))
  }
}

#do some more in depth tuning
# for (j in 1:length(tuneRanks)) {
#   #create empty lists
#   pcList_tuned <- list()
#   #loop through each of the dimensions
#   for (m in (tuneRanks)) {
#     dime2<- pcsLabelled[,1:m]
#     pcList_tuned<- append(pcList_tuned, list(dime2))
#   }
# }
# pcList is also for PCA machine learning-

##


####changes for PCA t-testing here####-update made a new fctn
#make a function to take in the list of dfs created above
# pca_function <- function(listPCAres, labs){
#   #dataframe <- list()
#   pvalList <- numeric()
#   for (pc in 1:ncol(listPCAres)){
#     pvalueResult <- t.test(y=labs, x=(listPCAres[,pc]))$p.value
#     pvalList <-c(pvalList, pvalueResult)
#   }
#   minPval <- min(pvalList)
#   mindim <- which.min(pvalList)
#   dataframe <- data.frame(pVals = minPval,
#                           Dimension = ncol(listPCAres),
#                           Feature=mindim,
#                           Algorithm = "PCA")
# }

#####just do all t-tests here ####


#### perform t-test on NMf result#### ignore this for now because have made a function to do it for both pca and nmf
#(all dimensions and obtain smallest p-value on each one)
# nmfMinPvals <- list()
#
# length(Wmatrices)
# #loop through the no. of Ws not the actual matrices
# for (mt in 1:length(Wmatrices)) {
#   #get ready list
#   pValues <- numeric()
#   #loop through each col
#   for (dim in 1:ncol(labelledGenesNMFRes[[mt]])-1) {
#     pResultNMF <- t.test(y= mouseHumanWithLabs$significant, x=labelledGenesNMFRes[[mt]][, dim])$p.value
#     boxplot(Wmatrices[[mt]][,dim]~smallGeneExpress2$significant, xlab="")
#     #nmfPvals <- c(nmfPvals, pResultNMF)
#     pValues <- c(pValues, pResultNMF)
#   }
#   #get min p val and their positions and make into list as DFs
#   minP <- min(pValues)
#   minPdim <- which.min(pValues)
#   nmfMinPvals[[mt]] <- data.frame(pVals = minP,
#                               Dimension = ncol(Wmatrices[[mt]]),
#                               Feature=minPdim,
#                               Algorithm = "NMF")
#
# }


# checking a single to-test p-value result here
tempor <- data.frame(Wmatrices[[1]])
t.test(y=as.logical(mouseHumanWithLabs$significant), x=tempor$X3)$p.value
t.test(y=as.logical(mouseHumanWithLabs$significant), x=tempor$X3)$statistic
t.test(y=as.logical(mouseHumanWithLabs$significant), x=tempor$X5)$statistic
t.test(y=as.logical(mouseHumanWithLabs$significant), x=basis_matrices[[2]]$X5)$statistic
t.test(y=as.logical(mouseHumanWithLabs$significant), x=basis_matrices[[2]]$X5)$p.value

#


performTtestFctn = function(listRes, labels, algorithm){ 
  if(algorithm == "PCA"){
    pvalList <- numeric()
    for (c in 1:ncol(listRes)){
      # df <- as.data.frame(merge(listRes[,c], labels, by = 0)); df$Row.names <-NULL
      # print(head(df))
      # df_T <- df %>% filter(y=="TRUE")
      # df_t <- c(df_T$x)
      # df_F <- df%>% filter(y=="FALSE")
      # df_f <- c(df_F$x)
      # 
      # pvalueResult <- t.test(df_f, df_t, var.equal=F)$p.value
      pvalueResult <- t.test(y=as.logical(labels) , x=listRes[,c])$p.value
      pvalList <-c(pvalList, pvalueResult)
    }
    minPval <- min(pvalList)
    mindim <- which.min(pvalList)
    dataframe <- data.frame(pVals = minPval,
                            Dimension = ncol(listRes),
                            Feature=mindim,
                            Algorithm = algorithm)

  }
  else if(algorithm == "NMF"){

    pvalList <- numeric()
    for (c in 1:ncol(listRes)){
      pvalueResult <- t.test((listRes[,c]) ~ labels)$p.value
      pvalList <-c(pvalList, pvalueResult)
    }
    minPval <- min(pvalList)
    mindim <- which.min(pvalList)
    dataframe <- data.frame(pVals = minPval,
                            Dimension = ncol(listRes),
                            Feature=mindim,
                            Algorithm = algorithm)

  }
}

# labelledGenesNMFRes <- lapply(basis_matrices)
unlabelledNMFgenes <- lapply(labelledGenesNMFRes, removeLabelsFctn)
pcaDimPValues <- lapply(X = pcList, FUN = performTtestFctn, labels = labelsFullDf$significant, algorithm="PCA")
 # pcaDimPValues2 <- lapply(X=pcList2, FUN = performTtestFctn, labels = mouseHumanWithLabs$significant, algorithm="PCA")
nmfMinPvals <- lapply(X=unlabelledNMFgenes, FUN = performTtestFctn, labels = mouseHumanWithLabs$significant, algorithm="NMF")
t.test(labelledGenesNMFRes[[5]][,4] ~ labelledGenesNMFRes[[5]][,ncol(labelledGenesNMFRes[[5]])])$p.value
ncol(labelledGenesNMFRes[[1]])



# genesNMFres <- lapply(labelledGenesNMFRes,removeLabelsFctn)
# resultofnmf <- lapply(genesNMFres, pca_function, labs = mouseHumanWithLabs$significant)


# pcaDfs <- lapply(pcList, FUN=pca_function, mouseHumanWithLabs$significant)
# pvalsPcaDf <- do.call(rbind, pcaDfs)
pvalsPcaDf<- do.call("rbind", pcaDimPValues)
ggplot(pvalsPcaDf, aes(x = Dimension, y=-log10(pVals), col = Algorithm)) +
  geom_point()+
  theme_bw(base_size = 18) +
  xlab("k dimensions") + ylab("-log10 P-Value")

pcaStrongFeat <- pvalsPcaDf[which.min(pvalsPcaDf$pVals),]
pcaFeature <- pvalsPcaDf[which.min(pvalsPcaDf$pVals),]$Feature
pcaDim <-pvalsPcaDf[which.min(pvalsPcaDf$pVals),]$Dimension
boxplot(pcList[[which.min(pvalsPcaDf$pVals)]][,pcaStrongFeat$Feature]
        ~labelsFullDf$significant)
pcaFeatForBoxplot <- data.frame(Feature = pcList[[which.min(pvalsPcaDf$pVals)]][,pcaStrongFeat$Feature],
                     labels = labelsFullDf$significant)
pcaFeatForBoxplot <- pcaFeatForBoxplot %>%mutate(Association = case_when(labels=="TRUE"~"Associated", labels=="FALSE"~"Not associated"))

ggplot(pcaFeatForBoxplot, aes(x=Association, y=Feature))+
  geom_boxplot()+
  ylab(paste0("PCA Feature ", pcaFeature, " (", "k =", pcaDim, ")", sep = " "))+
  xlab("")
####boxplot for feature giving strongest signal NMF####
#rbind list of dfs
NMFfeaturesMinP <- do.call("rbind", nmfMinPvals)
#boxplot for the feature with the most strongest signal:
nmfStrongFeat <- NMFfeaturesMinP[which.min(NMFfeaturesMinP$pVals),]
nmfFeature <- NMFfeaturesMinP[which.min(NMFfeaturesMinP$pVals),]$Feature
nmfDim <- NMFfeaturesMinP[which.min(NMFfeaturesMinP$pVals),]$Dimension
nmfFeatForBoxplot <- data.frame(Feature = labelledGenesNMFRes[[which.min(NMFfeaturesMinP$pVals)]][,nmfStrongFeat$Feature],
                                labels = labelsFullDf$significant)
nmfFeatForBoxplot <- nmfFeatForBoxplot %>%mutate(Association = case_when(labels=="TRUE"~"Associated", labels=="FALSE"~"Not associated"))
ggplot(nmfFeatForBoxplot, aes(x=Association, y=Feature))+
  geom_boxplot()+
  ylab(paste0("NMF Feature ", nmfFeature, " (", "k =", nmfDim, ")", sep = " "))+
  xlab("")
boxplot(labelledGenesNMFRes[[which.min(NMFfeaturesMinP$pVals)]][,nmfStrongFeat$Feature]
~mouseHumanWithLabs$significant)


#plot result
ggplot(NMFfeaturesMinP, aes(x = Dimension, y=-log10(pVals), col = Algorithm)) +
  geom_point()+
  theme_bw(base_size = 18) +
  xlab("k dimensions") + ylab("-log10 P-Value")


###plot for both algorithms on one page####
# pValsdf2 <- pValsdf
# colnames(pValsdf2)[1]<- "pVals"
allDimensionReduction <- rbind(NMFfeaturesMinP, pvalsPcaDf)
ggplot(allDimensionReduction, aes(x = Dimension, y=-log10(pVals), size=Algorithm, col=Algorithm, shape=Algorithm)) +
  geom_point() +
  theme_bw(base_size = 20)+
  xlab("k dimensions") +ylab(expression("-log"[10]* "(adj " *italic(P)* "-value)"))+
  scale_colour_manual(values = c("#91D1C2B2", "#E64B35B2"))+
  scale_size_manual(values = c(7, 5))
  # scale_x_continuous(breaks = seq(2,50, 2))

#now save this plot



labelledGenesPCARes <- lapply(pcList, FUN = getLabelledGenesFctn, labelsFullDf)
# labelledGenesNMFRes

# read.table("pro") - need to read in the new labels
# differentlyLabelledPCA <- map(.x = pcList, .f = getLabelledGenesFctn, newImpcGenes)


saveRDS(labelledGenesPCARes, "processed/pcaDataframes.rds")
saveRDS(labelledGenesNMFRes, "processed/nmfDataframes.rds")

##save result for model training
# saveRDS(basis_matrices, "processed/nmf_wMatrices.rds")


####i want to get all the unlabelled genes so can use the model to rank them later####
fullGeneExpressN <- as.data.frame(fullGeneExpressN)
noLabels <- row.names(fullGeneExpressN)[!row.names(fullGeneExpressN) %in%
                                          row.names(mouseHumanWithLabs)]
unlabelled <- fullGeneExpressN[noLabels,]
saveRDS(unlabelled, "processed/unStudiedGenes.rds")

#pc500 <- pcsLabelled[,1:500]
#pc500_Labelled <- merge(pc500, labelsFullDf, by=0);pc500_Labelled$Row.names<-NULL
# if you are going to call in the functions either do it in the script where the functions are or write out 
#the functions above
# split_pc500 <- split_processingData_fctn(pc500_Labelled, 0.8)
# pca500dimMLres <- justToTestRF(split_pc500, algorithm="PCA")


#--------------------------------------------------------------------------
