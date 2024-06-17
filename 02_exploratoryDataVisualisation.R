#02_exlploratoryDataVisualisation
library(tidyverse)
library(preprocessCore)
library(NMF)


# geneExpressionWithLabels <- read.table("processed/geneExpressionDataWithLabels.txt", header = T, stringsAsFactors = F)
geneExpressionWithLabels <- readRDS("processed/mouseHumanWithLabels.rds")
labels <- data.frame(geneExpressionWithLabels$significant)
row.names(labels) <- rownames(geneExpressionWithLabels)
colnames(labels)[1] <- "significant"

#just get the first 50 samples as a subset
smallGeneExpress <- geneExpressionWithLabels %>% select(1:50)

#remove rows that sum to 0
smallGeneExpress <- smallGeneExpress %>% filter(rowSums(across(where(is.numeric)))!=0)

filteredRows <- rownames(smallGeneExpress)
#change 0 values to 1: 
smallGeneExpress[smallGeneExpress == 0] <- 1
smallGeneExpress <- log2(smallGeneExpress)

#lowestVal <-abs(min(smallGeneExpress))
# add on this value to each value
#smallGeneExpress <- smallGeneExpress+lowestVal
#now let's normalise quantiles
smallGeneExpress <- normalize.quantiles(as.matrix(smallGeneExpress))
row.names(smallGeneExpress) <- filteredRows
#we don't really care about colnames anyways
# colnames(smallGeneExpress) <- colnames(geneExpressionWithLabels[,1:50])
smallGeneExpress <- as.data.frame(smallGeneExpress)
lowestVal <-abs(min(smallGeneExpress))
smallGeneExpress <- smallGeneExpress+lowestVal



#smallGeneExpress2 <- geneExpressionWithLabels[,1:50]
#smallGeneExpress2 <- normalize.quantiles(as.matrix(smallGeneExpress2))
#row.names(smallGeneExpress2) <- rownames(geneExpressionWithLabels)
#colnames(smallGeneExpress2) <- colnames(geneExpressionWithLabels[,1:50])
#smallGeneExpress2 <- as.data.frame(smallGeneExpress2)

#dim(smallGeneExpress3) 

#log2 transform the data
#smallGeneExpress3 <- log2(smallGeneExpress3)
#hist(smallGeneExpress3$t0h_3)
#merge back with the labels now, by row:
smallGeneExpress2 <- merge(smallGeneExpress, labels, by = 0)
row.names(smallGeneExpress2)<- smallGeneExpress2$Row.names; smallGeneExpress2$Row.names <- NULL
sum(smallGeneExpress2$significant == "TRUE")
sum(smallGeneExpress2$significant == "FALSE")
#now we have 1300 and 6578 positive and negative genes


#nmfSeed()

#determine the sequence
noRanks <- seq(2,10,2)
#and the number of ranks
noofruns <- 2
#set seed as well
res.multiRank <- nmf(smallGeneExpress[,1:50], rank = noRanks, nrun=noofruns, seed = 123456)


#nmf without logg transf
#res.multiRank2 <- nmf(smallGeneExpress[,1:50], rank = noRanks, nrun=noofruns, seed = 123456)

#look at performance measures of the factorisation
summary(res.multiRank)
consensusmap(res.multiRank, labCol = NA, labRow = NA)

#res.multi.method <- nmf(smallGeneExpress[,1:50], 2, seed =123456, list("brunet","lee", "ns"), .options= "t")
#compare(res.multi.method)

#maxbrunCoph <- max(res.multi.method$brunet$measures$cophenetic)
#maxNsCoph <- max(res.multi.method$nsNMF$measures$cophenetic)

#check sumary measures for each method
#compare(res.multi.method)
#look at error tracks for each method
#plot(res.multi.method)
#plot(res.multiRank)


res.multiRank$measures$rank

#idxCophenMax
#res.multiRank$measures$cophenetic
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

nmfMinPvals <- list()
length(Wmatrices)
#loop through the no. of Ws not the actual matrices
for (mt in 1:length(Wmatrices)) {
  #get ready list
  pValues <- numeric()
  #loop through each col
  for (dim in 1:ncol(Wmatrices[[mt]])) {
    pResultNMF <- t.test(Wmatrices[[mt]][, dim] ~ smallGeneExpress2$significant)$p.value
    boxplot(Wmatrices[[mt]][,dim]~smallGeneExpress2$significant, xlab="")
    #nmfPvals <- c(nmfPvals, pResultNMF)
    pValues <- c(pValues, pResultNMF)
  }
  #get min p val and their positions and make into list as DFs
  minP <- min(pValues)
  minPdim <- which.min(pValues)
  nmfMinPvals[[mt]] <- data.frame(pVals = minP,
                              Dimension = ncol(Wmatrices[[mt]]),
                              Feature=minPdim, 
                              Algorithm = "NMF")
  
}


#rbind list of dfs
NMFfeaturesMinP <- do.call("rbind", nmfMinPvals)
#boxplot for the feature with the most strongest signal:
nmfStrongFeat <- NMFfeaturesMinP[which.min(NMFfeaturesMinP$pVals),]
boxplot(Wmatrices[[which.min(NMFfeaturesMinP$pVals)]][,nmfStrongFeat$Feature] 
        ~smallGeneExpress2$significant)

##########

#make the current df have the dimensions as rownames:
tempDF <- NMFfeaturesMinP
rownames(tempDF) <- tempDF$Dimension

tempDF[,min(tempDF$pVals)]
ggplot(NMFfeaturesMinP, aes(x = Dimension, y=-log10(pVals), col = Algorithm)) +
  geom_point()+
  theme_bw(base_size = 18) +
  xlab("k dimensions") + ylab("-log10 P-Value")


# checking a single to-test p-value result here is consistent with table created above.
# tempor <- data.frame(Wmatrices[[1]])
# t.test(tempor[,2] ~ smallGeneExpress2$significant)$p.value
# t.test(tempor[,2]~ smallGeneExpress2[,51])$p.value
# 


####PCA####

#perform PCA. each of our data points needs to be a gene so not going to transpose
pca_res <- prcomp((smallGeneExpress2[,1:50]), scale =T)
pcs <- data.frame(pca_res$x)

tmp <- t.test(pcs[,"PC3"] ~smallGeneExpress2$significant)
tmp$p.value
tmp$statistic
#empty vector to store pvalue results:
#t_testPvaluesPC <- numeric()
#so I am going to loop through all the principal components: minus 1 b/c don't want to include last column
#for (pc in 1:(ncol(pcs)-1)) {
  #print(pc)
  #extract the result of t-test for each
 # tresult[[pc]] <- t.test(pcs[[pc]] ~ significant, data = pcs)
  #obtain the p-value for each
  #pVal[[pc]] <- tresult[[pc]]$p.value
  #then add it to empty numeric vector created above
  #t_testPvaluesPC <- c(t_testPvaluesPC, pVal[[pc]])
#}


pValsPCA <- numeric()
for (PC in 1:ncol(pcs)) {
  resultfortest <- t.test(pcs[,PC]~ smallGeneExpress2$significant)$p.value
  pValsPCA <- c(pValsPCA, resultfortest)
  #resultfortest_pval <- resultfortest$
  #smallestPval <- which.min(pValsPCA)
}
pValsdf <- as.data.frame(pValsPCA); pValsdf$Dimension<- 1:ncol(pcs)

###colnames(pValsdf)[1] <- "pVals"
pValsdf$Feature <- c(1:ncol(pcs))
pValsdf$Algorithm <- c('PCA')
#t.test(smallGeneExpress3$significant, pcs$PC50)$p.value
ggplot(pValsdf, aes(x=Dimension, y=-log10(pValsPCA), col = Algorithm)) +
  geom_point()+
  theme_bw(base_size = 18) +
  scale_x_continuous(breaks = seq(2,50, 2)) +
  xlab("k dimensions") + ylab("-log10 P-Value")

#also create a scree plot for pca
var_explained <- (pca_res$sdev)^2/sum((pca_res$sdev)^2)
varExpl_df <- data.frame(PC = 1:length(var_explained),
                         varianceExpl = (var_explained),
                         cumulativeVariance = cumsum(var_explained))

ggplot(varExpl_df[1:10,], aes(x = PC, y = 100*cumulativeVariance)) +
  geom_bar(stat = 'identity', col = 'black', fill="lightgrey") + 
  #geom_point()+
  #geom_line()+
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 95, col = "red", lty = "dashed")+
  scale_x_continuous(breaks = seq(1, nrow(varExpl_df[1:8,]), 1)) +
  xlab("PC index") + ylab("% explained variance by principal components")
  
  #xlim(0, 10)
####plot for both algorithms on one page####
pValsdf2 <- pValsdf
colnames(pValsdf2)[1]<- "pVals"
allDimensionReduction <- rbind(NMFfeaturesMinP, pValsdf2)
ggplot(allDimensionReduction, aes(x = Dimension, y=-log10(pVals), col=Algorithm)) +
  geom_point() +
  theme_bw(base_size = 20)+
  xlab("k dimensions") +ylab("-log10 P-Value")+
  scale_x_continuous(breaks = seq(2,50, 2))
  

##save result for model training
saveRDS(basis_matrices, "processed/nmf_wMatrices.rds")
saveRDS(smallGeneExpress2, "processed/geneExpression.rds")


##### changes for PCA here: #####
pcaSeq <- seq(2,10,2)


# reduceFctn <- function(pcDf){
#   for (n in 1:length(pcDf)){  
#     listDFS <- list()
#     for (j in pcaSeq) {
#       reduced <- pcDf[,1:j]
#       # listDFS[[n]] <- reduced
#       listDFS <- append(listDFS, list(reduced))
#     }
#   }
# }

#   listofdfs <- list()
#   theSeq <- seq(2,10,2)
#   for (j in theSeq){
#     reduced <- pcDf[,1:i]
#     listofdfs[[j]] <- reduced
#   }
#   
# }
# testingAgaintimes2 <- apply(pcs,2, FUN = reduceFctn)

#make the seq the same as nmf
pcaSeq <- seq(2,10,2)
#loop through length
#here creating a list for all dimensions when applied PCA
for (i in 1:length(pcaSeq)) {
  #create empty lists
  pcList <- list()

  # testingagain <- data.frame(pcs[,1:i])
  # testingagain <- data.frame(pcs[,i])
  #loop through each of the dimensions
  for (l in (pcaSeq)) {
    #grab each dimension
    dime <- pcs[,1:l]
    pcList <- append(pcList, list(dime))
    # testing[[i]] <- (testingagain)
  }
}

#make a function to take in the list of dfs created above
pca_function <- function(listPCAres){
  #dataframe <- list()
  pvalList <- numeric()
  for (pc in 1:ncol(listPCAres)){
    pvalueResult <- t.test(listPCAres[,pc]~smallGeneExpress2$significant)$p.value
    pvalList <-c(pvalList, pvalueResult)
  }
  minPval <- min(pvalList)
  mindim <- which.min(pvalList)
  dataframe <- data.frame(pVals = minPval,
                          Dimension = ncol(listPCAres),
                          Feature=mindim, 
                          Algorithm = "PCA")
}

#########test fccctn on 1 dataset
test <-pcList[[1]]
t <- numeric()
for (v in 1:ncol(test)) {
  
  pres <- t.test(test[,v]~smallGeneExpress2$significant)$p.value
  t <- c(t, pres)
}
min(t)
which.min(t)
dataf <- data.frame(value = min(t), dimen = ncol(test), f = which.min(t))
###############

pcaDfs <- lapply(pcList, FUN=pca_function)
pvalsPcaDf <- do.call(rbind, pcaDfs)
ggplot(pvalsPcaDf, aes(x = Dimension, y=-log10(pVals), col = Algorithm)) +
  geom_point()+
  theme_bw(base_size = 18) +
  xlab("k dimensions") + ylab("-log10 P-Value")


# # compressed_PCA <- pcs
# compressed_PCA$significant <- smallGeneExpress3$significant
# str(compressed_PCA)
# compressed_PCA$significant <- as.factor(compressed_PCA$significant)
