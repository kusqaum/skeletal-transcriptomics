#02_exlploratoryDataVisualisation
library(tidyverse)
library(NMF)


geneExpressionWithLabels <- read.table("processed/geneExpressionDataWithLabels.txt", header = T, stringsAsFactors = F)
#just get the first 50 samples as a subset
smallGeneExpress <- geneExpressionWithLabels %>% select(1:50)
#add on the labels
smallGeneExpress$significant <- geneExpressionWithLabels$significant

#remove genes(rows) that have a sum of 0 abundance
smallGeneExpress <- smallGeneExpress %>% filter(rowSums(across(where(is.numeric)))!=0)
#change to matrix
matrix <- as.matrix(smallGeneExpress)
nmfSeed()

#determine the sequence
noRanks <- seq(2,10,2)
#and the number of ranks
noofruns <- 2
#set seed as well
res.multiRank <- nmf(matrix[,1:50], rank = noRanks, nrun=noofruns, seed = 123456)
#can choose to plo just cophenetic + silhouette
plot(res.multiRank, 'cophenetic', 'silhouette')
#res.multiRank2 <- nmf(expressionDataMat[,1:1128], seq(2,200,5), seed = 123456)

#look at performance measures of the factorisation
summary(res.multiRank)
consensusmap(res.multiRank, labCol = NA, labRow = NA)

res.multi.method <- nmf(smallGeneExpress[,1:50], 2, seed =123456, list("brunet","lee", "ns"), .options= "t")
#compare(res.multi.method)

#maxbrunCoph <- max(res.multi.method$brunet$measures$cophenetic)
#maxNsCoph <- max(res.multi.method$nsNMF$measures$cophenetic)

#check sumary measures for each method
compare(res.multi.method)
#look at error tracks for each method
plot(res.multi.method)
plot(res.multiRank)


res.multiRank$measures
# choosing the best rank
# so when thte cophenetic correlation coefficient begins to fall
max_cophenetic <- max(res.multiRank$measures$cophenetic)
# loop through the position of each rank.
#no. of ranks we used:
length(res.multiRank$measures$rank)
for (i in 1:length(res.multiRank$measures$rank)){
  if (res.multiRank$measures$cophenetic[i]== max_cophenetic){
    idxCophenMax <- i #idxCophen is the position of the rank.
  }
}
idxCophenMax
res.multiRank$measures$cophenetic
#let's look at the silhouette:
res.multiRank$measures$silhouette.consensus

nmfFitClass  <- res.multiRank$fit[[idxCophenMax]]
#get the w and h matrices for the rank with highest cophenetic
h_matrix <- nmfFitClass@fit@H
w_matrix <- nmfFitClass@fit@W
###ggplot(as.data.frame(w_matrix), aes(x=V1, y=V2, col=smallGeneExpress$significant))+geom_point()
consensusmap(res.multiRank$fit[[idxCophenMax]])
consensusmap(res.multiRank)
res <- nmf(smallGeneExpress[,1:50], 2)
nmf.w <- basis(res)
ggplot(as.data.frame(h_matrix))


#perform PCA. each of our data points needs to be a gene so not going to transpose
pca_res <- prcomp((smallGeneExpress[,1:50]), scale. = T)
#look at variance
summary(pca_res)$importance
var_explained <- summary(pca_res)$importance[2,]*100
#get 1&2nd PCs to plot
pca_res$x[1:5,1:2]
ggplot(as.data.frame(pca_res$x[,1:2])) + 
  aes(x=PC1, y=PC2) + 
  geom_point(aes(col=smallGeneExpress[,51])) +
  xlab(paste("PC1", round(var_explained[1],1), "%")) +
  ylab(paste("PC2", round(var_explained[2],1), "%"))

0.710630 + 0.104510 +0.056000 
