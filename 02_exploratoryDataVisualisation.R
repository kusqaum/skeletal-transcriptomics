#02_exlploratoryDataVisualisation
library(tidyverse)
library(preprocessCore)
library(NMF)


geneExpressionWithLabels <- read.table("processed/geneExpressionDataWithLabels.txt", header = T, stringsAsFactors = F)
labels <- data.frame(geneExpressionWithLabels$significant)
row.names(labels) <- rownames(geneExpressionWithLabels)
colnames(labels)[1] <- "significant"

#just get the first 50 samples as a subset
smallGeneExpress <- geneExpressionWithLabels %>% select(1:50)
smallGeneExpress2 <- geneExpressionWithLabels[,1:50]
smallGeneExpress2 <- normalize.quantiles(as.matrix(smallGeneExpress2))
row.names(smallGeneExpress2) <- rownames(geneExpressionWithLabels)
colnames(smallGeneExpress2) <- colnames(geneExpressionWithLabels[,1:50])
smallGeneExpress2 <- as.data.frame(smallGeneExpress2)

smallGeneExpress3 <- smallGeneExpress2[apply(smallGeneExpress2!=0,1,all),]

#log2 transform the data
smallGeneExpress3 <- log2(smallGeneExpress3)

#merge back with the labels now, by row:
smallGeneExpress3 <- merge(smallGeneExpress3, labels, by = 0)
row.names(smallGeneExpress3)<- smallGeneExpress3$Row.names; smallGeneExpress3$Row.names <- NULL
sum(smallGeneExpress3$significant == "TRUE")
sum(smallGeneExpress3$significant == "FALSE")
#now we have 944 and 4426 positive and negative genes



smallGeneExpressMatrix <- as.matrix(smallGeneExpress)
smallGeneExpressMatrix <- normalize.quantiles(smallGeneExpressMatrix[,1:50])
#smallGeneExpressMatrix$significant <- smallGeneExpress$significant
row.names(smallGeneExpressMatrix) <- rownames(smallGeneExpress)
colnames(smallGeneExpressMatrix) <- colnames(smallGeneExpress[,1:50])
#remove rows that have value of 0:
smallGeneExpress_2 <- as.data.frame(smallGeneExpressMatrix)
smallGeneExpress_2$significant <- smallGeneExpress$significant

smallGeneExpress_2 <- smallGeneExpress_2[apply(smallGeneExpress_2!=0,1,all),]
smallGeneExpress_2$significant <- smallGeneExpress$significant
#log transform
smallGeneExpress_2 <- log2(smallGeneExpress_2)

#add on the labels
smallGeneExpress_2$significant <- smallGeneExpress$significant
#remove genes(rows) that have a sum of 0 abundance
smallGeneExpress <- smallGeneExpress %>% filter(rowSums(across(where(is.numeric)))!=0)

#change to matrix
matrix <- as.matrix(smallGeneExpress)
#nmfSeed()

#determine the sequence
noRanks <- seq(2,10,2)
#and the number of ranks
noofruns <- 2
#set seed as well
res.multiRank <- nmf(abs(smallGeneExpress3[,1:50]), rank = noRanks, nrun=noofruns, seed = 123456)
#can choose to plo just cophenetic + silhouette
plot(res.multiRank, 'cophenetic', 'silhouette')

#look at performance measures of the factorisation
summary(res.multiRank)
consensusmap(res.multiRank, labCol = NA, labRow = NA)

#res.multi.method <- nmf(smallGeneExpress[,1:50], 2, seed =123456, list("brunet","lee", "ns"), .options= "t")
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
nmfFitClass2 <- res.multiRank$fit[[1]]
head(nmfFitClass2@fit@W)
#get the w and h matrices for the rank with highest cophenetic
h_matrix <- nmfFitClass@fit@H
w_matrix <- nmfFitClass@fit@W #this one has metagenes as the cols. and gene as rows
###ggplot(as.data.frame(w_matrix), aes(x=V1, y=V2, col=smallGeneExpress$significant))+geom_point()
consensusmap(res.multiRank$fit[[idxCophenMax]])
consensusmap(res.multiRank)
res <- nmf(smallGeneExpress[,1:50], 2)
nmf.w <- basis(res)
ggplot(as.data.frame(h_matrix))


#perform PCA. each of our data points needs to be a gene so not going to transpose
pca_res <- prcomp((smallGeneExpress3[,1:50]), scale =T)
#look at variance
summary(pca_res)$importance
#var_explained <- summary(pca_res)$importance[2,]*100
#get 1&2nd PCs to plot
pcs <- data.frame(pca_res$x)
pcs$significant <- smallGeneExpress3$significant


#t_subset <- smallGeneExpress[,50:51]
#t.test(culture_fibroblast ~ significant, data = smallGeneExpress)

tmp <- t.test(PC1 ~ significant, data = pcs)
tmp$p.value
tmp$statistic
#empty vector to store pvalue results:
t_testPvaluesPC <- numeric()
#so I am going to loop through all the principal components: minus 1 b/c don't want to include last column
for (pc in 1:(ncol(pcs)-1)) {
  #print(pc)
  #extract the result of t-test for each
  tresult[[pc]] <- t.test(pcs[[pc]] ~ significant, data = pcs)
  #obtain the p-value for each
  pVal[[pc]] <- tresult[[pc]]$p.value
  #then add it to empty numeric vector created above
  t_testPvaluesPC <- c(t_testPvaluesPC, pVal[[pc]])
}
head(data.frame(tres))






