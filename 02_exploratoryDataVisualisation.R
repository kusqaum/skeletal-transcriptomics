#02_exlploratoryDataVisualisation
library(NMF)


geneExpressionWithLabels <- read.table("processed/geneExpressionDataWithLabels.txt", header = T, stringsAsFactors = F)
smallGeneExpress <- geneExpressionWithLabels %>% select(1:50)
smallGeneExpress$significant <- geneExpressionWithLabels$significant
#remove genes(rows) that have a sum of 0 abundance
smallGeneExpress <- smallGeneExpress %>% filter(rowSums(across(where(is.numeric)))!=0)
matrix <- as.matrix(smallGeneExpress)
expressionDataMat <- as.matrix(geneExpressionWithLabels)
nmfSeed()

noRanks <- seq(2,10,2)
noofruns <- 2
res.multiRank <- nmf(matrix[,1:50], rank = noRanks, nrun=noofruns, seed = 123456)
plot(res.multiRank)
#res.multiRank2 <- nmf(expressionDataMat[,1:1128], seq(2,200,5), seed = 123456)

consensusmap(res.multiRank, labCol = NA, labRow = NA)

#res.multi.method <- nmf(smallGeneExpress[,1:50], 2, seed =123456, list("brunet","lee", "ns"), .options= "t")
#compare(res.multi.method)

#maxbrunCoph <- max(res.multi.method$brunet$measures$cophenetic)
#maxNsCoph <- max(res.multi.method$nsNMF$measures$cophenetic)
compare(res.multi.method)
plot(res.multi.method)
plot(res.multiRank)
plot(res.multi.method$brunet$measures$cophenetic)
plot(res.multi.method$nsNMF$measures$cophenetic)

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
#let's look at the silhouette:
res.multiRank$measures$cophenetic
res.multiRank$measures$silhouette.consensus
nmfFitClass  <- res.multiRank$fit[[idxCophenMax]]
#get the w and h matrices for the rank with highest cophenetic
h_matrix <- nmfFitClass@fit@H
w_matrix <- nmfFitClass@fit@W
consensusmap(res.multiRank$fit[[idxCophenMax]])
consensusmap(res.multiRank)
res <- nmf(smallGeneExpress[,1:50], 2)
nmf.h <- basis(res)
#filename_nmf_pdf <- paste("NMF_rank_",idxCophenMax+1, nrow(matrix), "metagenes",ncol(matrix), "samples")
ggplot(as.data.frame(h_matrix))