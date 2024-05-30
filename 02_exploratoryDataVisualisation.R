#02_exlploratoryDataVisualisation

geneExpressionWithLabels <- read.table("processed/geneExpressionDataWithLabels.txt", header = T, stringsAsFactors = F)
smallGeneExpress <- geneExpressionWithLabels %>% select(1:50)
smallGeneExpress$significant <- geneExpressionWithLabels$significant
#remove genes(rows) that have a sum of 0 abundance
smallGeneExpress <- smallGeneExpress %>% filter(rowSums(across(where(is.numeric)))!=0)
matrix <- as.matrix(smallGeneExpress)
expressionDataMat <- as.matrix(geneExpressionWithLabels)

start <- 2
end <- 
seq(2,40,6)
res.multiRank <- nmf(matrix[,1:50], seq(2,40,8), nrun=2, seed = 123456)
plot(res.multiRank)
#res.multiRank2 <- nmf(expressionDataMat[,1:1128], seq(2,200,5), seed = 123456)

consensusmap(res.multiRank, labCol = NA, labRow = NA)
nmfSeed()

res.multi.method <- nmf(matrix, 2:8, nrun=5, seed =123456, list("brunet", "ns"), .options= "t")
compare(res.multi.method)
maxbrunCoph <- max(res.multi.method$brunet$measures$cophenetic)
maxNsCoph <- max(res.multi.method$nsNMF$measures$cophenetic)


plot(res.multi.method$brunet$measures$cophenetic)
plot(res.multi.method$nsNMF$measures$cophenetic)

res.multiRank$measures
# choosing the best rank
# so when thte cophenetic correlation coefficient begins to fall
max_cophenetic <- max(res.multiRank$measures$cophenetic)
# loop through the position of each rank.
for (i in 1:length(res.multiRank)){
  if (res.multiRank$measures$cophenetic[i]== max_cophenetic){
    idxCophenMax <- i #idxCophen is the position of the rank.
  }
}
idxCophenMax

res.multiRank$measures$silhouette.consensus
nmfFitClass  <- res.multiRank$fit[[idxCophenMax]]
h_matrix <- nmfFitClass@fit@H
w_matrix <- nmfFitClass@fit@W
consensusmap(res.multiRank$fit[[idx_cophenMax]])
