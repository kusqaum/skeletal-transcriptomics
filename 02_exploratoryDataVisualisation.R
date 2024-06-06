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


###ignore for now
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


nmfPvals <- numeric()
minPvals <- numeric()
for (mat in basis_matrices) {
  #pValue <- rep
  #View(mat)
  for (dim in 1:ncol(mat)) {
    resTtest <- t.test(smallGeneExpress3$significant, mat[, dim])$p.value
    #mat$pValue <- resTtest
    nmfPvals <- c(nmfPvals, resTtest)
    #minPvals <- which.min(nmfPvals)
    #basis_matrices <- append(basis_matrices[[mat]], resTtest[[mat]])
    
  }
  #minPvals <- c(minPvals)
  
}





####PCA####

#perform PCA. each of our data points needs to be a gene so not going to transpose
pca_res <- prcomp((smallGeneExpress3[,1:50]), scale =T)
#look at variance
#summary(pca_res)$importance
#get 1&2nd PCs to plot
pcs <- data.frame(pca_res$x)
#pcs$significant <- smallGeneExpress3$significant


#t_subset <- smallGeneExpress[,50:51]
#t.test(culture_fibroblast ~ significant, data = smallGeneExpress)

tmp <- t.test(smallGeneExpress3$significant, pcs[,"PC1"])
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
  resultfortest <- t.test(smallGeneExpress3$significant,  pcs[,PC])$p.value
  pValsPCA <- c(pValsPCA, resultfortest)
  #resultfortest_pval <- resultfortest$
  smallestPval <- which.min(pValsPCA)
}
pValsdf <- as.data.frame(pValsPCA); pValsdf$dimension<- colnames(pcs)
pValsdf$dim <- c(1:ncol(pcs))
pValsdf$Algorithm <- c('PCA')
#t.test(smallGeneExpress3$significant, pcs$PC50)$p.value
ggplot(pValsdf, aes(x=dim, y=-log10(pValsPCA), col = Algorithm)) +
  geom_point()+
  geom_smooth(method = 'loess') +
  theme_bw(base_size = 18) +
  scale_x_continuous(breaks = seq(2,50, 2)) +
  xlab("k dimensions") + ylab("-log10 P-Value")

#also create a scree plot for pca
var_explained <- (pca_res$sdev)^2/sum((pca_res$sdev)^2)
varExpl_df <- data.frame(PC = 1:length(var_explained), 
                                var_expl = var_explained,
                         cumulativeVariance = cumsum(var_expl))

ggplot(varExpl_df[1:8,], aes(x = PC, y = 100*cumulativeVariance)) +
  geom_bar(stat = 'identity', col = 'black', fill="lightgrey") + 
  #geom_point()+
  #geom_line()+
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 95, col = "red", lty = "dashed")+
  scale_x_continuous(breaks = seq(1, nrow(varExpl_df[1:8,]), 1)) +
  xlab("PC index") + ylab("% explained variance by principal components")
  
  #xlim(0, 10)

####rough of machine learning:####
#library(tidymodels)


compressed_PCA <- pcs
compressed_PCA$significant <- smallGeneExpress3$significant
str(compressed_PCA)
compressed_PCA$significant <- as.factor(compressed_PCA$significant)
