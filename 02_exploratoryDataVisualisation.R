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

lowestValue <- abs(min(fullGeneExpressforNMF)) # there are no negative values now

#these are what will be used for NMF and PCA
fullGeneExpressforNMF[1,1]+lowestValue
fullGeneExpressforNMF <- fullGeneExpressforNMF+lowestValue
head(fullGeneExpressforNMF[1:5,1:5])
fullGeneExpressforPCA <- as.data.frame(fullGeneExpressN)
head(fullGeneExpressforPCA[1:5,1:5])
saveRDS(fullGeneExpressforPCA, "processed/fullGeneExpressForPCA.rds")

# sum(smallGeneExpress2$significant == "TRUE")
# sum(smallGeneExpress2$significant == "FALSE")
#now we have 1300 and 6578 positive and negative genes


#nmfSeed()
####perform NMF####
#determine the sequence
ranks = c(5,10,50,100,150,200)
#noRanks <- seq(2,10,2)
#and the number of ranks
noofruns <- 2
#set seed as well
set.seed(1234)
Ranks <- c(5,10)

##generate shuffled data
shuffledNMF <- randomize(fullGeneExpressforNMF); row.names(shuffledNMF)<- row.names(fullGeneExpressforNMF)
# library(future)
# plan(multisession, workers = availableCores())

res.multiRank <- nmf(fullGeneExpressforNMF, rank = ranks, nrun=noofruns, seed = 123456)
saveRDS(res.multiRank, "processed/res.multiRank.rds")
# res.multiRank <- nmf(as.data.frame(shuffled)[1:10], rank = c(5,10), nrun=noofruns, seed = 123456)
#^ gives the same result for t test

# res.test <- nmf(fullGeneExpressforNMF[1:10], rank = c(5,10), seed=123456)

#nmf without logg transf
#res.multiRank2 <- nmf(smallGeneExpress[,1:50], rank = noRanks, nrun=noofruns, seed = 123456)

#look at performance measures of the factorisation
summary(res.multiRank)
consensusmap(res.multiRank, labCol = NA, labRow = 1)

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

#here, from NMF reduction result I will get up to only the genes that have labels from each matrix:
# so that can be used for machine learning:

getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}
#this is for NMF machine learning -cool
labelledGenesNMFRes <- lapply(Wmatrices, FUN=getLabelledGenesFctn, labelsFullDf)
#dim to check merged correctly just in case
dim(labelledGenesNMFRes[[1]])
head(labelledGenesNMFRes[[1]])
labelledGenesNMFResTest <- labelledGenesNMFRes
s<- lapply(Wmatrices, FUN=getLabelledGenesFctn, labelsFullDf)


####perform PCA####

#perform PCA. each of our data points needs to be a gene so not going to transpose
pca_res <- prcomp(fullGeneExpressforPCA, scale = T)
pcs <- data.frame(pca_res$x)
dim(pcs)
pcsLabelled <- merge(pcs, labelsFullDf, by=0)
dim(pcsLabelled)
row.names(pcsLabelled) <- pcsLabelled$Row.names; pcsLabelled$Row.names <- NULL
dim(pcsLabelled)

tmp <- t.test(y=as.logical(pcsLabelled$significant), x=pcsLabelled[,518])
tmp$p.value
tmp$statistic
t.test(pcsLabelled$PC2 ~ mouseHumanWithLabs$significant)$p.value

####do pca on shuffled data and then let's do machine learning
shuffledPCA<- randomize(fullGeneExpressforPCA); rownames(shuffledPCA)<- rownames(fullGeneExpressforPCA)
pcaRes_shuffled <- prcomp(shuffledPCA, scale. = T)
shuffPCs <- pcaRes_shuffled$x
#just test on dim 50 or 100
shuffPCsdim50 <- shuffPCs[,1:50]
#add on labels:
shuffPCsdim50_L <- merge(shuffPCsdim50, labelsFullDf, by=0)
rownames(shuffPCsdim50_L)<- shuffPCsdim50_L$Row.names ;shuffPCsdim50_L$Row.names<- NULL
head(shuffPCsdim50_L)

#split and create rec
shuff_Split <- split_processingData_fctn(shuffPCsdim50_L, proportion = 0.8)
shuff_MLRes <- justToTestRF(shuff_Split, algorithm = "PCA")

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
unlabelledNMFgenes <- lapply((labelledGenesNMFRes), removeLabelsFctn)
pcaDimPValues <- lapply(X = pcList, FUN = performTtestFctn, labels = labelsFullDf$significant, algorithm="PCA")
 # pcaDimPValues2 <- lapply(X=pcList2, FUN = performTtestFctn, labels = mouseHumanWithLabs$significant, algorithm="PCA")
nmfMinPvals <- lapply(X=unlabelledNMFgenes, FUN = performTtestFctn, labels = mouseHumanWithLabs$significant, algorithm="NMF")
t.test(labelledGenesNMFRes[[5]][,4] ~ labelledGenesNMFRes[[5]][,ncol(labelledGenesNMFRes[[5]])])$p.value
ncol(labelledGenesNMFRes[[1]])
head(pcList[[1]])
pl <- pcList[[1]][,2]
then <- data.frame(merge(pl, labelsFullDf$significant, by=0))


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
  xlab("k dimensions") +ylab("-log10 P-Value")+
  scale_colour_manual(values = c("mediumpurple", "tan2"))+
  scale_size_manual(values = c(7, 5))
  # scale_x_continuous(breaks = seq(2,50, 2))

#now save this plot



labelledGenesPCARes <- lapply(pcList, FUN = getLabelledGenesFctn, labelsFullDf)
# labelledGenesNMFRes

saveRDS(labelledGenesPCARes, "processed/pcaDataframes.rds")
saveRDS(labelledGenesNMFRes, "processed/nmfDataframes.rds")

##save result for model training
# saveRDS(basis_matrices, "processed/nmf_wMatrices.rds")
# saveRDS(smallGeneExpress2, "processed/geneExpression.rds")
# saveRDS(smallGeneExpress2, "processed/tempGeneExpress.rds")



####i want to get all the unlabelled genes####
fullGeneExpressN <- as.data.frame(fullGeneExpressN)
noLabels <- row.names(fullGeneExpressN)[!row.names(fullGeneExpressN) %in%
                                          row.names(mouseHumanWithLabs)]
unlabelled <- fullGeneExpressN[noLabels,]
saveRDS(unlabelled, "processed/unlabelledGenes.rds")




#--------------------------------------------------------------------------
#temp for go terms potentially
# just extract pc50 for now:
pc50 <- pcList[[3]]
pc50_1<-pc50
pc50_1$significant <- labelsFullDf$significant
goTerms <- goTerms %>% distinct(ensembl_gene_id, .keep_all = T)
rownames(goTerms)<- goTerms$ensembl_gene_id
pc50go <- merge(pc50, goTerms, by=0)
rownames(pc50go) <- pc50go$Row.names; pc50go$Row.names <- NULL; pc50go$ensembl_gene_id<- NULL

pc50go$significant <- labelsFullDf$significant
pc50go2 <- pc50go%>% filter(go_id!="")

GOsplit_processingData_fctn <- function(data, proportion){
  set.seed(123)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  rec <-recipe(significant~., data = trainData) %>%
    step_downsample(significant, under_ratio = 1, seed = 456) %>%
    step_dummy(go_id)
  #return(list(rec, trainData, testData))
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}

split <- GOsplit_processingData_fctn(pc50go2, .8)
split2<- split_processingData_fctn(pc50_1, 0.8)
plan(multisession, workers=availableCores())
resultMLWithGO<- justToTestRF(split, algorithm="PCA")

###
#function:
justToTestRF <- function(preprocessResult, algorithm){
  if(algorithm == "NMF"){
    model_RF <- rand_forest(trees = 500, mtry = sqrt(ncol(preprocessResult$train)), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=3)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      #trees(range = c(1,2000)),
      #mtry(range = c(1,limit)),
      min_n(range = c(1,limit)),
      levels = limit
    )
    #build workflow
    wkflow <- workflow() %>%
      add_recipe(preprocessResult$recipe) %>%
      add_model(model_RF)
    #tune model
    #plan(multisession, workers = 16)
    res <- tune_grid(
      wkflow,
      resamples = folds,
      grid = tuningGrid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc),
    )
    resDf <- res %>% collect_metrics() %>%
      mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    paramsPlot <- autoplot(res)
    paramsPlot2 <- res %>%
      collect_metrics() %>%
      ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
      geom_point()+
      geom_line()
    
    final_model <- res %>% select_best(metric = "roc_auc")
    
    final_fit <- finalize_workflow(wkflow, final_model) %>%
      parsnip::fit(data = preprocessResult$train) 
    
    importancePlot <- final_fit %>% extract_fit_parsnip()%>%
      vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
      theme_classic()
    
    importanceDf <- final_fit %>% extract_fit_parsnip() %>%
      vi() %>% as.data.frame()
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    
    roc_auc <- roc_auc(aug, significant, .pred_TRUE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_TRUE)
    rocCurve <- autoplot(two_classCurve)
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    # saveRDS(result, sprintf("%sGTEXtestingMLRes_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
    
  }
  else if(algorithm == "PCA"){
    model_RF <- rand_forest(trees = 500, mtry = sqrt(ncol(preprocessResult$train)), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=3)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      #trees(range = c(1,2000)),
      #mtry(range = c(1,limit)),
      min_n(range = c(1,limit)),
      levels = limit
    )
    #build workflow
    wkflow <- workflows::workflow() %>%
      add_recipe(preprocessResult$recipe) %>%
      add_model(model_RF)
    #tune model
    #plan(multisession, workers = 16)
    res <- tune_grid(
      wkflow,
      resamples = folds,
      grid = tuningGrid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc),
    )
    resDf <- res %>% collect_metrics() %>%
      mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    paramsPlot <- autoplot(res)
    paramsPlot2 <- res %>%
      collect_metrics() %>%
      ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
      geom_point()+
      geom_line()
    
    final_model <- res %>% select_best(metric = "roc_auc")
    
    final_fit <- finalize_workflow(wkflow, final_model) %>%
      parsnip::fit(data = preprocessResult$train) 
    
    importancePlot <- final_fit %>% extract_fit_parsnip()%>%
      vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
      theme_classic()
    
    importanceDf <- final_fit %>% extract_fit_parsnip() %>%
      vi() %>% as.data.frame()
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    
    roc_auc <- roc_auc(aug, significant, .pred_TRUE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_TRUE)
    rocCurve <- autoplot(two_classCurve)
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    #saveRDS(result, sprintf("processed/%sGTEXtestingMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "resDF"=resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
  }
  
}
