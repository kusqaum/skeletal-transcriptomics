#integrate network data with gene expression data
library(tidyverse)
#might move this to a different script
library(tidymodels)
library(themis)
library(vip)
library(doParallel)
library(foreach)
library(NMF)

#uncommcent from line 12 to 93
# #embeddings created using pecanpy in python
# embeddings <- read.table("processed/networkEdgeList.emb", sep = "\t", skip = 1)
# colnames(embeddings)
# 
# #read in gene expression data:
# # fullGeforPCA <- readRDS("processed/fullGeneExpressForPCA.rds")
# # head(fullGeforPCA[1:5,1:3])
# # #and labels:
# # labelsFullDf <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
# # head(labelsFullDf)
# #importing node2vec result now
# #split into 1000 dims because that is the dimensions used for pecanPy
# readNetworkData <- function(embData, dimens){
#   #dimens is the number of dimensions used for pecanpy
#   embedding <- read.table(embData, sep = "\t", skip = 1)
#   emb <- embedding %>%
#     separate(colnames(embedding), c("hgnc_symbol", paste0("D", 1:dimens)), sep = " ")
#   #rownames(emb) <- emb$D0; emb$D0 <- NULL
#   return(emb)
# }
# 
# #just separate because all has been put into one column
# # embTest <- embeddings %>%
# #   separate(colnames(embeddings), c("hgnc_symbol", paste0("D", 1:500)), sep = " ")
# 
# # embTest <- mutate_all(embTest, function(x) as.numeric(as.character(x)))
# # rownames(emb) <- emb$D0
# 
# 
# emb <- readNetworkData("processed/networkEdgeList.emb", dimens = 500)
# #origEmb <-readNetworkData("processed/origNetworkEdgeList.emb", 1000)
# #dim(origEmb)
# dim(emb)
# #head(origEmb[1:4,1:6])
# #head(origEmb[,999:1001])
# #now need to convert the genes to ensembl IDs
# netGenes <- data.frame(hgnc_symbol = emb$hgnc_symbol)
# dim(netGenes)
# networkGenes <- data.frame(hgnc_symbol = emb$hgnc_symbol)
# write.table(netGenes, "processed/networkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)
# #write.table(networkGenes,"processed/AllnetworkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)
# 
# ##now read in the mapped file
# 
# networkGenesWithEnsembl <- read.table("raw/allNetworkGenesMapped.txt", sep = "\t", header = T)
# # networkGenesWithEnsembl <- networkGenesWithEnsembl %>%
# #   distinct(hgnc_symbol, .keep_all = T)
# #embTest$hgnc_symbol <- embTest$D0
# #rownames(networkGenesWithEnsembl)<- networkGenesWithEnsembl$hgnc_symbol
# # netData_2 <- merge(networkGenesWithEnsembl, embTest, by="hgnc_symbol")
# networkData <- merge(networkGenesWithEnsembl, emb, by="hgnc_symbol")
# head(networkData[1:4,1:3])
# rownames(networkData)<- networkData$ensembl_gene_id; networkData$hgnc_symbol<- NULL; networkData$ensembl_gene_id<-NULL
# str(networkData$D2)
# dim(networkData)
# head(networkData[,488:500])
# 
# #now need to merge with gene exp
# networkData <- mutate_all(networkData, function(x) as.numeric(as.character(x)))
# 
# sVGeneExpress <- readRDS("processed/fullGeneExpressForNMF.rds")
# geneExpressWithNet <- merge(sVGeneExpress, networkData, by = 0)
# head(geneExpressWithNet[1:3,1:3])
# 
# dim(geneExpressWithNet)
# 
# rownames(geneExpressWithNet)<- geneExpressWithNet$Row.names; geneExpressWithNet$Row.names <- NULL
# minVal <- abs(min(geneExpressWithNet))
# minVal
# geneExpressWithNetoffSet <- geneExpressWithNet+minVal
# ranks = c(5,10,50,100,150,200,500)
# noofruns <- 2
# set.seed(1234)
# 
# # cl <- makeCluster(25, outfile="processed/nmfSvWithNetwork.txt")
# # registerDoParallel(cl)
# # # integratedDataPcaRes <- prcomp(fullGPCA, scale=T)
# # print("running NMF")
# # resInteg.multiRank <- nmf(geneExpressWithNetoffSet, rank = ranks, nrun=noofruns, seed = 123456)
# # #
# 
# # saveRDS(resInteg.multiRank,"processed/resInteg.multiRank.rds")
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
labelsForNet <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
labelsForNet$significant <- as.factor(labelsForNet$significant)
#might put this in separate script?

sVwithNetworkLabelled <- lapply(matricesWithNet, getLabelledGenesFctn, labelsForNet)
head(sVwithNetworkLabelled[[1]][1:5,1:5])
dim(sVwithNetworkLabelled[[1]])
#function for preprocessing
split_processingData_fctn <- function(data, proportion){
  set.seed(123)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  rec <-recipe(significant~., data = trainData) %>%
    #step_adasyn(over_ratio = 1, seed = 456) %>%
    step_downsample(significant, under_ratio = 1, seed = 456)
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}

sVwithNetworkSplit <- lapply(sVwithNetworkLabelled, split_processingData_fctn, 0.8)

randForest_fctnWithNetwork <- function(preprocessResult, algorithm){
  print(paste0("Number of predictors = ", ncol(preprocessResult$train)-1))
  if(algorithm == "NMF"){
    
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", 
                                                                    importance = TRUE)
    
    set.seed(234)
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 3)
    levels <- NULL
    if(ncol(preprocessResult$train)-1 < 50){
      levels <- ncol(preprocessResult$train)-1
    }
    else if (ncol(preprocessResult$train)-1 >=50){
      levels <- 18
    }
    print(paste0("Number of levels = ", levels))
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      trees(range = c(1,2000)),
      mtry(range = c(1,limit)),
      min_n(range = c(1,limit)),
      levels = levels
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
      control = control_grid(save_pred = TRUE, verbose = TRUE),
      metrics = metric_set(roc_auc),
    )
    # rfbayes <- tune_bayes(wkflow, folds,
    #                       initial = res,
    #                       control = control_bayes(verbose = TRUE, save_pred=TRUE))
    
    #make a df for correct R
    df <- data.frame(k = res$id2, r = res$id)
    metrics <- res$.metrics
    getBestMet <- function(metr){
      bestRows <- metr[which.max(metr$.estimate),]
      return(bestRows)
    }
    bestMet <- map(.x = metrics, .f = getBestMet)
    allbestMet <- do.call("rbind", bestMet)
    df$values <- allbestMet$.estimate
    df <- df%>% mutate(model = "RF")
    df$k <- as.integer(substr(df$k, nchar(df$k), nchar(df$k)))
    df$r <- as.integer(substr(df$r, nchar(df$r), nchar(df$r)))
    
    dof <- length(res$splits)-1
    alpha <- 0.05
    tScore <- qt(p = alpha, df = dof, lower.tail = F)
    resDf <- res %>% collect_metrics() %>%
      mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm, marginError = std_err* tScore) %>%
      mutate(lowerbound = mean - marginError, upperbound = mean + marginError)
    
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
    
    roc_auc <- roc_auc(aug, significant, .pred_FALSE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_FALSE)
    rocCurve <- autoplot(two_classCurve)
    dimension <- ncol(preprocessResult$train)-1
    
    result <- list("workflow" = wkflow,"res"= res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                   "dim" = dimension, "dfForCorrectR" = df)
    saveRDS(result, sprintf("processed/%sRFResNetworkWithSV_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, 
                "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                "dim" = dimension, "dfForCorrectR" = df))
    
  }
  else if(algorithm == "PCA"){
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    print(ncol(preprocessResult$train))
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 3)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      trees(range = c(1,2000)),
      mtry(range = c(1,limit)),
      min_n(range = c(1,limit)),
      levels = limit
    )
    #build workflow
    wkflow <- workflow() %>%
      add_recipe(preprocessResult$recipe) %>%
      add_model(model_RF)
    #tune model
    res <- tune_grid(
      wkflow,
      resamples = folds,
      grid = tuningGrid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc),
    )
    
    
    
    df <- data.frame(k = res$id2, r = res$id)
    metrics <- res$.metrics
    getBestMet <- function(metr){
      bestRows <- metr[which.max(metr$.estimate),]
      return(bestRows)
    }
    bestMet <- map(.x = metrics, .f = getBestMet)
    allbestMet <- do.call("rbind", bestMet)
    df$values <- allbestMet$.estimate
    df <- df%>% mutate(model = "RF")
    df$k <- as.integer(substr(df$k, nchar(df$k), nchar(df$k)))
    df$r <- as.integer(substr(df$r, nchar(df$r), nchar(df$r)))
    
    rf_metrics <- res%>% collect_metrics()
    dof <- length(res$splits)-1
    alpha <- 0.05
    tScore <- qt(p = alpha, df = dof, lower.tail = F)
    resDf <- res %>% collect_metrics() %>%
      mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm, marginError = std_err* tScore) %>%
      mutate(lowerbound = mean - marginError, upperbound = mean + marginError)
    paramsPlot <- autoplot(res)
    paramsPlot2 <- res %>%
      collect_metrics() %>%
      ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
      geom_point()+
      geom_line()
    
    final_model <- res %>% select_best(metric = "roc_auc")
    
    final_fit <- finalize_workflow(wkflow, final_model) %>%
      fit(data = preprocessResult$train) 
    
    importancePlot <- final_fit %>% extract_fit_parsnip()%>%
      vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
      theme_classic()
    
    importanceDf <- final_fit %>% extract_fit_parsnip() %>%
      vi() %>% as.data.frame()
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    
    roc_auc <- roc_auc(aug, significant, .pred_FALSE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_FALSE)
    rocCurve <- autoplot(two_classCurve)
    dimension <- ncol(preprocessResult$train)-1
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "dfForCorrectR" = df,"resDf" = resDf, "finalMod" = final_model, 
                   "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                   "dim" = dimension)
    saveRDS(result, sprintf("processed/%sMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "dfForCorrectR" = df,"resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                "dim" = dimension))
  }
  
}

cl <- makeCluster(25, outfile="processed/200dimNMF_NetworkRF.txt")
registerDoParallel(cl)
print("randforest on 200 dim NMF SkeletalVis integrated with network:")
nmfNetworkdim200Res <- randForest_fctnWithNetwork(sVwithNetworkSplit[[6]], algorithm = "NMF")




