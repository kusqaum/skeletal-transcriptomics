#integrate network data with gene expression data


#embeddings created using pecanpy in python
embeddings <- read.table("processed/networkEdgeList.emb", sep = "\t", skip = 1)
colnames(embeddings)

#read in gene expression data:
fullGeforPCA <- readRDS("processed/fullGeneExpressForPCA.rds")
head(fullGeforPCA[1:5,1:3])
#and labels:
labelsFullDf <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
head(labelsFullDf)
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
embTest <- embeddings %>%
  separate(colnames(embeddings), c("hgnc_symbol", paste0("D", 1:500)), sep = " ")

# embTest <- mutate_all(embTest, function(x) as.numeric(as.character(x)))
# rownames(emb) <- emb$D0


emb <- readNetworkData("processed/networkEdgeList.emb", dimens = 500)
origEmb <-readNetworkData("processed/origNetworkEdgeList.emb", 1000)
dim(origEmb)
head(origEmb[1:4,1:6])
head(origEmb[,999:1001])
#now need to convert the genes to ensembl IDs 
netGenes <- data.frame(hgnc_symbol = origEmb$hgnc_symbol)
dim(netGenes)
networkGenes <- data.frame(hgnc_symbol = emb$hgnc_symbol)
write.table(netGenes, "processed/networkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)
write.table(networkGenes,"processed/AllnetworkGenesExportForMapping.txt", sep = "\t", row.names = F, quote = F)

##now read in the mapped file 
# networkGenesWithEnsembl <- read.table("raw/networkGenesWithEnsembl.txt", sep = "\t", header = T)
networkGenesWithEnsembl <- read.table("raw/networkGenesWithEnsembl.txt", sep = "\t", header = T)
# networkGenesWithEnsembl <- networkGenesWithEnsembl %>% 
#   distinct(hgnc_symbol, .keep_all = T)
#embTest$hgnc_symbol <- embTest$D0
#rownames(networkGenesWithEnsembl)<- networkGenesWithEnsembl$hgnc_symbol
# netData_2 <- merge(networkGenesWithEnsembl, embTest, by="hgnc_symbol")
networkData <- merge(networkGenesWithEnsembl, origEmb, by="hgnc_symbol")
head(networkData[1:4,1:3])
rownames(networkData)<- networkData$ensembl_gene_id; networkData$hgnc_symbol<- NULL; networkData$ensembl_gene_id<-NULL
str(networkData$D2)
dim(networkData)
head(networkData[,996:1000])

#now need to merge with gene exp
networkData <- mutate_all(networkData, function(x) as.numeric(as.character(x)))

fullGPCA <- merge(fullGeforPCA, networkData, by = 0)
head(fullGPCA[1:3,1:3])
dim(fullGPCA)
head(fullGPCA[15342:15346,3310:3312])
rownames(fullGPCA)<- fullGPCA$Row.names; fullGPCA$Row.names <- NULL
plan(multisession, workers=availableCores())
dim(fullGPCA)
head(fullGPCA[15339:15346, 3309:3312])
integratedDataPcaRes <- prcomp(fullGPCA, scale=T)
integPCs <- as.data.frame(integratedDataPcaRes$x)
#newpcares <- prcomp(fullGPCA, scale. = T)

integ_200dim <- as.data.frame(integPCs)[,1:200]
integ_50dim <- as.data.frame(integPCs)[,1:50]

integlabelled_50dim <- merge(integ_50dim, labelsFullDf, by=0); integlabelled_50dim$Row.names<-NULL
integlabelled_200dim <- merge(integ_200dim, labelsFullDf, by=0); integ_100dim$Row.names<-NULL

#split data and create rec
integsplit200 <- split_processingData_fctn(integlabelled_200dim, 0.8)
integsplit50<- split_processingData_fctn(integlabelled_50dim, 0.8)

#ml
plan(multisession, workers=availableCores())
plan(multisession, workers=availableCores())
plan(multisession, workers = availableCores())
mlresultWithPPI_200dim <- justToTestRF(integsplit200, algorithm = "PCA")
mlresultWithPPI_50dim <- justToTestRF(integsplit50, algorithm = "PCA")

##




#here are the functions used:
split_processingData_fctn <- function(data, proportion){
  set.seed(123)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  rec <-recipe(significant~., data = trainData) %>%
    step_downsample(significant, under_ratio = 1, seed = 456)
  #return(list(rec, trainData, testData))
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}


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
#
