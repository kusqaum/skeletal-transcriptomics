library(tidyverse)
library(tidymodels)
library(themis)
library(skimr)
library(randomForest)
library(vip)
library(xgboost)
library(foreach)
library(doParallel)
# library(preprocessCore)
# library(future)
# library(future)

getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}
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

geneExpress <- readRDS("processed/unlabelledNMFDfs.rds")
head(geneExpress[[1]][1:5,1:3])
mortalityLabels <- read.table("processed/processedMortalityLabels.txt", header = T, sep = "\t")
mortalityLabels$significant <- as.factor(mortalityLabels$significant)
predictMortGenes <- lapply(geneExpress, getLabelledGenesFctn, mortalityLabels)

head(predictMortGenes[[5]][1:4,5:7])
dim(predictMortGenes[[5]])

mortSplit <- lapply(predictMortGenes, split_processingData_fctn , .8)

randomForestFctn_positiveControl <- function(preprocessResult, algorithm){
  print(ncol(preprocessResult$train))
  if(algorithm == "NMF"){
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    print(paste0("number of predictors", ncol(preprocessResult$train)-1))
    set.seed(234)
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 3)
    levels <- NULL
    if(ncol(preprocessResult$train)-1 <= 50){
      levels <- 5
    }
    else if (ncol(preprocessResult$train)-1 >50){
      levels <- 15
    }
    limit <- (ncol(preprocessResult$train))-1
    print(paste0("Levels = ", levels))
    tuningGrid <- grid_regular(
      trees(range = c(1, limit)),
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
    # df <- data.frame(k = res$id2, r = res$id)
    # metrics <- res$.metrics
    # getBestMet <- function(metr){
    #   bestRows <- metr[which.max(metr$.estimate),]
    #   return(bestRows)
    # }
    # bestMet <- map(.x = metrics, .f = getBestMet)
    # allbestMet <- do.call("rbind", bestMet)
    # df$values <- allbestMet$.estimate
    # df <- df%>% mutate(model = "RF")
    # df$k <- as.integer(substr(df$k, nchar(df$k), nchar(df$k)))
    # df$r <- as.integer(substr(df$r, nchar(df$r), nchar(df$r)))
    # 
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
    
    result <- list("workflow" = wkflow,"res"= res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                   "dim" = dimension)
    saveRDS(result, sprintf("processed/%sRFResSVmortality_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, 
                "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                "dim" = dimension))
    
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
    #saveRDS(result, sprintf("processed/%sMLResSVmortality_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "dfForCorrectR" = df,"resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                "dim" = dimension))
  }
  
}

clus <- makeCluster(25, outfile = "processed/50dimMortalityRFSV.txt")
registerDoParallel(clus)
print("Run RF 50 dim to predict mortality")
# tempRes <- randomForestFctn_positiveControl(mortSplit[[3]], algorithm = "NMF")
# tempRes2 <- randomForestFctn_positiveControl(mortSplit[[2]], algorithm = "NMF")
res50dim <- randomForestFctn_positiveControl(mortSplit[[3]], algorithm = "NMF")

print("finished")


##
#now read in the ML results:

mortFiles <- list.files("processed", pattern = "ResSVmortality", full.names = T)
mortalityRes <- lapply(mortFiles, readRDS)
mortalityRes <- mortalityRes[order(sapply(mortalityRes, function(x) x$dim))]

mortPred <- lapply(mortalityRes, function(x){
  x$aug
})
mortPredDf <- bind_rows(mortPred)
mortalityROC <- mortPredDf %>% group_by(dim)%>%
  roc_curve(truth = significant, .pred_FALSE)%>%
  ggplot( aes(x = 1-specificity, y = sensitivity, colour = as.factor(dim)))+
  geom_path(linewidth = 0.7)+
  geom_abline(slope = 1, intercept = 0, size = 0.4, lty = "dashed")+
  theme(panel.border = element_rect(colour = "black", linewidth = 1, fill = "white"))+
  theme_bw(base_size = 28)+
  scale_colour_npg()+
  facet_wrap(~ factor(paste0(dim, " NMF dimensions"), c("50 NMF dimensions", "100 NMF dimensions",
                                                       "150 NMF dimensions","200 NMF dimensions"))) +
  theme(legend.position = "none")+
  theme(strip.text = element_text(size =20))

mortalityAUC <- lapply(mortalityRes, function(x){
  x$AUC
})
mortalityAUC <- do.call("rbind",mortalityAUC)
  
ggsave("output/mortalityCurves.png", mortalityROC, width = 14.5, height = 15)
# gtex <- read.delim("raw/rnaSeqGtex/gtex_Analysis.gct", header = T,skip = 2, sep = "\t")
# 
# 
# 
# #just retain protein coding genes
# # temp <- gtex[1:5,1:3]
# # 
# # split the rownames
# geneIDsgtex <- as.data.frame(do.call('rbind', strsplit(as.character(gtex$Name),
#                                    '.', fixed = T)))
# 
# #remove first 2 columns now
# #gtex <- gtex[,-c(1:2)]
# 
# length(unique(geneIDsgtex$V1))
# 
# humanProteinCoding <- read.table("processed/human_coding_genes.txt", header = T,
#                                  stringsAsFactors = F, sep = "\t")
# #
# 
# row.names(gtex) <- geneIDsgtex$V1
# head(gtex[1:3,1:3])
# #also get rid of name and description columns:
# gtex <- gtex[-c(1,2)]
# 
# gtex$ensembl_gene_id <- rownames(gtex)
# gtexProteinCoding <- merge(humanProteinCoding, gtex, by = "ensembl_gene_id")
# row.names(gtexProteinCoding)<- gtexProteinCoding$ensembl_gene_id; gtexProteinCoding$ensembl_gene_id <- NULL
# 
# 
# #now we have all the gene expression data in ensembl ID format
# filtgtex <- gtexProteinCoding[rowSums(gtexProteinCoding)>0,]
# rowsFiltgtex <- rownames(filtgtex)
# head(filtgtex[1:5,1:5])
# # gtexNor <- normalize.quantiles(as.matrix(filtgtex))
# # rownames(gtexNor)<- rowsFiltgtex
# # gtexNor<- as.data.frame(gtexNor)
# # 
# # 
# # #perform PCA
# # pca_gtex <- prcomp(gtexProteinCoding, scale. = T)
# # pcaGtexpcs <- pca_gtex$x
# # for (i in 1:length(ranks)) {
# #   #create empty lists
# #   gtexpcList <- list()
# #   #loop through each of the dimensions
# #   for (l in (ranks)) {
# #     #grab each dimension
# #     dime <- pcsLabelled[,1:l]
# #     gtexpcList <- append(gtexpcList, list(dime))
# #   }
# # }
# #add labels function
# getLabelledGenesFctn <- function(matrixList, knownLabels){
#   merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
#   return(merged)
# }
# # 
# # #apply to result
# # labelledgtex <- map(.x = gtexpcList, .f = getLabelledGenesFctn, knownLabels = labelsFullDf)
# 
# # gtextTest <- map(.x = gtexpcList,.f = performTtestFctn, labels=labelsFullDf$significant, algorithm = "PCA")
# # 
# # gtex_5dim <- labelledgtex[[5]]
# #temporary RF function not using too much tuning or repeats of cross-validation just for testing
# gtexRFfctn <- function(preprocessResult, algorithm){
#   if(algorithm == "NMF"){
#     model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
#                             mode = "classification") %>% set_engine("randomForest", importance = TRUE)
#     set.seed(234)
#     
#     folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 5)
#     limit <- (ncol(preprocessResult$train))-1
#     tuningGrid <- grid_regular(
#       trees(range = c(1,2000)),
#       mtry(range = c(1,limit)),
#       min_n(range = c(1,limit)),
#       levels = limit
#     )
#     #build workflow
#     wkflow <- workflow() %>%
#       add_recipe(preprocessResult$recipe) %>%
#       add_model(model_RF)
#     #tune model
#     #plan(multisession, workers = 16)
#     res <- tune_grid(
#       wkflow,
#       resamples = folds,
#       grid = tuningGrid,
#       control = control_grid(save_pred = TRUE),
#       metrics = metric_set(roc_auc),
#     )
#     resDf <- res %>% collect_metrics() %>%
#       mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
#     paramsPlot <- autoplot(res)
#     paramsPlot2 <- res %>%
#       collect_metrics() %>%
#       ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
#       geom_point()+
#       geom_line()
#     
#     final_model <- res %>% select_best(metric = "roc_auc")
#     
#     final_fit <- finalize_workflow(wkflow, final_model) %>%
#       parsnip::fit(data = preprocessResult$train) 
#     
#     importancePlot <- final_fit %>% extract_fit_parsnip()%>%
#       vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
#       theme_classic()
#     
#     importanceDf <- final_fit %>% extract_fit_parsnip() %>%
#       vi() %>% as.data.frame()
#     
#     aug <- augment(final_fit, preprocessResult$test)
#     aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
#     
#     roc_auc <- roc_auc(aug, significant, .pred_1)
#     two_classCurve <- roc_curve(aug, truth = significant,
#                                 .pred_1)
#     rocCurve <- autoplot(two_classCurve)
#     
#     #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
#     result <- list("workflow" = wkflow, "res" = res, "resDf"= resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
#                    "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
#     # saveRDS(result, sprintf("%sGTEXMLRes_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
#     
#     return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
#                 "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
#     
#   }
#   else if(algorithm == "PCA"){
#     model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
#                             mode = "classification") %>% set_engine("randomForest", importance = TRUE)
#     set.seed(234)
#     
#     folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 5)
#     limit <- (ncol(preprocessResult$train))-1
#     tuningGrid <- grid_regular(
#       trees(range = c(1,2000)),
#       mtry(range = c(1,limit)),
#       min_n(range = c(1,limit)),
#       levels = limit
#     )
#     #build workflow
#     wkflow <- workflows::workflow() %>%
#       add_recipe(preprocessResult$recipe) %>%
#       add_model(model_RF)
#     #tune model
#     #plan(multisession, workers = 16)
#     res <- tune_grid(
#       wkflow,
#       resamples = folds,
#       grid = tuningGrid,
#       control = control_grid(save_pred = TRUE),
#       metrics = metric_set(roc_auc),
#     )
#     resDf <- res %>% collect_metrics() %>%
#       mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
#     paramsPlot <- autoplot(res)
#     paramsPlot2 <- res %>%
#       collect_metrics() %>%
#       ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
#       geom_point()+
#       geom_line()
#     
#     final_model <- res %>% select_best(metric = "roc_auc")
#     
#     final_fit <- finalize_workflow(wkflow, final_model) %>%
#       parsnip::fit(data = preprocessResult$train) 
#     
#     importancePlot <- final_fit %>% extract_fit_parsnip()%>%
#       vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
#       theme_classic()
#     
#     importanceDf <- final_fit %>% extract_fit_parsnip() %>%
#       vi() %>% as.data.frame()
#     
#     aug <- augment(final_fit, preprocessResult$test)
#     aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
#     
#     roc_auc <- roc_auc(aug, significant, .pred_1)
#     two_classCurve <- roc_curve(aug, truth = significant,
#                                 .pred_1)
#     rocCurve <- autoplot(two_classCurve)
#     
#     #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
#     result <- list("workflow" = wkflow, "res" = res, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
#                    "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
#     saveRDS(result, sprintf("processed/%sGTEXallMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
#     
#     return(list("workflow" = wkflow, "res" = res, "resDF"=resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
#                 "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
#   }
#   
# }
# 
# 
# ##----------------------------------------------------------------------------------------------------------------------------
# 
# sampleNames <- (colnames(gtexProteinCoding))
# names <- data.frame(id = sampleNames)
# # colnames(namesC)[1]<- id
# # df <- (do.call("rbind", strsplit(as.character(sampleNames$colnames(gtexProteinCoding), '.', fixed=T))))
# 
# # just get the first part of ensemblID
# id <- names %>% mutate(first = unlist(lapply(strsplit(names$id, '\\.'), function(x)x[1])))#
# id <- id %>% mutate(second = unlist(lapply(strsplit(names$id, '\\.',), function(x)x[2]))) 
# #make it as the same format as the sex is
# id$integrate <- paste(id$first, id$second, sep = '-')
# dim(id)
# colnames(gtexProteinCoding)
# #now make these the columnames
# colnames(filtgtex) <- id$integrate
# 
# gtexPhenotype <- read.table("raw/GTEX_v7_SubjectPhenotypeDS.txt", header = T) # 1 is male, 2 female
# length(unique(gtexPhenotype$SUBJID))
# # rownames(gtexPhenotype) <- gtexPhenotype$SUBJID
# # gtexPhenotypeNew <- data.frame(gtexPhenotype$SEX)
# # rownames(gtexPhenotypeNew)<- gtexPhenotype$SUBJID
# to_match <- colnames(filtgtex)
# 
# findMatch <- match(to_match,gtexPhenotype$SUBJID)
# sex <- gtexPhenotype[findMatch, 2]
# 
# #gtexTranspose <- as.data.frame(t(filtgtex))                       
# filtgtex_t <- as.data.frame(t(filtgtex))
# #filtgtex_t$sex <- sexdont do this yet because perform pca first!
# 
# filtgtex_tL <- log2(filtgtex_t+1)
# head(filtgtex[1:3,1:4])
# #normalise quantiles
# # filtgtexN <- normalize.quantiles(as.matrix(filtgtex_tL))
# # rownames(filtgtexN) <- rownames(filtgtex_tL)
# # filtgtexN <- as.data.frame(filtgtexN)
# # saveRDS(filtgtexN, "processed/filtgtexN.rds")
# filtgtexN <- readRDS("processed/filtgtexN.rds")
# ranks <- c(5,10,50,100,150,200)
# ## now let's perform PCA 
# # plan(multisession, workers=availableCores())
# print("performing pca")
# gtexPCres <- prcomp(filtgtexN, scale. = T)
# print("finished PCA")
# gtexPCs <- as.data.frame(gtexPCres$x)
# gtexPCs$sex <- sex
# for (g in 1:length(ranks)) {
#   #create empty lists
#   gtexpcList <- list()
#   #loop through each of the dimensions
#   for (l in (ranks)) {
#     #grab each dimension
#     dime <- gtexPCs[,1:l]
#     gtexpcList <- append(gtexpcList, list(dime))
#   }
# }
# #now i need to add on the labels for sex!!!
# 
# # rownames(sexDf) <- colnames(filtgtex)
# 
# gtexSex <- data.frame(significant = as.factor(gtexPCs[,ncol(gtexPCs)]))
# head(gtexSex)
# rownames(gtexSex)<- rownames(gtexPCs)
# 
# gtexPCsLabelled <- map(.x = gtexpcList , .f =  getLabelledGenesFctn, knownLabels = gtexSex)
# ("print split data")
# splitGTEXdata <- function(data, proportion){
#   set.seed(123)
#   dataSplit <- initial_split(data, prop = proportion)
#   trainData <- training(dataSplit)
#   testData <- testing(dataSplit)
#   #create recipe
#   rec <-recipe(significant~., data = trainData) %>%
#     step_downsample(significant, under_ratio = 1, seed = 456)
#   #return(list(rec, trainData, testData))
#   return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
# }
# print("create preprocessing recipe now")
# gtexProcessed <- lapply(X = gtexPCsLabelled, FUN = splitGTEXdata, 0.8)
# print("done creating recipe")
# print("run ML now")
# gtexAllMlRes <- lapply(X=gtexProcessed, FUN = gtexRFfctn, algorithm="PCA")
# print("finished ML")

