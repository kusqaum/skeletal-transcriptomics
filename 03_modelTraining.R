####rough of machine learning:####
library(tidymodels)
library(themis)
library(skimr)
library(tidyverse)
#library(parsnip) don't thin I actually need this cause gets loaded in with tidymodels
library(randomForest)
library(vip)
library(xgboost)
# library(future)

#will have to change these files to load in the latest files btw..
# nmfMatrices <- readRDS("processed/nmf_wMatrices.rds")
# geneExpress <- readRDS("processed/geneExpression.rds")
#these below are the full dataframes for each dimension 5to100
pcaDfs <- readRDS("processed/pcaDataframes.rds")
nmfDfs <- readRDS("processed/nmfDataframes.rds")
# str(pcaDfs[[1]]$significant)



####real thing####
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
  #return(list(rec, trainData, testData))
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}
print("creating preprocessing recipe")
#apply split data function to all matrices
# nmfProcessed <- lapply(X=nmfDfs, FUN = split_processingData_fctn, proportion =0.8)
nmfProcessed <- map(.x = nmfDfs, .f = split_processingData_fctn, proportion=0.8)
# pcaProcessed <- map(.x = pcaDfs, .f = split_processingData_fctn, proportion = 0.8)

randForest_fctn <- function(preprocessResult, algorithm){
  if(algorithm == "NMF"){
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 5)
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
    #plan(multisession, workers = 16)
    res <- tune_grid(
      wkflow,
      resamples = folds,
      grid = tuningGrid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc),
    )
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
    result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    saveRDS(result, sprintf("processed/%sMLRes_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                "dim" = dimension, "dfForCorrectR" = df))
    
  }
  else if(algorithm == "PCA"){
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 5)
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
    #plan(multisession, workers = 16)
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
      fit(data = preprocessResult$train) 
    
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
    dimension <- ncol(preprocessResult$train)-1
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    saveRDS(result, sprintf("processed/%sMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "dfForCorrectR" = df,"resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "algorithm" = algorithm,
                "dim" = dimension))
  }
  
}
# availableCores()
# plan(multisession, workers = availableCores())
print("performing randforest function on NMF datasets")
nmfDataMlResList <- lapply(nmfProcessed, FUN=randForest_fctn, algorithm="NMF")

# pcaDataMlResList <- lapply(pcaProcessed, FUN=randForest_fctn, algorithm="PCA")
# pcaNewList <- list()
# pcaNewList[[1]] <- pcaProcessed[[3]]
# pcaNewList[[2]] <- pcaProcessed[[4]]
# pcaNewList[[3]] <- pcaProcessed[[5]]
# pcaNewList[[4]] <- pcaProcessed[[6]]

# lapply(pcaNewList, FUN = randForest_fctn, algorithm = "PCA")
# saveRDS(pcaDataMlResList, "processed/pcaDataMlResList.rds")
# saveRDS(mlResList, "processed/test_03.rds")
# saveRDS(nmfDataMlResList, "processed/tempMLResults.rds")
# pcaProcessed[[1]]
xgboost_fctn <- function(preprocessResult, algorithm) {
  if(algorithm == "NMF"){
  model_xgb <- boost_tree(trees = 1000, tree_depth = tune(), mtry=tune(), min_n = tune(), learn_rate = tune(),
                          sample_size = tune(),
                          loss_reduction = tune(), mode = "classification") %>%
    set_engine("xgboost")
  set.seed(234)
  folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 5)
  
  #use grid latin ihiypercube because this covers
  xgbGrid <- grid_latin_hypercube(
    tree_depth(),
    min_n(),
    learn_rate(),
    loss_reduction(),
    sample_size = sample_prop(),
    finalize(mtry(), preprocessResult$train),
    size = 100
  )
  
  wkflow<- workflow() %>%
    add_recipe(preprocessResult$recipe) %>%
    add_model(model_xgb)
  
  res <- tune_grid(
    wkflow,
    resamples = folds,
    grid = xgbGrid,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )
  
  crossValMetrics <- do.call("rbind", res$.metrics)
  
  metrics_xgb <- res %>% collect_metrics()
  final_model <- res %>% select_best(metric = "roc_auc")
  
  finalised_wf <- finalize_workflow(wkflow, final_model)
  
  final_fit <- finalize_workflow(wkflow, final_model) %>%
    fit(data = preprocessResult$train) 
  
  aug <- augment(final_fit, preprocessResult$test)
  aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
  
  roc_auc <- roc_auc(aug, significant, .pred_FALSE)
  two_classCurve <- roc_curve(aug, truth = significant,
                              .pred_FALSE)
  rocCurve <- autoplot(two_classCurve)

  return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, 
              "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
              "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat, 
              "aug"= aug_m,"confMat"= confMat, "algorithm"=algorithm, "dim" = dimension))
  }
}


nmf150_reduced <- nmfDataframes[[5]] %>% select(V110, V128, V40, V73, V86, significant)
splitnmf150 <- split_processingData_fctn(nmf150_reduced, .8)
plan(multisession, workers=30)
reduceddfMLRES <- justToTestRF(splitnmf150, "NMF")
plan(multisession, workers=availableCores())
xgbRES <- justtotestXGB(split10dim, "NMF")
#-----------------------------------------------------------------------------------------------------------------------

# ######this all for 04_modelInterpretation script#####
# 
# # need to read in results from last script
# # nmfRes <- readRDS("processed/tempMLResults.rds")
# # read in also the list of different dimensions:
# # nmfDataframes <- readRDS("processed/nmfDataframes.rds")
# 
# #### get all the auc scores ####
# 
# 
# auc <- lapply(nmfDataMlResList, function(x){
#   x$AUC$.estimate
# })
# 
# #convert each one to df
# df_auc <- lapply(auc, function(x){
#   df <- data.frame(auc_scores = x)
# })
# 
# 
# aucDf <-do.call("rbind", df_auc)
# 
# 
# 
# #remember the dimensions that were input to model to give these outputs:
# dimensions <- data.frame(dimensions =c(5,10), Algorithm = "NMF")
# resDf <- cbind(aucDf, dimensions)
# ggplot(resDf, aes(x=dimensions, y=auc_scores))+
#   geom_point() +
#   theme_classic()
# 
# bestDimension <- resDf[which.max(resDf$auc),]$dimensions
# ##here need to change because want tot loop through the nmf matrices
# #not the ml results
# for (d in 1:length(nmfRes)){
#   #print(nmfRes[[d]])
#   if(ncol(nmfRes[[d]]) == bestDimension){
#     print(nmfRes[[d]])
#     bestDf <- as.data.frame(nmfRes[[d]])
#   }
# }
# 
# # bestDf <- pcList[[5]]
# 
# vipDf <- data.frame(Variable = c("PC5", "PC1", "PC4","PC3","PC2"), Importance= c(4.76,4.37,4.20,-1.10,-1.23))
# # row.names(tempVipDf) <- tempVipDf$Variable
# 
# vec <- c(rownames(vipDf))
# vec2 <- c(vipDf$Variable)
# relevant <- bestDf %>% select(all_of(c(vipDf$Variable)))
# 
# # ^these are for GSEA
# 
# ####now part 2 of this script is to make predictioins on the unlabelled genes:
# # first need to read in the unlabelled genes
# unlabelledGenes <- readRDS("processed/unlabelledGenes.rds")
# unlabelledGenes$prediction <- predict(final_fit, unlabelledGenes)
# 
# # ctrl+shift+C
# 
# 
# #### some testing on just one dataset####
# 
# 
# 
# str()
# temp <- labelledGenesPCARes[[1]]
# #temp <- temp %>% select(1,4,5,6,7,10)
# 
# #convert outcome to a factor
# #temp$significant <- as.factor(temp$significant)
# str(temp)
# skim(temp)
# #set seed when creating split
# dataSplit <- initial_split(temp, prop = 0.8)
# tempTrain <- training(dataSplit)
# tempTest <- testing(dataSplit)
# 
# #look at before sampling
# ggplot(tempTrain, aes(factor(significant)))+
#   geom_bar(aes(y=after_stat(count)/sum(after_stat(count))), colour = "black", fill = "lightgrey")+
#   scale_y_continuous(labels = percent)+
#   xlab("") + ylab("% of genes") +
#   theme_minimal(base_size = 20)
# 
# 
# temp_rec <- recipe(significant ~ ., data = tempTrain) %>%
#   step_downsample(significant, under_ratio =1)
# 
# temp_rec %>% prep() %>% bake(NULL)
# 
# 
# 
# #how does undersampling look after
# temp_rec %>%
#   prep() %>%
#   juice() %>%
#   ggplot(aes(factor(significant))) +
#   geom_bar(aes(y = (after_stat(count))/sum(after_stat(count))), colour="black",fill="lightgrey") +
#   scale_y_continuous(labels = percent) +
#   xlab("") + ylab("% of genes") +
#   theme_minimal(base_size = 20)
# 
# 
# # set seed for CV:
# ?yardstick::metric_set
# cv_folds <- vfold_cv(data = tempTrain, v = 5, repeats = 5)
# 
# #model specification: RF
# 
# show_engines("rand_forest")
# rf_mod <- rand_forest(trees = 500,
#                       mtry = tune(),min_n = tune(),
#                       mode = "classification") %>%
#   set_engine("randomForest", importance = TRUE)
# 
# rf_tune_grid <- grid_regular(
#   mtry(range = c(1,6)),
#   min_n(range = c(1,2)),
#   levels = 6
# )
# ##build workflow:
# temp_wf <- workflow() %>%
#   add_recipe(temp_rec) %>%
#   add_model(rf_mod)
# 
# #doParallel::registerDoParallel()
# 
# #tune model:
# tune_rf <- tune_grid(
#   temp_wf,
#   resamples = cv_folds,
#   grid = rf_tune_grid,
#   control = control_grid(save_pred = TRUE), parallel_over="everything",
#   metrics = metric_set(roc_auc)
# ) #had to install package randomForest for this to work
# ##now fit model on
# #can look at the roc  for each combination of parameters
# autoplot(tune_rf)
# #the exact same plot above:
# plot <- tune_rf%>%
#   collect_metrics()%>%
#   ggplot(aes(x=mtry, y=mean))+
#   geom_point()+
#   geom_line()
# plot
# 
# 
# best <- tune_rf %>% collect_metrics() %>%
#   arrange(.metric)
# #collect_metrics(tune_rf)
# #select_best(tune_rf, metric = "roc_auc")
# #############################fit_mod <- best_fit(temp_wf, split = dataSplit, metric_set(roc_auc))
# best_mod <- tune_rf %>% select_best(metric = "roc_auc")
# 
# #fit the final model with best params on training
# final_fit <- finalize_workflow(temp_wf, best_mod) %>%
#   fit(data = tempTrain)
# 
# 
# 
# #final_fit%>% extract_fit_parsnip()%>%tidy()
# #ff <- finalize_workflow(temp_wf, best_mod)
# ffff <- best_mod %>% last_fit(split = dataSplit)
# 
# 
# 
# final_fit %>% extract_fit_parsnip() %>%
#   vip(geom = 'col',horizontal=T, num_features=5,
#       aesthetics = list(colour="black", fill ="lightblue", alpha=0.7)) +
#   theme_classic()
# #also can get a tibble of variable importances! but will change to df
# final_fit %>% extract_fit_parsnip() %>%
#   vi() %>% as.data.frame()
# 
# 
# 
# #get the test performance
# temp_aug <- augment(final_fit, tempTest)
# head(temp_aug)
# roc_auc(temp_aug, significant, .pred_TRUE)
# conf_mat(temp_aug, significant, .pred_class)
# 
# #plot ROC curve:
# two_classCurve <- roc_curve(temp_aug, truth = significant, .pred_TRUE)
# autoplot(two_classCurve)
# roc_auc(temp_aug, significant, .pred_FALSE)
# 
