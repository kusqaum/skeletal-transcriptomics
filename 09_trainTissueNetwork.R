library(tidyverse)
library(tidymodels)
library(themis)
library(vip)
library(doParallel)
library(foreach)
library(RColorBrewer)

##### read in data sets####
networkAlone <- readRDS("processed/processedNetworkEmb.rds")
networkAloneLabelled <- readRDS("processed/processedNetworkEmbLabelled.rds")
sVwithNetworkLabelled <- readRDS("processed/SVwithNetworkLabelled.rds")

labelsAlone <- sVwithNetworkLabelled[[1]] %>% select(significant)

head(sVwithNetworkLabelled[[1]][1:5,1:5])
dim(sVwithNetworkLabelled[[1]])
#### split data and create preprocessing rec ####
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
networkAloneSplit <- split_processingData_fctn(networkAloneLabelled, 0.8)
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
    else if (ncol(preprocessResult$train)-1 <500){
      levels <- 18
    }
    else if(ncol(preprocessResult$train)-1 >=500){
      levels <- 15
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

xgboost_fctnNetworkAlone <- function(preprocessResult){
  paste0("Number of predictors = ",print(ncol(preprocessResult$train)-1))

  model_xgb <- boost_tree(trees = tune(), tree_depth = tune(),
                          mtry=tune(), min_n = tune(), 
                          learn_rate = tune(),
                          sample_size = tune(),
                          loss_reduction = tune(), mode = "classification") %>%
    set_engine("xgboost")
  
  set.seed(234)
  folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 3)
  #use grid latin ihiypercube because this covers
  xgbGrid <- grid_latin_hypercube(
    trees(),
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
    control = control_grid(save_pred = TRUE, verbose = TRUE),
    metrics = metric_set(roc_auc)
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
  df <- df%>% mutate(model = "XGB")
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
  metrics_xgb <- res %>% collect_metrics()
  finalised_wf <- finalize_workflow(wkflow, final_model)
  
  aug <- augment(final_fit, preprocessResult$test)
  aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
  confMat <- aug %>%
    conf_mat(truth = significant, estimate= .pred_class)
  roc_auc <- roc_auc(aug, significant, .pred_FALSE)
  two_classCurve <- roc_curve(aug, truth = significant,
                              .pred_FALSE)
  rocCurve <- autoplot(two_classCurve)
  dimension <- ncol(preprocessResult$train)-1
  result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
                 "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                 "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
                 "aug"= aug_m, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df)
  
  saveRDS(result, sprintf("processed/xgbResNetworkOnly_%s.rds", ncol(preprocessResult$train)-1))
  
  return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
              "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
              "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
              "aug"= aug_m,"confMat"= confMat, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df))
  }

xgboost_fctnNetwork <- function(preprocessResult, algorithm) {
  paste0("Number of predictors = ",print(ncol(preprocessResult$train)-1))
  if(algorithm == "NMF"){
    
    model_xgb <- boost_tree(trees = tune(), tree_depth = tune(), 
                            mtry=tune(), min_n = tune(), 
                            learn_rate = tune(),
                            sample_size = tune(),
                            loss_reduction = tune(), mode = "classification") %>%
      set_engine("xgboost")
    set.seed(234)
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 3)
    
    #use grid latin ihiypercube because this covers
    xgbGrid <- grid_latin_hypercube(
      trees(),
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
      control = control_grid(save_pred = TRUE, verbose = TRUE),
      metrics = metric_set(roc_auc)
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
    df <- df%>% mutate(model = "XGB")
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
    
    metrics_xgb <- res %>% collect_metrics()
    final_model <- res %>% select_best(metric = "roc_auc")
    
    finalised_wf <- finalize_workflow(wkflow, final_model)
    
    final_fit <- finalize_workflow(wkflow, final_model) %>%
      parsnip::fit(data = preprocessResult$train)
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    
    confMat <- aug %>%
      conf_mat(truth = significant, estimate= .pred_class)
    
    roc_auc <- roc_auc(aug, significant, .pred_FALSE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_FALSE)
    rocCurve <- autoplot(two_classCurve)
    dimension <- ncol(preprocessResult$train)-1
    
    result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
                   "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
                   "aug"= aug_m, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df)
    saveRDS(result, sprintf("processed/%sxgbResSV_%s.rds", ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
                "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
                "aug"= aug_m,"confMat"= confMat, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df))
  }
  
  else if(algorithm == "PCA"){
    model_xgb <- boost_tree(trees = tune(), tree_depth = tune(), mtry=tune(), min_n = tune(),
                            learn_rate = tune(),
                            sample_size = tune(),
                            loss_reduction = tune(), mode = "classification") %>%
      set_engine("xgboost")
    set.seed(234)
    folds <- vfold_cv(data = preprocessResult$train, v=5, repeats = 3)
    
    #use grid latin ihiypercube because this covers
    xgbGrid <- grid_latin_hypercube(
      trees(),
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
      control = control_grid(save_pred = TRUE, verbose = TRUE),
      metrics = metric_set(roc_auc)
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
    df <- df%>% mutate(model = "XGB")
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
    
    metrics_xgb <- res %>% collect_metrics()
    
    finalised_wf <- finalize_workflow(wkflow, final_model)
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train)-1, Algorithm = algorithm)
    
    confMat <- aug %>%
      conf_mat(truth = significant, estimate= .pred_class)
    
    roc_auc <- roc_auc(aug, significant, .pred_FALSE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_FALSE)
    rocCurve <- autoplot(two_classCurve)
    dimension <- ncol(preprocessResult$train)-1
    result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
                   "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
                   "aug"= aug_m, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df)
    saveRDS(result, sprintf("processed/%sxgbRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
                "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
                "aug"= aug_m,"confMat"= confMat, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df))
  }
}


cl <- makeCluster(25, outfile="processed/500dimNMF_NetworkRF.txt")
registerDoParallel(cl)
print("randforest on 500 dim NMF SkeletalVis integrated with network:")
nmfNetworkdim500Res <- randForest_fctnWithNetwork(sVwithNetworkSplit[[7]], algorithm = "NMF")
# 
# c <- makeCluster(25, outfile="processed/xgbNetworkOnly.txt")
# registerDoParallel(c)
# print("xgb on network data alone")
# netRes <- xgboost_fctnNetwork(networkAloneSplit ,algorithm = "none")

#-------------------------------------------------------------------------------------------
# read in the results of network alone!
networkOnlyXGBres <- readRDS("processed/xgbResNetworkOnly_500.rds")
networkOnlyXGBres$AUC
#0.515
networkOnlyXGBres$roc_curve

###read in the results of tissue network
rfNetworkPlusSVnmfFiles <- list.files("processed", full.names = T, pattern = "ResNetworkWithSV_")
rfNetworkPlusSVnmfResList <- lapply(rfNetworkPlusSVnmfFiles, readRDS)

rfNetworkPlusSVnmfResList <- rfNetworkPlusSVnmfResList[order(sapply(rfNetworkPlusSVnmfResList, function(x) x$dim))]
cvRFNetSV <- lapply(rfNetworkPlusSVnmfResList, function(x){
  x$resDf
})
#nonnetwork for comparison
cvNonetwork <- readRDS("processed/cvAllRf_noNetworkRes.rds")
cvNonetwork <- cvNonetwork %>% mutate(Feature = "SkeletalVis only")

cvNetSv <- do.call("rbind", cvRFNetSV) %>% mutate(model = "RF", Feature = "SkeletalVis With network")
networkandNoNetwork <- bind_rows(cvNonetwork, cvNetSv)
networkSVcvBP <- ggplot(networkandNoNetwork, aes(x= as.factor(dim), y = mean, fill=Feature)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#4DBBD5B2", "#7E6148B2"))+
  labs(x="Number of NMF dimensions", y="Cross-validation AUC", fill= "") +
  theme_cowplot(font_size = 26)+
  facet_wrap(~"RF")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4), 
        strip.background = element_rect(color = "black",  linewidth = 0.4),
        legend.position = "bottom")
networkSVcvBP

p1 <- readRDS("output/p1_RFvsXGB.rds")
p2 <- readRDS("output/p2_pcaVsNmf.rds")

bot <- plot_grid(p1,p2, labels = c("b","c"), label_size = 16)
bot
# allp <- plot_grid(networkSVcvBP / (p1+p2), labels = "auto", label_size = 20)
allp <- plot_grid(networkSVcvBP, bot, ncol = 1, labels = "auto", label_size = 16)
allp
ggsave("output/allp.png", allp, height = 10, width = 11)
#
#
cvRfMetrNET_SV <- lapply(rfNetworkPlusSVnmfResList, function(x){
  r <- x$resDf
  maxAuc <- as.data.frame(r[which.max(r$mean),])
})
maxAucNetDims <- do.call("rbind", cvRfMetrNET_SV) %>%
  select(mean, dim, lowerbound, upperbound, Algorithm)
## although the cv metrics are quite better on the training data, they do not reflect as 
#well on the test data:...



predSVWithNet <- lapply(rfNetworkPlusSVnmfResList, function(x){
  x$aug
})

aucWithNet <-lapply(rfNetworkPlusSVnmfResList, function(x){ # where the input is a list containing multiple ML results
  x$AUC$.estimate
})
dfAUCwNet <- lapply(aucWithNet, function(x){
  df<- data.frame(auc_scores = x)
})
dfAUCwNet <- do.call("rbind", dfAUCwNet)
dimens <- lapply(rfNetworkPlusSVnmfResList, function(x){
  x$dim
})
dimens <- unlist(dimens)
dfAUCwNet_text <- data.frame(dim = dimens, 
                             auc = dfAUCwNet$auc_scores, Algorithm = "NMF")
rfSVWithNetpredictionsMetr <- bind_rows(predSVWithNet) %>%
  group_by(dim) %>%
  roc_curve(truth = significant, .pred_TRUE, event_level = "second") %>%
  mutate(alg = paste(dim, "NMF dimensions"))

str(rfSVWithNetpredictionsMetr$alg)
alg <- "NMF dimensions"
svWnetROCrf <- ggplot(rfSVWithNetpredictionsMetr, aes(x=1-specificity, y=sensitivity, colour=as.factor(dim)))+
  geom_line(linewidth=1.5, show.legend = F)+ 
  geom_abline(slope = 1, intercept = 0, linewidth=0.4, lty="dashed", alpha = 0.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 1.0, fill="white"),
        aspect.ratio = 1, legend.position="none")+
  theme_minimal_grid(font_size = 30)+# theme(legend.position = "none")+
  scale_colour_npg()+
  # facet_wrap(~ paste0(dim, " NMF dimensions"))+
  facet_wrap(~ factor(paste0(dim, " NMF dimensions"), c("50 NMF dimensions", "100 NMF dimensions","150 NMF dimensions",
                                                        "200 NMF dimensions", "500 NMF dimensions")))+
  theme(strip.text = element_text(size=32), legend.position = "none")+
  geom_text(data = dfAUCwNet_text, mapping = aes(x=0.3, y=0.85, label = paste0("AUC = ", round(auc, 3))), 
            size=12)

svWnetROCrf
ggsave("output/SVnetworkROCcurveRF.png", svWnetROCrf, width = 20, height = 12)


posBestFitNet <- which.max(dfAUCwNet_text$auc)
posBestFitNet
bestDimensionNet <- dfAUCwNet_text[which.max(dfAUCwNet_text$auc),]$dim
bestDimensionNet
dataForTraining <- sVwithNetworkLabelled
for (n in 1:length(rfNetworkPlusSVnmfResList)){
  if(ncol(dataForTraining[[n]])-1 == bestDimensionNet){#allnmfmatrix is a list containing all the nmf reduced matrices +labels
    bestDfwNet <- as.data.frame(dataForTraining[[n]])
  }
}
dim(bestDfwNet)
bestFitNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$finalFit
bestFitNet

correctRdfBestDimNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$dfForCorrectR

bestModMetricsNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$resDf
finalModMetricsNet <- bestModMetricsNet[which.max(bestModMetricsNet$mean),]




#lets look at the variable importance
bestModVarImportNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$importanceDf
rfNetworkPlusSVnmfResList[[posBestFitNet]]$importancePlot
head(bestModVarImportNet)

bestModVarImportNet <- bestModVarImportNet %>% 
  mutate(sign = case_when(Importance<0 ~"negative", TRUE~"positive"))
head(bestModVarImportNet)
dim(bestModVarImportNet)
top5FeatsNet <- bestModVarImportNet[1:5,]
fiveVarsNet <- top5FeatsNet$Variable



top5FeatsModeldfNet <- bestDfwNet %>% dplyr::select(all_of(fiveVarsNet))
head(top5FeatsModeldfNet)
colnames(top5FeatsModeldfNet)<- sub("X", "Feature", colnames(top5FeatsModeldfNet))
head(top5FeatsModeldfNet)


dim(top5FeatsModeldfNet)
# labelsFullDf <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
# str(labelsFullDf$significant)
labelsAlone$significant<- as.factor(labelsAlone$significant)
all(rownames(top5FeatsModeldfNet)%in% rownames(labelsAlone))

all(rownames(top5FeatsModeldfNet)== rownames(labelsAlone))
#they're in order so can just do cbind
top5FeatsModeldfLabelledNet <- cbind(top5FeatsModeldfNet, labelsAlone)

# df <- top5FeatsModeldfLabelledNet[order(top5FeatsModeldfLabelledNet$Feature76, decreasing = T),]
# ggplot(df[1:15,], aes(x=significant, y = Feature76))+
#   geom_boxplot()

for (m in 1:(ncol(top5FeatsModeldfLabelledNet)-1)) {
  #print(m)
  feature <- colnames(top5FeatsModeldfLabelledNet)[m]
  n <- sprintf("%s_nmfSvWithNet_%s",m,colnames(top5FeatsModeldfLabelledNet)[m])
  p <- ggplot(top5FeatsModeldfLabelledNet)+ 
    aes(y = top5FeatsModeldfLabelledNet[,m], x= significant, fill=significant)+
    geom_boxplot()+
    xlab("") + ylab(paste0("NMF " ,colnames(top5FeatsModeldfLabelledNet)[m]))+
    theme_cowplot(font_size = 16)+
    theme(panel.background = element_rect(colour = "black", fill = NA, linewidth = 0.4))+
    # scale_fill_bmj()+
    scale_fill_manual(values = c("#69BE28B2", "#E37222B2"))+ 
    scale_x_discrete(labels = c("Associated", "Not associated"))+
    theme(legend.position = "none")
  # ggsave(paste0("output/", n,".png"), p, height = 3, width = 4)
  print(n)
  print(p)
}





top5FeatsModeldfNet <- top5FeatsModeldfNet %>% arrange(desc(Feature76))
hGenesSymbs <- read.table("processed/human_coding_genes.txt", sep = "\t", header = T)
hGenesSymbs <- hGenesSymbs %>% 
  filter(hgnc_symbol!=""); rownames(hGenesSymbs) <- hGenesSymbs$ensembl_gene_id
hGenesSymbs$ensembl_gene_id <-NULL
head(hGenesSymbs) 
# so now gonna merge with top feats DF
head(top5FeatsModeldfLabelledNet)
topFeatsMapped <- merge(top5FeatsModeldfNet, hGenesSymbs, by=0); rownames(topFeatsMapped) <- topFeatsMapped$hgnc_symbol; topFeatsMapped$Row.names<-NULL; topFeatsMapped$hgnc_symbol <-NULL
head(topFeatsMapped)
topFeatsMapped <- topFeatsMapped %>% arrange(desc(Feature76))

matN <- as.matrix(topFeatsMapped)
head(matN)
matN <- apply(matN, 2, rank)
head(matN)
brewer.pal.info
heatmapN <- pheatmap::pheatmap(matN[1:15,], border_color = "white",
                              cluster_rows = F, 
                              cluster_cols = F, 
                              show_rownames = T, 
                              fontsize = 20, 
                              color = brewer.pal(8, "Reds")
                              
)

heatmapN
ggsave("output/heatmap.png", heatmapN, width = 10, height = 12)
ggsave("output/heatmap.png", heatmapN, width = 11, height = 6)

####clusterprofiler code####

getEachFeature_fctn <- function(dataframeOfFeats){
  single <- list()
  for (j in 1:ncol(dataframeOfFeats)){ 
    print(j)
    single[[j]] <- data.frame(dataframeOfFeats[,j],
                              row.names = rownames(dataframeOfFeats))
    colnames(single[[j]]) <- colnames(dataframeOfFeats)[j]
    single[[j]] <- single[[j]] %>% arrange(desc(colnames(single[[j]])))
  }
  return(single)
} 

eachNetFeat <- getEachFeature_fctn(top5FeatsModeldfNet)
saveRDS(eachNetFeat, "processed/top5NetFeats.rds") # do the gsea on my laptop cause clusterprofiler 
#not installing




############
colnamesNet <- paste0("Feature", 1:ncol(bestDfwNet))
unlabelledGenes <- readRDS("processed/unStudiedGenes.rds")
head(unlabelledGenes[1:4,1:5])
colnames(unlabelledGenes) <- sub('V', 'X', colnames(unlabelledGenes))
head(unlabelledGenes[1:4,1:3])
unlabelledGenesPredNet <- augment(bestFitNet, unlabelledGenes)


unlabelledGenesPreddfNet <- as.data.frame(unlabelledGenesPredNet); rownames(unlabelledGenesPreddfNet) <- rownames(unlabelledGenesPredNet)
head(unlabelledGenesPreddfNet[1:3,1:4])
length(which(unlabelledGenesPreddfNet$.pred_class=="TRUE"))
length(which(unlabelledGenesPreddfNet$.pred_class!= "TRUE"))


########

head(unlabelledGenesPreddfNet[1:4,1:5])
mgiGenes <- read.table("processed/annotatedMGIgenes.txt", header = T, sep = "\t")
dim(mgiGenes)
#just remove duplicates since they're all the same level anyways
mgiGenes <- mgiGenes %>% distinct(ensembl_gene_id, .keep_all = T)
dim(mgiGenes)
str(mgiGenes$ensembl_gene_id)
# make gene names rownames
rownames(mgiGenes)<- mgiGenes$ensembl_gene_id ; mgiGenes$ensembl_gene_id<-NULL
head(mgiGenes)
mgiGenes$significant<- as.factor(mgiGenes$significant)
head(unlabelledGenes[1:3,1:5])
# let's just get the genes in both datasets

length(which(rownames(unlabelledGenesPreddfNet) %in% rownames(mgiGenes)))
# so 882 of the mgi genes are known to be associated to a skeletal phenotype
nrow(unlabelledGenesPreddfNet) - (length(which(rownames(unlabelledGenesPreddfNet) %in% rownames(mgiGenes))))
# and 8119 genes are unstudied
#which are those genes then? that aren't in the mgi genes (i.e., unstudied)
totallyUnstudied <- rownames(unlabelledGenesPreddfNet)[!rownames(unlabelledGenesPreddfNet) %in% rownames(mgiGenes)]
class(totallyUnstudied)
length(totallyUnstudied)

unstudied <- data.frame(significant = rep(c("FALSE"), times = length(totallyUnstudied)), row.names = totallyUnstudied)
head(unstudied) # nice
tail(unstudied)
str(unstudied$significant)
unstudied$significant <- as.factor(unstudied$significant)
nrow(unstudied)
nrow(mgiGenes) # cool
mgiGenesPlusUnstudied <- rbind(mgiGenes, unstudied)
# now let's retain only those genes in mgiGenesplus unstudied DF that are also in our unstudied/unlabelled genes
unlabelledGenesWithMGIannot <- merge(unlabelledGenesPreddfNet,mgiGenesPlusUnstudied, by=0); rownames(unlabelledGenesWithMGIannot)<- unlabelledGenesWithMGIannot$Row.names; unlabelledGenesWithMGIannot$Row.names<- NULL
# cool
head(unlabelledGenesWithMGIannot[1:3,1:3])
dim(unlabelledGenesWithMGIannot)

# unlabelledGenesWithMGIannotTEST <- unlabelledGenesWithMGIannot
# unlabelledGenesWithMGIannotTEST$significant <- unlabelledGenesWithMGIannotTEST$.pred_class
mgiCurve <-roc_curve(unlabelledGenesWithMGIannot, truth = significant, .pred_FALSE)

head(mgiCurve)
mgiCurve <- mgiCurve %>% mutate(database = "MGI")
head(mgiCurve)
rocaucMGI <- roc_auc(unlabelledGenesWithMGIannot, truth = significant, .pred_FALSE)
rocaucMGI <- rocaucMGI %>% mutate(database = "MGI")
autoplot(mgiCurve)

# now looking at human phenotype ontology- genes that are annotated to a skeletal abnormality
# lets quickly do this 
hpOnt <- read.table("processed/humphenotOntGenes.txt", header = T, sep = "\t")
head(hpOnt)
length(unique(hpOnt$ensembl_gene_id))

rownames(hpOnt)<- hpOnt$ensembl_gene_id; hpOnt$ensembl_gene_id<-NULL
hpOnt$significant <- as.factor(hpOnt$significant)
length(which(rownames(unlabelledGenesPreddfNet) %in% rownames(hpOnt)))
# we have 1415 genes that are unlabelled (unstudied in impc) to be known to have skel phenotype according to hpo
length(rownames(unlabelledGenesPreddfNet)[!rownames(unlabelledGenesPreddfNet) %in% rownames(hpOnt)])
# the remainder of these 7586genes^ are not known to have a skeletal abnormality in the human phenotype ont
unstGene <- rownames(unlabelledGenesPreddfNet)[!rownames(unlabelledGenesPreddfNet) %in% rownames(hpOnt)]
# so lets get all their gene names and make them as a dataframe and call them false (i.e., not associated with a phenotype)
unstGeneDf <- data.frame(significant = rep(c("FALSE"), times = length(unstGene)), row.names = unstGene)
head(unstGeneDf)
dim(unstGeneDf)
#ncie
dim(hpOnt)
# now let's merge the positively associated genes from hpo and the 'negatively'associated genes#
#we have that are not labelled to be associated with HPO
hpoGenesPlusUnstudied <- rbind(hpOnt, unstGeneDf)
# so now we have + and - labels

head(hpoGenesPlusUnstudied)
dim(hpoGenesPlusUnstudied)

unlabelledGenesWithHPOannot <- merge(unlabelledGenesPreddfNet, hpoGenesPlusUnstudied, by= 0); rownames(unlabelledGenesWithHPOannot) <- unlabelledGenesWithHPOannot$Row.names; unlabelledGenesWithHPOannot$Row.names <- NULL
head(unlabelledGenesWithHPOannot[1:3,1:5])
unlabelledGenesWithHPOannot$significant <-as.factor(unlabelledGenesWithHPOannot$significant)
#cool looks good i guess
hpoCurve <- roc_curve(unlabelledGenesWithHPOannot, truth = significant, .pred_FALSE)
hpoCurve <- hpoCurve %>% mutate(database = "HPO")
rocaucHPO <- roc_auc(unlabelledGenesWithHPOannot, significant, .pred_FALSE)
rocaucHPO <- rocaucHPO %>% mutate(database = "HPO")
rocaucHPO

hpoMgi <- rbind(hpoCurve, mgiCurve)
#let's plot mgi and hpo now
ext_text <- rbind(rocaucHPO, rocaucMGI)

mgihpoROC <- ggplot(hpoMgi, aes(x=1-specificity, y=sensitivity, colour=database)) +
  geom_path(linewidth=0.9)+
  geom_abline(slope = 1, intercept = 0, size=0.4, lty="dashed", alpha = 0.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill="white"),
        aspect.ratio = 1)+
  theme_bw(base_size = 22)+
  theme(legend.position = "none")+
  scale_color_bmj()+
  facet_wrap(~database) + 
  geom_text(data = ext_text, mapping = aes(x=0.3, y=0.87, label = paste0("AUC = ", round(.estimate, 3))),
            size = 8)

mgihpoROC



#### which are the top ranked genes?

orderedGenesNet <- unlabelledGenesPreddfNet[order(unlabelledGenesPreddfNet$.pred_FALSE, decreasing = T),]
top20genesNet <- orderedGenesNet[1:20,]
top20genesWsymbsNet <- merge(top20genesNet, hGenesSymbs, by=0)
rownames(top20genesWsymbsNet) <- top20genesWsymbsNet$hgnc_symbol; top20genesWsymbsNet$Row.names <-NULL; top20genesWsymbsNet$hgnc_symbol<-NULL
top20genesWsymbsNet <- top20genesWsymbsNet[order(top20genesWsymbsNet$.pred_TRUE, decreasing = T),]
head(top20genesWsymbsNet[1:3,1:6])
# rearrange again cause not in order
top20genesWsymbsNet <- top20genesWsymbsNet[order(top20genesWsymbsNet$.pred_FALSE, decreasing = T),]
head(top20genesWsymbsNet[16:20,1:5])
## can save it if you like...
write.table(as.data.frame(top20genesWsymbsNet[,1:3]), "processed/temporary.txt", col.names = T, sep = "\t", quote = F)
