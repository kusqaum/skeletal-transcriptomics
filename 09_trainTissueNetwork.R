library(tidyverse)
library(tidymodels)
library(themis)
library(vip)
library(doParallel)
library(foreach)


networkAlone <- readRDS("processed/processedNetworkEmb.rds")
networkAloneLabelled <- readRDS("processed/processedNetworkEmbLabelled.rds")
sVwithNetworkLabelled <- readRDS("processed/SVwithNetworkLabelled.rds")

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

xgboost_fctnNetwork <- function(preprocessResult, algorithm) {
  paste0("Number of predictors = ",print(ncol(preprocessResult$train)-1))
  if(algorithm == "none"){
    
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
    saveRDS(result, sprintf("processed/xgbResNetworkOnly_%s.rds", ncol(preprocessResult$train)-1))
    
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
    saveRDS(result, sprintf("processed/%sxgbRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model,
                "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat,
                "aug"= aug_m,"confMat"= confMat, "algorithm"=algorithm, "dim" = dimension, "dfForCorrectR" = df))
  }
}


# cl <- makeCluster(25, outfile="processed/200dimNMF_NetworkRF.txt")
# registerDoParallel(cl)
# print("randforest on 200 dim NMF SkeletalVis integrated with network:")
# nmfNetworkdim200Res <- randForest_fctnWithNetwork(sVwithNetworkSplit[[6]], algorithm = "NMF")
# 
c <- makeCluster(25, outfile="processed/xgbNetworkOnly.txt")
registerDoParallel(c)
print("xgb on network data alone")
netRes <- xgboost_fctnNetwork(networkAloneSplit ,algorithm = "none")
