library(tidyverse)
library(tidymodels)
library(themis)
library(skimr)
library(randomForest)
library(vip)
library(xgboost)
library(foreach)
library(doParallel)
# this is just as a sanity check (+ve control) to make sure that the pipeline can predict 
# another phenotype from the IMPC 
#function to label genes
getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}
#function to split data and apply preprocessing steps
split_processingData_fctn <- function(data, proportion){
  set.seed(123)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  rec <-recipe(significant~., data = trainData) %>%
    step_downsample(significant, under_ratio = 1, seed = 456)
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}

geneExpress <- readRDS("processed/unlabelledNMFDfs.rds")
head(geneExpress[[1]][1:5,1:3])
#read in the mortality/aging associated genes
mortalityLabels <- read.table("processed/processedMortalityLabels.txt", header = T, sep = "\t")
mortalityLabels$significant <- as.factor(mortalityLabels$significant)
#add on the labels to the nmf dataset
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


res50dim <- randomForestFctn_positiveControl(mortSplit[[3]], algorithm = "NMF")
res100dim <- randomForestFctn_positiveControl(mortSplit[[4]], algorithm = "NMF")
res150dim <- randomForestFctn_positiveControl(mortSplit[[5]], algorithm = "NMF")
res200dim <- randomForestFctn_positiveControl(mortSplit[[6]], algorithm = "NMF")
# print("finished")


##
#now read in the ML results:

mortFiles <- list.files("processed", pattern = "ResSVmortality", full.names = T)
mortalityRes <- lapply(mortFiles, readRDS)
mortalityRes <- mortalityRes[order(sapply(mortalityRes, function(x) x$dim))]

#look at their results on held-out set
mortPred <- lapply(mortalityRes, function(x){
  x$aug
})
#and plot:
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

# get the auc values to add to the plot as well
mortalityAUC <- lapply(mortalityRes, function(x){
  x$AUC
})
mortalityAUC <- do.call("rbind",mortalityAUC)
  
ggsave("output/mortalityCurves.png", mortalityROC, width = 14.5, height = 15)
