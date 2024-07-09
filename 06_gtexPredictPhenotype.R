#this script is testing to see if pipeline can predict skeletal phenotypes in
#a different gene expression dataset

library(tidyverse)

# library(preprocessCore)
# library(future)
# library(future)

gtex <- read.delim("raw/rnaSeqGtex/gtex_Analysis.gct", header = T,skip = 2, sep = "\t")

#just retain protein coding genes
# 
# split the rownames
geneIDsgtex <- as.data.frame(do.call('rbind', strsplit(as.character(gtex$Name),
                                                       '.', fixed = T)))

#remove first 2 columns now
#gtex <- gtex[,-c(1:2)]

length(unique(geneIDsgtex$V1))

humanProteinCoding <- read.table("processed/human_coding_genes.txt", header = T,
                                 stringsAsFactors = F, sep = "\t")
#

row.names(gtex) <- geneIDsgtex$V1
head(gtex[1:3,1:3])
#also get rid of name and description columns:
gtex <- gtex[-c(1,2)]

gtex$ensembl_gene_id <- rownames(gtex)
gtexProteinCoding <- merge(humanProteinCoding, gtex, by = "ensembl_gene_id")
row.names(gtexProteinCoding)<- gtexProteinCoding$ensembl_gene_id; gtexProteinCoding$ensembl_gene_id <- NULL


#now we have all the gene expression data in ensembl ID format
filtgtex <- gtexProteinCoding[rowSums(gtexProteinCoding)>0,]
rowsFiltgtex <- rownames(filtgtex)
gt <- filtgtex[, colSums(filtgtex)>0]
gtex_L <- log2(gt+1)
#head(filtgtex[1:5,1:5])
gtexNor <- normalize.quantiles(as.matrix(gtex_L))
rownames(gtexNor)<- rowsFiltgtex
gtexNor<- as.data.frame(gtexNor)
# 
# 
# #perform PCA
plan(multisession, workers=availableCores())
pca_gtexT <- prcomp((gtexNor), scale. = T)
pcaGtextpcs <- pca_gtexT$x
dim(pcaGtextpcs)
# pcaGtextpcsLabelledPhenotype 
for (i in 1:length(ranks)) {
  #create empty lists
  gtexpcList_t <- list()
  #loop through each of the dimensions
  for (l in (ranks)) {
    dime <- pcaGtextpcs[,1:l]
    gtexpcList_t <- append(gtexpcList_t, list(dime))
  }
}
#copy labels function
getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}
#
#apply to result
labelledgtexPhenotype <- map(.x = gtexpcList_t, .f = getLabelledGenesFctn, knownLabels = labelsFullDf)

gtex_150dim <- labelledgtexPhenotype[[5]]

#temporary RF function not using too much tuning or repeats of cross-validation just for testing- 
#not the function for real results basically
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
plan(multisession, workers = availableCores())
# gtexTemporaryList <- list()
gtexPrePR <- split_processingData_fctn(data = gtex_150dim,proportion = 0.8)
# 
# gtexTemporaryList[[1]] <- gtexPrePR
# #test whether can detect skeletal phenotype
gtex_150dim_res <- justToTestRF(gtexPrePR, algorithm = "PCA")
# gtexMLRES <- map(.x = gtexTemporaryList, .f = justToTestRF, algorithm="PCA")
# #looks like it doesn't
