#this script is testing to see if pipeline can predict skeletal phenotypes in
#a different gene expression dataset
# but don't really need now cause i have another positive control

library(tidyverse)

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
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE, 
                                                                    event_level="second" )
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=3, repeats = 2)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      trees(range = c(1,1000)),
      mtry(range = c(1,2)),
      min_n(range = c(1,2)),
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
    options(event_LEVEL = FALSE)
    
    confMat <- aug %>%
      conf_mat(truth = significant, estimate= .pred_class)
    
    roc_auc <- roc_auc(aug, significant, .pred_TRUE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_TRUE, event_level="second")
    rocCurve <- autoplot(two_classCurve)
    dimension <- ncol(preprocessResult$train)-1
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    # saveRDS(result, sprintf("%sGTEXtestingMLRes_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
    return(list("workflow" = wkflow, "res" = res, "dfForCorrectR" = df,"resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve,"confMat" = confMat, "aug"= aug_m,"confMat"= confMat, "algorithm"=algorithm, "dim" = dimension))
    
  }
  else if(algorithm == "PCA"){
    model_RF <- rand_forest(trees = 500, mtry = sqrt(ncol(preprocessResult$train)), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=2)
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
    confMat <- aug %>%
      conf_mat(truth = significant, estimate= .pred_class)
    
    roc_auc <- roc_auc(aug, significant, .pred_FALSE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_FALSE,event_level="first")
    rocCurve <- autoplot(two_classCurve) 
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    #saveRDS(result, sprintf("processed/%sGTEXtestingMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    dimension <- ncol(preprocessResult$train)-1
    return(list("workflow" = wkflow, "res" = res, "dfForCorrectR" = df, "resDF"=resDf, "finalMod" = final_model, 
                "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m, "confMat" = confMat, 
                "algorithm"=algorithm, "dim" = dimension))
  }
  
}
#
plan(multisession, workers = availableCores())
# gtexTemporaryList <- list()
gtexPrePR <- split_processingData_fctn(data = gtex_150dim,proportion = 0.8)
# 
plan(multisession, workers=20)
nmfDataSplit <- map(.x = nmfData, .f = split_processingData_fctn, 0.8)
nmResMLall <- map(.x = nmfDataSplit, .f = justToTestRF, algorithm = "NMF")
nmf_res_150 <- justToTestRF(nmfDataSplit[[5]], algorithm = "NMF")
# gtexTemporaryList[[1]] <- gtexPrePR
# #test whether can detect skeletal phenotype
gtex_150dim_res <- justToTestRF(gtexPrePR, algorithm = "PCA")
# gtexMLRES <- map(.x = gtexTemporaryList, .f = justToTestRF, algorithm="PCA")
# #looks like it doesn't



# ####use gtex to predict sex####
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
