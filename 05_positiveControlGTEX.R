library(tidyverse)
library(tidymodels)
library(themis)
library(skimr)
library(randomForest)
library(vip)
# library(preprocessCore)
# library(future)
# library(future)

gtex <- read.delim("raw/rnaSeqGtex/gtex_Analysis.gct", header = T,skip = 2, sep = "\t")

#just retain protein coding genes
# temp <- gtex[1:5,1:3]
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
head(filtgtex[1:5,1:5])
# gtexNor <- normalize.quantiles(as.matrix(filtgtex))
# rownames(gtexNor)<- rowsFiltgtex
# gtexNor<- as.data.frame(gtexNor)
# 
# 
# #perform PCA
# pca_gtex <- prcomp(gtexProteinCoding, scale. = T)
# pcaGtexpcs <- pca_gtex$x
# for (i in 1:length(ranks)) {
#   #create empty lists
#   gtexpcList <- list()
#   #loop through each of the dimensions
#   for (l in (ranks)) {
#     #grab each dimension
#     dime <- pcsLabelled[,1:l]
#     gtexpcList <- append(gtexpcList, list(dime))
#   }
# }
#add labels function
getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}
# 
# #apply to result
# labelledgtex <- map(.x = gtexpcList, .f = getLabelledGenesFctn, knownLabels = labelsFullDf)

# gtextTest <- map(.x = gtexpcList,.f = performTtestFctn, labels=labelsFullDf$significant, algorithm = "PCA")
# 
# gtex_5dim <- labelledgtex[[5]]
#temporary RF function not using too much tuning or repeats of cross-validation just for testing
gtexRFfctn <- function(preprocessResult, algorithm){
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
    
    roc_auc <- roc_auc(aug, significant, .pred_1)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_1)
    rocCurve <- autoplot(two_classCurve)
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "resDf"= resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    # saveRDS(result, sprintf("%sGTEXMLRes_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
    return(list("workflow" = wkflow, "res" = res, "resDf" = resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
    
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
    
    roc_auc <- roc_auc(aug, significant, .pred_1)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_1)
    rocCurve <- autoplot(two_classCurve)
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    saveRDS(result, sprintf("processed/%sGTEXallMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "resDF"=resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
  }
  
}
#
# plan(multisession, workers = 20)
# gtexTemporaryList <- list()
# gtexPrePR <- split_processingData_fctn(data = gtex_5dim,proportion = 0.8)
# # testFctn <- split_processingData_fctn(data = gtex_5dim, proportion = .8)
# 
# gtexTemporaryList[[1]] <- gtexPrePR
# #test whether can detect skeletal phenotype
# gtexMLRES <- map(.x = gtexTemporaryList, .f = gtexRFfctn, algorithm="PCA")
# #looks like it doesn't

##----------------------------------------------------------------------------------------------------------------------------

sampleNames <- (colnames(gtexProteinCoding))
names <- data.frame(id = sampleNames)
# colnames(namesC)[1]<- id
# df <- (do.call("rbind", strsplit(as.character(sampleNames$colnames(gtexProteinCoding), '.', fixed=T))))

# just get the first part of ensemblID
id <- names %>% mutate(first = unlist(lapply(strsplit(names$id, '\\.'), function(x)x[1])))#
id <- id %>% mutate(second = unlist(lapply(strsplit(names$id, '\\.',), function(x)x[2]))) 
#make it as the same format as the sex is
id$integrate <- paste(id$first, id$second, sep = '-')
dim(id)
colnames(gtexProteinCoding)
#now make these the columnames
colnames(filtgtex) <- id$integrate

gtexPhenotype <- read.table("raw/GTEX_v7_SubjectPhenotypeDS.txt", header = T) # 1 is male, 2 female
length(unique(gtexPhenotype$SUBJID))
# rownames(gtexPhenotype) <- gtexPhenotype$SUBJID
# gtexPhenotypeNew <- data.frame(gtexPhenotype$SEX)
# rownames(gtexPhenotypeNew)<- gtexPhenotype$SUBJID
to_match <- colnames(filtgtex)

findMatch <- match(to_match,gtexPhenotype$SUBJID)
sex <- gtexPhenotype[findMatch, 2]

#gtexTranspose <- as.data.frame(t(filtgtex))                       
filtgtex_t <- as.data.frame(t(filtgtex))
#filtgtex_t$sex <- sexdont do this yet because perform pca first!

filtgtex_tL <- log2(filtgtex_t+1)
head(filtgtex[1:3,1:4])
#normalise quantiles
# filtgtexN <- normalize.quantiles(as.matrix(filtgtex_tL))
# rownames(filtgtexN) <- rownames(filtgtex_tL)
# filtgtexN <- as.data.frame(filtgtexN)
# saveRDS(filtgtexN, "processed/filtgtexN.rds")
filtgtexN <- readRDS("processed/filtgtexN.rds")
ranks <- c(5,10,50,100,150,200)
## now let's perform PCA 
# plan(multisession, workers=availableCores())
print("performing pca")
gtexPCres <- prcomp(filtgtexN, scale. = T)
print("finished PCA")
gtexPCs <- as.data.frame(gtexPCres$x)
gtexPCs$sex <- sex
for (g in 1:length(ranks)) {
  #create empty lists
  gtexpcList <- list()
  #loop through each of the dimensions
  for (l in (ranks)) {
    #grab each dimension
    dime <- gtexPCs[,1:l]
    gtexpcList <- append(gtexpcList, list(dime))
  }
}
#now i need to add on the labels for sex!!!

# rownames(sexDf) <- colnames(filtgtex)

gtexSex <- data.frame(significant = as.factor(gtexPCs[,ncol(gtexPCs)]))
head(gtexSex)
rownames(gtexSex)<- rownames(gtexPCs)

gtexPCsLabelled <- map(.x = gtexpcList , .f =  getLabelledGenesFctn, knownLabels = gtexSex)
("print split data")
splitGTEXdata <- function(data, proportion){
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
print("create preprocessing recipe now")
gtexProcessed <- lapply(X = gtexPCsLabelled, FUN = splitGTEXdata, 0.8)
print("done creating recipe")
print("run ML now")
gtexAllMlRes <- lapply(X=gtexProcessed, FUN = gtexRFfctn, algorithm="PCA")
print("finished ML")




##below is just testing the ML pipeline on a smalll portion of GTEX dataset
#just get a small subset of the gtex data
# smallgtex <- filtgtex[,1:790]

# smallgtexTransposed <- as.data.frame(t(smallgtex))
#add sex labels
# smallgtexLabelled <- merge(smallgtexTransposed, gtexPhenotypeNew, by=0)
# rownames(smallgtexLabelled)<- smallgtexLabelled$Row.names; smallgtexLabelled$Row.names <- NULL

#18738 - this is sex column
# colnames(smallgtexLabelled)[18786] <- "significant"
# colnames(smallgtexLabelled)[18786]

# smallgtexLabelled <- smallgtexLabelled[, colSums(smallgtexLabelled)>0]
# gtexColname<- colnames(smallgtexLabelled[,1:18737])
# gtexRowname<- rownames(smallgtexLabelled)
# dim(smallgtexLabelled)
# #log transf
# smallgtexLabelledL <- log2(smallgtexLabelled[,1:18737]+1)
# smallgtexLabelsForSex <- smallgtexLabelled$significant
# 
# # ggplot(sma)
# #normalise
# smallgtexLabelledN <- normalize.quantiles(as.matrix(smallgtexLabelledL))
# max(smallgtexLabelledN)
# min(smallgtexLabelledN)
# # gtexColname<- colnames(smallgtexLabelledL)
# # gtexRowname<- rownames(smallgtexLabelled)
# smallgtexLabelledN <- as.data.frame(smallgtexLabelledN)
# rownames(smallgtexLabelledN)<- gtexRowname
# colnames(smallgtexLabelledN) <- gtexColname
# ggplot(smallgtexLabelledN, aes( y=smallgtexLabelledN[,5]))+geom_boxplot()
# pca_for_gtex <- prcomp(smallgtexLabelledN, scale. = T)
# resultgtex <- pca_for_gtex$x
# dimension100 <- as.data.frame(resultgtex[,1:100])
# dimension100$significant <- as.factor(smallgtexLabelsForSex)
# dimension200 <- as.data.frame(resultgtex[,1:200])
# dimension200$significant <- as.factor(smallgtexLabelsForSex)
# # var_explained <- summary(pca_for_gtex)$importance[2,]*100
# withoutpcaDimension100 <- smallgtexLabelledN[,1:100]
# withoutpcaDimension100$significant <- as.factor(smallgtexLabelsForSex)
# # ggplot(resultgtex, aes(x=PC1, y=PC2))+
# #   geom_point(aes(colour = dimension100$significant))+
# #   xlab(paste("PC1", var_explained[1], "%")) +
# #   ylab(paste("PC2", var_explained[2], "%"))
# 
# 
# splitRes100gtex <- split_processingData_fctn(dimension100, 0.8)
# splitRes200gtex <- split_processingData_fctn(dimension200, 0.8)
# fullTemp <- fullGeneExpressforPCA[,1:100]
# fullTemp <- merge(fullTemp, labelsFullDf, by=0)
# fullTemp$Row.names<- NULL
# splitTemp <- split_processingData_fctn(fullTemp,.8)
# 
# plan(multisession, workers=availableCores())
# rfTemp <- gtexRFfctn(splitTemp, "PCA")
# rfgtex100<-gtexRFfctn(splitRes100gtex, "PCA")
# saveRDS(rfgtex100, "processed/GTEX_sex_MLRes100.rds")
# rfgtex200 <- gtexRFfctn(splitRes200gtex, "PCA")
# saveRDS(rfgtex200, "processed/GTEX_sex_MLRes200.rds")
# augBound <- list()
# augBound[[1]] <- rfgtex$aug
# augBound[[2]] <- rfgtex200$aug
# augDF <- bind_rows(augBound)
# augDF %>% group_by(dim) %>%
#   roc_curve(truth = significant, .pred_1) %>%
#   ggplot(aes(x=1-specificity, y=sensitivity, colour = as.factor(dim)))+
#   geom_path(size=1.1) +
#   geom_abline(slope = 1, intercept = 0, size=0.4, lty="dashed")+
#   theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill=NA),
#         aspect.ratio = 1)+
#   theme_bw(base_size = 16)
# 
# class(smallgtexLabelledN)
# str(smallgtexLabelledN)
# 
# 

#--------------
# upsampled <- gtexRFfctn(pcaProcessed[[4]], algorithm = "NMF")
# i don't think upsampling helps