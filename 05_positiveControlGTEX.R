library(preprocessCore)
library(future)
library(future)

gtex <- read.delim("raw/rnaSeqGtex/gtex_Analysis.gct", header = T,skip = 2, sep = "\t")

#just retain protein coding genes
# temp <- gtex[1:5,1:3]
# 
# 
# temp4 <- as.data.frame(do.call('rbind', strsplit(as.character(gtex$Name),
#                                    '.', fixed = T)))

#remove first 2 columns now
gtex <- gtex[,-c(1:2)]

length(unique(temp4$V1))

humanProteinCoding <- read.table("processed/human_coding_genes.txt", header = T,
                                 stringsAsFactors = F, sep = "\t")
#

row.names(gtex) <- temp4$V1


gtex$ensembl_gene_id <- rownames(gtex)
gtexProteinCoding <- merge(humanProteinCoding, gtex, by = "ensembl_gene_id")
row.names(gtexProteinCoding)<- gtexProteinCoding$ensembl_gene_id; gtexProteinCoding$ensembl_gene_id <- NULL


#now we have all the gene expression data in ensembl ID format
filtgtex <- gtexProteinCoding[rowSums(gtexProteinCoding)>0,]
rowsFiltgtex <- rownames(filtgtex)
head(filtgtex[1:5,1:5])
gtexNor <- normalize.quantiles(as.matrix(filtgtex))
rownames(gtexNor)<- rowsFiltgtex
gtexNor<- as.data.frame(gtexNor)


#perform PCA
pca_gtex <- prcomp(gtexProteinCoding, scale. = T)
pcaGtexpcs <- pca_gtex$x
for (i in 1:length(ranks)) {
  #create empty lists
  gtexpcList <- list()
  #loop through each of the dimensions
  for (l in (ranks)) {
    #grab each dimension
    dime <- pcsLabelled[,1:l]
    gtexpcList <- append(gtexpcList, list(dime))
  }
}
#add labels function
getLabelledGenesFctn <- function(matrixList, knownLabels){
  merged <- merge(matrixList, knownLabels, by=0);rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  return(merged)
}  

#apply to result
labelledgtex <- map(.x = gtexpcList, .f = getLabelledGenesFctn, knownLabels = labelsFullDf)

performTtestFctn = function(listRes, labels, algorithm){
  if(algorithm == "PCA"){
    pvalList <- numeric()
    for (c in 1:ncol(listRes)){
      pvalueResult <- t.test(y=as.logical(labels), x=(listRes[,c]))$p.value
      pvalList <-c(pvalList, pvalueResult)
    }
    minPval <- min(pvalList)
    mindim <- which.min(pvalList)
    dataframe <- data.frame(pVals = minPval,
                            Dimension = ncol(listRes),
                            Feature=mindim, 
                            Algorithm = algorithm)
    
  }
  else if(algorithm == "NMF"){
    
    pvalList <- numeric()
    for (c in 1:ncol(listRes)){
      pvalueResult <- t.test(y=as.logical(labels), x=(listRes[,c]))$p.value
      pvalList <-c(pvalList, pvalueResult)
    }
    minPval <- min(pvalList)
    mindim <- which.min(pvalList)
    dataframe <- data.frame(pVals = minPval,
                            Dimension = ncol(listRes),
                            Feature=mindim, 
                            Algorithm = algorithm)
    
  }
}
gtextTest <- map(.x = gtexpcList,.f = performTtestFctn, labels=labelsFullDf$significant, algorithm = "PCA")

gtex_5dim <- labelledgtex[[5]]
#temporary RF function not using too much tuning or repeats of cross-validation just for testing
justToTestRF <- function(preprocessResult, algorithm){
  if(algorithm == "NMF"){
    model_RF <- rand_forest(trees = 500, mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=3)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      #trees(range = c(1,2000)),
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
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("workflow" = wkflow, "res" = res, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    saveRDS(result, sprintf("%sMLRes_%s.rds",algorithm, (ncol(preprocessResult$train)-1)))
    
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
      # mtry(range = c(1,limit)),
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
      fit(data = preprocessResult$train) 
    
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
    saveRDS(result, sprintf("processed/%sMLRes_%s.rds",algorithm, ncol(preprocessResult$train)-1))
    
    return(list("workflow" = wkflow, "res" = res, "resDF"=resDf, "finalMod" = final_model, "tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
  }
  
}
#
plan(multisession, workers = 20)
gtexTemporaryList <- list()
gtexPrePR <- split_processingData_fctn(data = gtex_5dim,proportion = 0.8)
# testFctn <- split_processingData_fctn(data = gtex_5dim, proportion = .8)

gtexTemporaryList[[1]] <- gtexPrePR
#test whether can detect skeletal phenotype
gtexMLRES <- map(.x = gtexTemporaryList, .f = justToTestRF, algorithm="PCA")
#looks like it doesn't



sampleNames <- (colnames(gtexProteinCoding))
names <- data.frame(id = sampleNames)
colnames(sampleNames)[1]<- id
# df <- (do.call("rbind", strsplit(as.character(sampleNames$colnames(gtexProteinCoding), '.', fixed=T))))

# just get the first part of ensemblID
id <- names %>% mutate(first = unlist(lapply(strsplit(names$id, '\\.'), function(x)x[1])))#
id <- id %>% mutate(second = unlist(lapply(strsplit(names$id, '\\.',), function(x)x[2]))) 
#make it as the same format as the sex is
id$unique <- paste(id$first, id$second, sep = '-')
dim(id)
colnames(gtexProteinCoding)
#now make these the columnames
colnames(filtgtex) <- id$unique

gtexPhenotype <- read.table("raw/GTEX_v7_SubjectPhenotypeDS.txt", header = T)
length(unique(gtexPhenotype$SUBJID))
rownames(gtexPhenotype) <- gtexPhenotype$SUBJID
gtexPhenotypeNew <- data.frame(gtexPhenotype$SEX)
rownames(gtexPhenotypeNew)<- gtexPhenotype$SUBJID
#gtexTranspose <- as.data.frame(t(filtgtex))                       

#just get a small subset of the gtex data
smallgtex <- filtgtex[,1:790]

smallgtexTransposed <- as.data.frame(t(smallgtex))
#add sex labels
smallgtexLabelled <- merge(smallgtexTransposed, gtexPhenotypeNew, by=0)
rownames(smallgtexLabelled)<- smallgtexLabelled$Row.names; smallgtexLabelled$Row.names <- NULL

#18738 - this is sex column
colnames(smallgtexLabelled)[18786] <- "significant"
colnames(smallgtexLabelled)[18786]

smallgtexLabelled <- smallgtexLabelled[, colSums(smallgtexLabelled)>0]
gtexColname<- colnames(smallgtexLabelled[,1:18737])
gtexRowname<- rownames(smallgtexLabelled)
dim(smallgtexLabelled)
#log transf
smallgtexLabelledL <- log2(smallgtexLabelled[,1:18737]+1)
smallgtexLabelsForSex <- smallgtexLabelled$significant

# ggplot(sma)
#normalise
smallgtexLabelledN <- normalize.quantiles(as.matrix(smallgtexLabelledL))
max(smallgtexLabelledN)
min(smallgtexLabelledN)
# gtexColname<- colnames(smallgtexLabelledL)
# gtexRowname<- rownames(smallgtexLabelled)
smallgtexLabelledN <- as.data.frame(smallgtexLabelledN)
rownames(smallgtexLabelledN)<- gtexRowname
colnames(smallgtexLabelledN) <- gtexColname
ggplot(smallgtexLabelledN, aes( y=smallgtexLabelledN[,5]))+geom_boxplot()
pca_for_gtex <- prcomp(smallgtexLabelledN, scale. = T)
resultgtex <- pca_for_gtex$x
dimension100 <- as.data.frame(resultgtex[,1:100])
dimension100$significant <- as.factor(smallgtexLabelsForSex)
dimension200 <- as.data.frame(resultgtex[,1:200])
dimension200$significant <- as.factor(smallgtexLabelsForSex)
# var_explained <- summary(pca_for_gtex)$importance[2,]*100
withoutpcaDimension100 <- smallgtexLabelledN[,1:100]
withoutpcaDimension100$significant <- as.factor(smallgtexLabelsForSex)
# ggplot(resultgtex, aes(x=PC1, y=PC2))+
#   geom_point(aes(colour = dimension100$significant))+
#   xlab(paste("PC1", var_explained[1], "%")) +
#   ylab(paste("PC2", var_explained[2], "%"))


splitRes100gtex <- split_processingData_fctn(dimension100, 0.8)
splitRes200gtex <- split_processingData_fctn(dimension200, 0.8)
fullTemp <- fullGeneExpressforPCA[,1:100]
fullTemp <- merge(fullTemp, labelsFullDf, by=0)
fullTemp$Row.names<- NULL
splitTemp <- split_processingData_fctn(fullTemp,.8)

plan(multisession, workers=availableCores())
rfTemp <- justToTestRF(splitTemp, "PCA")
rfgtex100<-justToTestRF(splitRes100gtex, "PCA")
rfgtex200 <- justToTestRF(splitRes200gtex, "PCA")
augBound <- list()
augBound[[1]] <- rfgtex$aug
augBound[[2]] <- rfgtex200$aug
augDF <- bind_rows(augBound)
augDF %>% group_by(dim) %>%
  roc_curve(truth = significant, .pred_1) %>%
  ggplot(aes(x=1-specificity, y=sensitivity, colour = as.factor(dim)))+
  geom_path(size=1.1) +
  geom_abline(slope = 1, intercept = 0, size=0.4, lty="dashed")+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill=NA),
        aspect.ratio = 1)+
  theme_bw(base_size = 16)

class(smallgtexLabelledN)
str(smallgtexLabelledN)



#--------------
# upsampled <- justToTestRF(pcaProcessed[[4]], algorithm = "NMF")
# i don't think upsampling helps