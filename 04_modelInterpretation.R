######this all for 04_modelInterpretation script#####

library(tidyverse)
library(correctR)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# need to read in results from last script
# nmfRes <- readRDS("processed/tempMLResults.rds")
# read in also the list of different dimensions:
# nmfDataframes <- readRDS("processed/nmfDataframes.rds")


#read in random forest file names
randForestfileNames <- list.files("processed", full.names = T, pattern = "RFResSV_")
randForestfileNames

nmfRfSVResList <- lapply(randForestfileNames, readRDS)
nmfRfSVResList[[1]]$dim
nmfRfSVResList[[2]]$dim
nmfRfSVResList <- nmfRfSVResList[order(sapply(nmfRfSVResList, function(x) x$dim))]
nmfRfSVResList[[2]]$dim
#remove the 5 dimension one
nmfRfSVResList[[1]] <- NULL
# read in xgboost file names
xgboostFileNames <- list.files("processed", full.names = T, pattern = "xgbResSV")
xgbResList <- lapply(xgboostFileNames, readRDS)



cvRFSV <- lapply(nmfRfSVResList, function(x){
  x$resDf
})
cvXGBsV <- lapply(xgbResList, function(x){
  x$resDf
})
cvAll <- do.call("rbind", cvRFSV) %>% mutate(model = "RF")
cvAllxgb <- do.call("rbind", cvXGBsV) %>% mutate( model = "XGB")
cvRfandXGB <- bind_rows(cvAll, cvAllxgb)
ggplot(cvRfandXGB, aes(x= as.factor(dim), y = mean, fill=model)) +
  geom_boxplot()+
  scale_fill_manual(values =c("#56B4E9", "tomato", "pink", "tan2", "grey", "turquoise4"))+
  labs(x="", y="AUC", fill= "Model") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 30))+
  theme_classic(base_size = 18)+
  scale_x_discrete(labels = c("10 dimensions", "50 dimensions", "100 dimensions",
                              " 150 dimensions", "200 dimensions", "500 dimensions"),
                   guide = guide_axis(angle = 15))


predictionsAll <- lapply(nmfRfSVResList, function(x){
  x$aug
})

predictionsXGB <- lapply(xgbResList, function(x){
  x$aug
})

rfAllPredictions <- bind_rows(predictionsAll)
xgbAllPred <- bind_rows(predictionsXGB)
everything <- bind_rows(rfAllPredictions, xgbAllPred)
rocCurve <- function(predictions){
  predictions %>% group_by(dim) %>%
  roc_curve(truth = significant, .pred_FALSE) %>%
  ggplot(aes(x=1-specificity, y=sensitivity, colour=as.factor(dim)))+
  geom_path(linewidth=0.9)+
  geom_abline(slope = 1, intercept = 0, size=0.4, lty="dashed", alpha = 0.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill="white"),
        aspect.ratio = 1)+
  theme_bw(base_size = 20) +
  guides(colour = guide_legend(title = "Dimension"))+
  scale_colour_manual(values = c("pink3", "tomato", "tan2", "purple4", "turquoise3", "darkgreen"))
  } 
  
noNetworkROCcurveRF <- rocCurve(predictions = rfAllPredictions)
noNetworkROCcurveRF +
  facet_wrap(~"RF")
ggsave("processed/nonetworkROCcurveRF.pdf", noNetworkROCcurveRF, width = 10, height = 8)
rocCurveXGB <- xgbAllPred %>% group_by(dim) %>%
  roc_curve(truth = significant, .pred_FALSE) %>%
  ggplot(aes(x=1-specificity, y=sensitivity, colour=as.factor(dim)))+
  geom_path(linewidth=0.9)+
  geom_abline(slope = 1, intercept = 0, size=0.4, lty="dashed", alpha = 0.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill="white"),
        aspect.ratio = 1)+
  theme_bw(base_size = 20) +
  guides(colour = guide_legend(title = "Dimension"))+
  scale_colour_manual(values = c("turquoise3", "darkgreen"))+
  facet_wrap(~"XGB")
  
rocCurveXGB
#### get all the auc scores ####
aucRF <-lapply(nmfRfSVResList, function(x){ # where the input is a list containing multiple ML results
  x$AUC$.estimate
})
df_auc <- lapply(aucRF, function(x){
  df<- data.frame(auc_scores = x)
})
head(df_auc)
# here just binding all the scores for each model

aucDf <- do.call("rbind", df_auc)
head(aucDf)
dims <- lapply(nmfRfSVResList, function(x){
  x$dim
})
listDims <- unlist(dims)
listDims
dimension <- data.frame(dimensions = listDims, Algorithm = "NMF")
resDf <- cbind(aucDf, dimension)
head(resDf)

#finding which dimension the model that gave the highest auc score was trained on
posBestFit <- which.max(resDf$auc)
posBestFit
bestDimension <- resDf[which.max(resDf$auc),]$dimensions
bestDimension
datausedForML <- readRDS("processed/nmfDataframes.rds")
datausedForML[[1]] <-NULL # need to get rid of 5 dimension since not using anymore
for (d in 1:length(nmfRfSVResList)){
  #print(nmfRes[[d]])
  if(ncol(datausedForML[[d]])-1 == bestDimension){#allnmfmatrix is a list containing all the nmf reduced matrices +labels
    # print(allNmfMatrix[[d]])
    bestDf <- as.data.frame(datausedForML[[d]])
  }
}
dim(bestDf)
#find the final fit for that model
bestFit <- nmfRfSVResList[[posBestFit]]$finalFit
rfCVdf_forCR <- nmfRfSVResList[[posBestFit]]$dfForCorrectR

#need to also get xgb mod 500 dim
for(e in 1:length(xgbResList)){
  if (xgbResList[[e]]$dim == bestDimension){
    dfxgb <- xgbResList[[e]]$dfForCorrectR
  }
}
dfCR <-rbind(rfCVdf_forCR, dfxgb)
testRes <- repkfold_ttest(dfCR, n1=80, n2=20, k = 5, r = 3)

#lets look at the variable importance
bestModVarImport <- nmfRfSVResList[[posBestFit]]$importanceDf
nmfRfSVResList[[posBestFit]]$importancePlot
head(bestModVarImport)

bestModVarImport <- bestModVarImport %>% 
  mutate(sign = case_when(Importance<0 ~"negative", TRUE~"positive"))
head(bestModVarImport)
dim(bestModVarImport)
top5Feats <- bestModVarImport[1:5,]
fiveVars <- top5Feats$Variable

#subset for the ones that are positive sign
# impFeats <- bestModVarImport %>% filter(sign=="positive")
# head(impFeats)
# dim(impFeats)
# #then extract the variables that contribute to model's predictions
# modelFeats <- impFeats$Variable
# #find those features 
# modelFeatsDf <- bestDf %>% dplyr::select(all_of(modelFeats))# these are gonna be input for GSEA 
# colnames(modelFeatsDf) <- sub('V', 'Feature', colnames(modelFeatsDf))
# head(modelFeatsDf)
top5FeatsModeldf <- bestDf %>% dplyr::select(all_of(fiveVars))
head(top5FeatsModeldf)
colnames(top5FeatsModeldf)<- sub("V", "Feature", colnames(top5FeatsModeldf))
head(top5FeatsModeldf)

top5FeatsModeldf<- top5FeatsModeldf %>% 
  dplyr::mutate(TotalWeight = rowSums(top5FeatsModeldf))
# x <- as.data.frame(top5FeatsModeldf$sum)

#maybe put the top important feature in order???
top5FeatsModeldf <- top5FeatsModeldf%>% arrange(desc(Feature272))
# 
# #let's just focus on one feature:
# firstFeature <- data.frame(top5FeatsModeldf[,1], row.names = rownames(top5FeatsModeldf))
# # nameFeat <- colnames(top5FeatsModeldf)[1]
# colnames(firstFeature)[1] <- colnames(top5FeatsModeldf)[1]
# firstFeat <- firstFeature %>% dplyr::arrange(desc(Feature272)) 
# ##
library(gplots)
?heatmap.2()
# heatmap(mat)
mat <- as.matrix(top5FeatsModeldf)
library(future)
library(gplots)

plan(multisession, workers = availableCores())
library(pheatmap)
#pdf("processed/testFig.pdf", width = 10, height = 10)
# just select the top 15 genes in the most important feature
heatmap <- pheatmap::pheatmap(mat[1:15,], border_color = "white",
                   cluster_rows = F, 
                   cluster_cols = F, 
                   )
ggsave("processed/heatmap.pdf",heatmap, height = 5, width = 10)
# hm <- heatmap.2(x = mat[1:15,], 
#           col = RColorBrewer::brewer.pal(9, c("RdBu")), 
#           #col="bluered",
#           dendrogram = "none", 
#           Rowv = F, 
#           Colv = F,
#           tracecol = NA,
#           rowsep = 1:nrow(mat),
#           #colsep = 1:ncol(mat)-1,
#           trace = 'none')
# hm

# ggplot(modelFeatsDf, aes(x = , y= ))+geom_
colnames <- paste0("Feature", 1:ncol(bestDf))
unlabelledGenes <- readRDS("processed/unStudiedGenes.rds")
unlabelledGenes$prediction <- predict(bestFit, unlabelledGenes)
head(unlabelledGenes[,2312:2313])
length(which(is.na(unlabelledGenes)))
length(which(unlabelledGenes$prediction=="TRUE"))
length(which(unlabelledGenes$prediction!= "TRUE"))


# vip(bestModVarImport, geom = 'point', mapping=aes(colour = sign))
##


####clusterprofiler code####
# enrichKEGG(gene = row.names(bestDf))

getEachFeature_fctn <- function(dataframeOfFeats){
  sing <- list()
  for (j in 1:(ncol(dataframeOfFeats)-1)){ # minus 1 cause don't want last col
    print(j)
    sing[[j]] <- data.frame(dataframeOfFeats[,j],
                            row.names = rownames(dataframeOfFeats))
    colnames(sing[[j]]) <- colnames(dataframeOfFeats)[j]
    sing[[j]] <- sing[[j]] %>% arrange(desc(colnames(sing[[j]])))
    #sing[[j]]$Gene <- rownames(dataframeOfFeats)
  }
  return(sing)
} 
# every time I successfully do a for loop I feel 100 times smarter
eachFeat <- getEachFeature_fctn(top5FeatsModeldf)
saveRDS(eachFeat, "processed/top5Feats.rds") # do the gsea on my laptop cause clusterprofiler 
#not installing





##external validation part 3
humanPhenotype <- read.delim("genes_for_HP_0000924", col.names = c("geneID", "geneSymbol"))
externalDatabase <- read.delim("https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt", header = F, 
                               sep = "\t")
externalDatabase$V6 <- NULL
colnames(externalDatabase)<- c("humanMarkerSymbol", "humanEntrezGeneID", "mouseMarkerSymbol", "mgiMarkerID", "mammalianPhenotypeID")
head(externalDatabase)
#search for skeletal phenotype and 
patterns <- c("MP:0005390", "MP:0005371")

exDb <- externalDatabase %>% filter(grepl(paste(patterns, collapse = '|'), mammalianPhenotypeID))
