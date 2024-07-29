######this all for 04_modelInterpretation script#####

library(tidyverse)
library(tidymodels)
library(correctR)
library(cowplot)
library(pheatmap)
library(ggsci)
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

nmfRfSVResList<- nmfRfSVResList[which(sapply(nmfRfSVResList, function(x) x$dim >5))]
#be mode diverse
#nmfRfSVResList <- nmfRfSVResList[which(sapply(nmfRfSVResList, '[[', 13 )>5)]
# read in xgboost file names
xgboostFileNames <- list.files("processed", full.names = T, pattern = "xgbResSV_")
length(xgboostFileNames)
xgbResList <- lapply(xgboostFileNames, readRDS)
xgbResList <- xgbResList[order(sapply(xgbResList, function(x) x$dim))]
xgbResList[[2]]<- NULL
length(xgbResList)
# so let's read in the PCA feature sets 

pcaXGBfileNames <- list.files("processed", full.names = T, pattern = "PCAxgbRes_")
length(pcaXGBfileNames)
pcaXGBResList <- lapply(pcaXGBfileNames, readRDS)
pcaXGBResList <- pcaXGBResList[order(sapply(xgbResList, function(x) x$dim))]
length(pcaXGBResList)

cvRFSV <- lapply(nmfRfSVResList, function(x){
  x$resDf
})
cvXGBsV <- lapply(xgbResList, function(x){
  x$resDf
})
cvAllrf <- do.call("rbind", cvRFSV) %>% mutate(model = "RF")
cvAllxgb <- do.call("rbind", cvXGBsV) %>% mutate( model = "XGB")
#cvRfandXGB <- bind_rows(cvAll, cvAllxgb)
nmfRFcvBP <- ggplot(cvAllrf, aes(x= as.factor(dim), y = mean, fill=model)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#4DBBD5B2"))+
  #scale_fill_npg()+
  #scale_fill_manual(values =c("#56B4E9", "tomato", "pink", "tan2", "grey", "turquoise4"))+
  labs(x="Number of NMF dimensions", y="Cross-validation AUC", fill= "Model") +
  theme_cowplot(font_size = 22)+
  facet_wrap(~"RF")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4), 
        strip.background = element_rect(color = "black",  linewidth = 0.4),
        legend.position = "none")
nmfRFcvBP
ggsave("output/nmfRFcvBP.png", nmfRFcvBP, width = 7, height = 5)
head(cvAllrf)
dim(cvAllrf)
cvRFFilt <- cvAllrf %>% filter(dim == "50" | dim == "500")
dim(cvRFFilt)
# so we know which AUC the final RF model from the cross validation was as well!
cvRfMetr <- lapply(nmfRfSVResList, function(x){
  r <- x$resDf
  maxAuc <- as.data.frame(r[which.max(r$mean),])
})
maxRFAucAllDims <- do.call("rbind", cvRfMetr) %>%
  select(mean, dim, lowerbound, upperbound, Algorithm)
# xgb
cvxgbMetr <- lapply(xgbResList, function(x){
  r <- x$resDf
  maxAuc <- as.data.frame(r[which.max(r$mean),])
})
maxXGBAucAllDims <- do.call("rbind", cvxgbMetr) %>%
  select(mean, dim, lowerbound, upperbound, Algorithm)

dim(cvAllxgb)
toPlotTwoFeatSetsCV <- bind_rows(cvAllxgb, cvRFFilt)

bestDimsRFvsXGBbplot <- ggplot(toPlotTwoFeatSetsCV, aes(x= as.factor(dim), y = mean, fill=model)) +
  geom_boxplot()+
  scale_fill_manual(values =c( "#4DBBD5B2","#CD202CB2")) + #, "tomato", "pink", "tan2", "grey", "turquoise4"))+
  labs(x="Number of NMF dimensions", y="Cross-validation AUC", fill= "Algorithm") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 30))+
  theme_cowplot(font_size = 20)+
  theme_cowplot(font_size = 18, line_size = 0.4)+
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA), 
        legend.position = "bottom")

#need to save above plot:
# ggsave("output/cccc", )
####how about PCA feature sets performance? we've only trained xgb for this...
# in the interest of time....
cvXGBsVpca <- lapply(pcaXGBResList, function(x){
  x$resDf
})
cvPcaAllxgb <- do.call("rbind", cvXGBsVpca) %>% mutate(model = "XGB")
cvPcaAllxgbFilt <- cvPcaAllxgb %>% filter(dim=="50"|dim=="500")

plotPCAXGBvsNMFxgbDf <- rbind(cvAllxgb, cvPcaAllxgbFilt)
pcaXGBcvBoxplot <- ggplot(cvPcaAllxgb, aes(x = as.factor(dim), y = mean, fill = model))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values= c("#8491B4B2"))+
  labs(x="Number of PCA dimensions", y="Cross-validation AUC")+#, fill= "Model") +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 30), legend.position = "none")+
  theme_cowplot(font_size = 22, line_size = 0.4)+
  #theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA))+
  facet_wrap(~"XGB")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4), 
        strip.background = element_rect(color = "black", size = 0.4))
pcaXGBcvBoxplot
ggsave("output/pcaXGBcvBP.png", pcaXGBcvBoxplot, width = 7, height = 5)

# instead let's just compare the NMF AND PCA XGB results:
pcaVSnmfXgbBoxplot <- ggplot(plotPCAXGBvsNMFxgbDf, aes(x = as.factor(dim), y = mean, fill = Algorithm))+
  geom_boxplot()+
  scale_fill_manual(values= c("#CD202CB2", "#8491B4B2"))+
  labs(x="Number of dimensions", y="Cross-validation AUC")+#, fill= "Model") +
  theme_cowplot(font_size = 22, line_size = 0.4)+
  #theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA))+
  facet_wrap(~"XGB")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4), 
        strip.background = element_rect(color = "black", size = 0.4), 
        legend.position = "bottom")

pcaVSnmfXgbBoxplot
ggsave("output/pcaVSnmfXgbBoxplot.png", pcaVSnmfXgbBoxplot ,width = 8, height = 6)
cvpcaXGBMetr <- lapply(pcaXGBResList, function(x){
  r <- x$resDf
  maxAuc <- as.data.frame(r[which.max(r$mean),])
})
maxPCAaucAllDims <- do.call("rbind", cvpcaXGBMetr)%>%
  select(mean, lowerbound, upperbound, dim, Algorithm)
# ok so we can see that the PCA feature sets when an xgb model is trained
# on them they give poorer performance than the NMF feature sets

#############--
predictionsAllRf <- lapply(nmfRfSVResList, function(x){
  x$aug
})

# emp <- lapply(nmfRfSVResList, function(x){
#   list("pred" = x$aug, "AUC"= x$AUC$.estimate)
# })

predictionsXGB <- lapply(xgbResList, function(x){
  x$aug
})

rfAllPredictions <- bind_rows(predictionsAllRf)
xgbAllPred <- bind_rows(predictionsXGB)

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

df_text <- data.frame(dim = listDims, x = 0.2, y = 0.8,
                      label = paste0("AUC = ",round(aucDf$auc_scores, 3)))
head(df_text)
# now we can plot the predictions
rfMetrics <- rfAllPredictions %>% 
  group_by(dim) %>%
  roc_curve(truth = significant, .pred_FALSE)


noNetworkROCcurveRF <-ggplot(rfMetrics, aes(x=1-specificity, y=sensitivity, colour=as.factor(dim)))+
  geom_line(linewidth=1.5, show.legend = F)+ 
  geom_abline(slope = 1, intercept = 0, linewidth=0.4, lty="dashed", alpha = 0.5)+
  
  theme(panel.border = element_rect(colour = "black", linewidth = 1.0, fill="white"),
         aspect.ratio = 1, legend.position="none")+
  theme_bw(base_size = 20)+# theme(legend.position = "none")+
  scale_colour_npg()+
  facet_wrap(~ paste0(dim, " NMF dimensions"))+ 
  theme(strip.text = element_text(size=20))+
  geom_text(data = df_text, mapping = aes(x = 0.3,  y=0.85, label = label),
            size= 8, show.legend = F)# +

noNetworkROCcurveRF
ggsave("output/nonetworkROCcurveRF.png", noNetworkROCcurveRF, width = 12, height = 8)


rocCurveXGB <- xgbAllPred %>% group_by(dim) %>%
  roc_curve(truth = significant, .pred_FALSE) %>%
  ggplot(aes(x=1-specificity, y=sensitivity, colour=as.factor(dim)))+
  geom_path(linewidth=1.5)+
  geom_abline(slope = 1, intercept = 0, size=0.4, lty="dashed", alpha = 0.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill="white"),
        aspect.ratio = 1)+
  theme_bw(base_size = 26) +
  guides(colour = guide_legend(title = "No. of NMF dimensions"))+
  facet_wrap(~"XGB")+
  theme(strip.text = element_text(size=24), legend.position = "bottom")+
  scale_colour_manual(values = c("#E64B35B2", "#3C5488B2"))

rocCurveXGB
ggsave("output/nonetworkROCcurveXGBnmf.png", rocCurveXGB, width = 7.5, height = 8)


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
# an also for that model, what was the area under the curve??/??/
bestModMetrics <- nmfRfSVResList[[posBestFit]]$resDf
# finalModMetrics <- which.max(bestModMetrics$mean)
finalModMetrics <- bestModMetrics[which.max(bestModMetrics$mean), ]

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
dim(top5FeatsModeldf)
labelsFullDf <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
str(labelsFullDf$significant)
labelsFullDf$significant<- as.factor(labelsFullDf$significant)
all(rownames(top5FeatsModeldf)%in% rownames(labelsFullDf))

all(rownames(top5FeatsModeldf)== rownames(labelsFullDf))
#they're in order so can just do cbind
top5FeatsModeldfLabelled <- cbind(top5FeatsModeldf, labelsFullDf)
# q1 <- list()
for (c in 1:(ncol(top5FeatsModeldfLabelled)-1)) {
  #print(c)
  feat <- colnames(top5FeatsModeldfLabelled)[c]
  name <- sprintf("%s_nmf_%s",c,colnames(top5FeatsModeldfLabelled)[c])
  pl <- ggplot(top5FeatsModeldfLabelled)+ 
         aes(y = top5FeatsModeldfLabelled[,c], x= significant, fill=significant)+
    geom_boxplot()+
    xlab("") + ylab(paste0("NMF " ,colnames(top5FeatsModeldfLabelled)[c]))+
    theme_cowplot(font_size = 16)+
    theme(panel.background = element_rect(colour = "black", fill = NA, linewidth = 0.4))+
    # scale_fill_bmj()+
    scale_fill_manual(values = c("#69BE28B2", "#E37222B2"))+ 
    scale_x_discrete(labels = c("Associated", "Not associated"))+
    theme(legend.position = "none")
  ggsave(paste0("output/", name,".png"), pl, height = 3, width = 4)
  print(name)
  print(pl)
  # q1[[feat]]<- pl
}
# l <- marrangeGrob(q1, ncol = 2,nrow = 4)
# ggsave("output/test.pdf", l, height = 20, width = 20)

# top5FeatsModeldf<- top5FeatsModeldf %>% 
#   dplyr::mutate(TotalWeight = rowSums(top5FeatsModeldf))
# x <- as.data.frame(top5FeatsModeldf$sum)

#maybe put the top important feature in order???
top5FeatsModeldf <- top5FeatsModeldf%>% arrange(desc(Feature272))
humanGenesSymbs <- read.table("processed/human_coding_genes.txt", sep = "\t", header = T)
humanGenesSymbs <- humanGenesSymbs%>% filter(hgnc_symbol!=""); rownames(humanGenesSymbs) <- humanGenesSymbs$ensembl_gene_id
humanGenesSymbs$ensembl_gene_id <-NULL
head(humanGenesSymbs) 
# so now gonna merge with top feats DF
head(top5FeatsModeldfLabelled)
topFeatsWithSymbols <- merge(top5FeatsModeldf, humanGenesSymbs, by=0); rownames(topFeatsWithSymbols) <- topFeatsWithSymbols$hgnc_symbol; topFeatsWithSymbols$Row.names<-NULL; topFeatsWithSymbols$hgnc_symbol <-NULL
head(topFeatsWithSymbols)
#topFeatsWithSymbols <- topFeatsWithSymbols %>% arrange(desc(Feature272))
mat <- as.matrix(topFeatsWithSymbols)
head(mat)
mat <- apply(mat, 2, rank)
head(mat)
#pdf("processed/testFig.pdf", width = 10, height = 10)
# just select the top 15 genes in the most important feature
heatmap <- pheatmap::pheatmap(mat[1:20,], border_color = "white",
                   cluster_rows = F, 
                   cluster_cols = F, show_rownames = T, 
                   )
ggsave("processed/heatmap.pdf",heatmap, height = 5, width = 10)

####clusterprofiler code####

getEachFeature_fctn <- function(dataframeOfFeats){
  single <- list()
  for (j in 1:ncol(dataframeOfFeats)){ 
    print(j)
    single[[j]] <- data.frame(dataframeOfFeats[,j],
                            row.names = rownames(dataframeOfFeats))
    colnames(single[[j]]) <- colnames(dataframeOfFeats)[j]
    single[[j]] <- single[[j]] %>% arrange(desc(colnames(single[[j]])))
    #single[[j]]$Gene <- rownames(dataframeOfFeats)
  }
  return(single)
} 

eachFeat <- getEachFeature_fctn(top5FeatsModeldf)
saveRDS(eachFeat, "processed/top5Feats.rds") # do the gsea on my laptop cause clusterprofiler 
#not installing




# ggplot(modelFeatsDf, aes(x = , y= ))+geom_
colnames <- paste0("Feature", 1:ncol(bestDf))
unlabelledGenes <- readRDS("processed/unStudiedGenes.rds")
head(unlabelledGenes[1:4,1:5])
unlabelledGenesPred <- augment(bestFit, unlabelledGenes)
# so because we are returned a tibble, when converting to df the rownames disappear
# temp <- as.data.frame(unlabelledGenesPred); rownames(temp) <- rownames(unlabelledGenesPred)
head(unlabelledGenesPred[,2312:2313])

unlabelledGenesPreddf <- as.data.frame(unlabelledGenesPred); rownames(unlabelledGenesPreddf) <- rownames(unlabelledGenesPred)
head(unlabelledGenesPreddf[1:3,1:4])
length(which(unlabelledGenesPreddf$.pred_class=="TRUE"))
length(which(unlabelledGenesPreddf$.pred_class!= "TRUE"))





##external validation part 3 
#- MGI
head(unlabelledGenesPreddf[1:4,1:5])
mgiGenes <- read.table("processed/annotatedMGIgenes.txt", header = T, sep = "\t")
#just remove duplicates since they're all the same level anyways
mgiGenes <- mgiGenes %>% distinct(ensembl_gene_id, .keep_all = T)
str(mgiGenes$ensembl_gene_id)
# make gene names rownames
rownames(mgiGenes)<- mgiGenes$ensembl_gene_id ; mgiGenes$ensembl_gene_id<-NULL
head(mgiGenes)
mgiGenes$significant<- as.factor(mgiGenes$significant)
head(unlabelledGenes[1:3,1:5])
# let's just get the genes in both datasets

head(unlabelledGenesPreddf[1:5,1:5])
length(which(rownames(unlabelledGenesPreddf) %in% rownames(mgiGenes)))
# so 1001 of the mgi genes are known to be associated to a skeletal phenotype
nrow(unlabelledGenesPreddf) - (length(which(rownames(unlabelledGenesPreddf) %in% rownames(mgiGenes))))
# and 8000 genes are unstudied
#which are those genes then? that aren't in the mgi genes (i.e., unstudied)
totallyUnstudied <- rownames(unlabelledGenesPreddf)[!rownames(unlabelledGenesPreddf) %in% rownames(mgiGenes)]
class(totallyUnstudied)
length(totallyUnstudied)
# so those genes that aren't in the mgi DB, assuming they're not associated, let's label them as that
unstudied <- data.frame(significant = rep(c("FALSE"), times = length(totallyUnstudied)), row.names = totallyUnstudied)
head(unstudied) # nice
tail(unstudied)
str(unstudied$significant)
unstudied$significant <- as.factor(unstudied$significant)
nrow(unstudied)
nrow(mgiGenes) # cool
mgiGenesPlusUnstudied <- rbind(mgiGenes, unstudied)
# now let's retain only those genes in mgiGenesplus unstudied DF that are also in our unstudied/unlabelled genes
unlabelledGenesWithMGIannot <- merge(unlabelledGenesPreddf,mgiGenesPlusUnstudied, by=0); rownames(unlabelledGenesWithMGIannot)<- unlabelledGenesWithMGIannot$Row.names; unlabelledGenesWithMGIannot$Row.names<- NULL
# cool
head(unlabelledGenesWithMGIannot[1:3,1:3])

mgiCurve <-roc_curve(unlabelledGenesWithMGIannot, truth = significant, .pred_TRUE,
              event_level = "first")
head(mgiCurve)
mgiCurve <- mgiCurve %>% mutate(database = "MGI")
head(mgiCurve)
rocaucMGI <- roc_auc(unlabelledGenesWithMGIannot, significant, .pred_TRUE)
rocaucMGI <- rocaucMGI %>% mutate(database = "MGI")


# now looking at human phenotype ontology- genes that are annotated to a skeletal abnormality
# lets quickly do this 
hpOnt <- read.table("processed/humphenotOntGenes.txt", header = T, sep = "\t")
head(hpOnt)
length(unique(hpOnt$ensembl_gene_id))

rownames(hpOnt)<- hpOnt$ensembl_gene_id; hpOnt$ensembl_gene_id<-NULL
hpOnt$significant <- as.factor(hpOnt$significant)
length(which(rownames(unlabelledGenesPreddf) %in% rownames(hpOnt)))
# we have 1415 genes that are unlabelled (unstudied in impc) to be known to have skel phenotype according to hpo
length(rownames(unlabelledGenesPreddf)[!rownames(unlabelledGenesPreddf) %in% rownames(hpOnt)])
# the remainder of these 7586genes^ are not known to have a skeletal abnormality in the human phenotype ont
unstGene <- rownames(unlabelledGenesPreddf)[!rownames(unlabelledGenesPreddf) %in% rownames(hpOnt)]
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

unlabelledGenesWithHPOannot <- merge(unlabelledGenesPreddf, hpoGenesPlusUnstudied, by= 0); rownames(unlabelledGenesWithHPOannot) <- unlabelledGenesWithHPOannot$Row.names; unlabelledGenesWithHPOannot$Row.names <- NULL
head(unlabelledGenesWithHPOannot[1:3,1:5])
unlabelledGenesWithHPOannot$significant <-as.factor(unlabelledGenesWithHPOannot$significant)
#cool looks good i guess
hpoCurve <- roc_curve(unlabelledGenesWithHPOannot, truth = significant, .pred_TRUE)
hpoCurve <- hpoCurve %>% mutate(database = "HPO")

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
ggsave("output/ROCmgihpo.png", mgihpoROC, width = 8, height = 4.5)  
 # theme_cowplot()
  #theme_minimal_grid(font_size = 17)
  #theme_bw(base_size = 20)
#w what is the area under curve???????????
rocaucHPO <- roc_auc(unlabelledGenesWithHPOannot, significant, .pred_TRUE)
rocaucHPO <- rocaucHPO %>% mutate(database = "HPO")
rocaucHPO

## part 4? what are the top rankeed genes?

# let's see

orderedGenes <- unlabelledGenesPreddf[order(unlabelledGenesPreddf$.pred_TRUE, decreasing = T),]
top20genes <- orderedGenes[1:20,]
top20genesWsymbs <- merge(top20genes, humanGenesSymbs, by=0)
rownames(top20genesWsymbs) <- top20genesWsymbs$hgnc_symbol; top20genesWsymbs$Row.names <-NULL; top20genesWsymbs$hgnc_symbol<-NULL
top20genesWsymbs <- top20genesWsymbs[order(top20genesWsymbs$.pred_TRUE, decreasing = T),]
## can save it if you like...
write.table(as.data.frame(top20genesWsymbs[,1:3]), "processed/temporary.txt", col.names = T, sep = "\t", quote = F)
