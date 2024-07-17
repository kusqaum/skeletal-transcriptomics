######this all for 04_modelInterpretation script#####

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
# need to read in results from last script
# nmfRes <- readRDS("processed/tempMLResults.rds")
# read in also the list of different dimensions:
# nmfDataframes <- readRDS("processed/nmfDataframes.rds")


#read in random forest file names
randForestfileNames <- list.files("processed", full.names = T, pattern = "RFResSV")
randForestfileNames
rfResList <- lapply(randForestfileNames, readRDS)

# read in xgboost file names
xgboostFileNames <- list.files("processed", full.names = T, pattern = "xgbResSV")
xgbResList <- lapply(xgboostFileNames, readRDS)

#### get all the auc scores ####
nmfDataMlResList <- list()
nmfDataMlResList[[1]]<- NMFMLRes_5
nmfDataMlResList[[2]]<- NMFMLRes_10
auc <-lapply(nmfDataMlResList, function(x){ # where the input is a list containing multiple ML results
  x$AUC$.estimate
})

df_auc <- lapply(auc, function(x){
  df<- data.frame(auc_scores = x)
})
View(df_auc)
# here just binding all the scores for each model
#so instead of getting just the best AUC metric.. get the metrics for all models trained: and then plot boxplots of them
t <- nmf50ml$res$.metrics
auc_50 <- do.call("rbind", t)
auc_50 <- auc_50 %>% mutate(dim = 50, Algorithm = "NMF")
ggplot(auc_50, aes(x=as.factor(dim),y=.estimate, col = Algorithm))+
  geom_boxplot()


aucDf <- do.call("rbind", df_auc)
head(aucDf)
dimension <- data.frame(dimensions = c(5,10), Algorithm = "NMF")
resDf <- cbind(aucDf, dimension)
head(resDf)
pcaAucDf <- data.frame(auc_scores = c(0.5, 0.509), dimensions = c(5,10), Algorithm = 'PCA')
bothAUC <- rbind(resDf, pcaAucDf)
head(bothAUC)
#then plotting - need to change this because usually for aucs, should plot all of it as a boxplot to show spread
ggplot(bothAUC, aes(x=as.factor(dimensions), y=auc_scores, col=Algorithm))+
  geom_point(size=4.5) +
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill=NA),
        panel.background = element_blank(),
        #strip.text = element_text(),
        #legend.position = "none",
        #legend.text = element_text("dimensionality"),
        #legend.title = element_text("dimension"),
        aspect.ratio = 1)+
  scale_colour_manual(values=c("tan2", "tomato"))+
  xlab(expression(italic(" k") * " dimensions"))+
  ylab("AUC") +
  guides(colour = guide_legend(title = "Algorithm"))+
  theme_bw(base_size = 16)
# facet_wrap(~Algorithm)

#finding which dimension the model that gave the highest auc score was trained on
posBestFit <- which.max(resDf$auc)
bestDimension <- resDf[which.max(resDf$auc),]$dimensions
for (d in 1:length(nmfDataMlResList)){
  #print(nmfRes[[d]])
  if(ncol(labelledGenesNMFRes[[d]])-1 == bestDimension){#allnmfmatrix is a list containing all the nmf reduced matrices +labels
    # print(allNmfMatrix[[d]])
    bestDf <- as.data.frame(labelledGenesNMFRes[[d]])
  }
}
#find the final fit for that model
bestFit <- nmfDataMlResList[[posBestFit]]$finalFit
#lets look at the variable importance
bestModVarImport <- nmfDataMlResList[[posBestFit]]$importanceDf
bestModVarImport <- bestModVarImport %>% 
  mutate(sign = case_when(Importance<0 ~"negative", TRUE~"positive"))

#subset for the ones that are positive sign
impFeats <- bestModVarImport %>% filter(sign=="positive")
#then extract the variables that contribute to model's predictions
modelFeats <- impFeats$Variable
#find those features 
modelFeatsDf <- bestDf %>% dplyr::select(all_of(modelFeats))# these are gonna be input for GSEA 
colnames(modelFeatsDf) <- sub('V', 'Feature', colnames(modelFeatsDf))


##
# modelFeatsDf <- modelFeatsDf %>% dplyr::mutate(sum = rowSums(modelFeatsDf))
library(gplots)
?heatmap.2()
# heatmap(mat)
mat <- as.matrix(modelFeatsDf)
library(future)
library(gplots)

plan(multisession, workers = availableCores())
pheatmap::pheatmap(mat)
hm <- heatmap.2(x = mat, 
          col = RColorBrewer::brewer.pal(9, c("RdBu")), 
          #col="bluered",
          dendrogram = "none", 
          Rowv = F, 
          Colv = F,
          tracecol = NA,
          rowsep = 1:nrow(mat),
          #colsep = 1:ncol(mat)-1,
          trace = 'none')
# ggplot(modelFeatsDf, aes(x = , y= ))+geom_
colnames <- paste0("Feature", 1:ncol(bestDf))
unlabelledGenes$prediction <- predict(bestFit, unlabelledGenes)
ggplot(bestModVarImport, aes(x=Variable, y=Importance, col = sign))+
  geom_point(size=4)+
  theme_bw(base_size = 18)+
  theme(legend.position="none")+
  xlab("")

# vip(bestModVarImport, geom = 'point', mapping=aes(colour = sign))
##

# get importance

getImportanceFctn <- function(mlResList, algorithm){
  
  impFeatResult <- mlResList$importanceDf
  importanceDfs <- do.call("bind_rows",impFeatResult)
  return(impFeatResult)
}#)
doCall <- do.call(bind_cols, resultingDf)

resultingDf<- lapply(nmfDataMlResList, FUN = getImportanceFctn, algorithm = "NMF")

####clusterprofiler code####
# enrichKEGG(gene = row.names(bestDf))
geneList <- data.frame(geneID = row.names(modelFeatsDf), Weight = modelFeatsDf$V8)
geneList_2<- as.data.frame(t(geneList))
genes <- c(row.names(modelFeatsDf))
weights <- as.vector(modelFeatsDf$V8)
weights <- sort(weights,decreasing = T)
weights <- as.vector(weights)
names(weights) <- genes
class(weights)
#genelis <- as.vector(geneList)
# class(genelis)
# enr
# enrichKEGG(gene = weights, 
#            organism = "hsa", keyType = )

gseGO <- gseGO(geneList = weights,
               ont = "BP",
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               maxGSSize=2000,
               eps=0,
               #pAdjustMethod="BH"
)
gseGO@result%>%
  ggplot(aes())
p<- dotplot(gseGO)
p$data %>% #filter(p.adjust<0.3)%>%
  ggplot(aes(x=GeneRatio,y=forcats::fct_reorder(Description, GeneRatio)))+#, colour = p.adjust, size =Count))+
  geom_segment(aes(xend=0, yend = Description))+
  geom_point(aes(colour=p.adjust, size = Count))+
  
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(1, 10)) +
  #scale_color_gradientn(colours = c(pal))+
  # scale_colour_gradientn(colours=c("coral", "tomato", "firebrick2", "firebrick3","slateblue4","plum4", "rosybrown", "plum3" ))+
  theme_bw(base_size = 14)+
  ylab("")

gseKEGG <- gseKEGG(geneList = weights,
                   organism = 'hsa',
)
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
