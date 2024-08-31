library(tidymodels)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(cowplot)
library(ggsci)
library(correctR)
#and execute the scripts needed to carry out validation of model on external dataset
# source("code/processHPOgenes.R") # only problem is biomart package is required for gene map
# source("code/processMGIgenes.R")

# read in dataused for training the models
dataForTraining <- readRDS("processed/SVwithNetworkLabelled.rds")

# read in the results of network alone!
networkOnlyXGBres <- readRDS("processed/xgbResNetworkOnly_500.rds")
networkOnlyXGBres$AUC
#0.515
networkOnlyXGBres$roc_curve
# and the 2nd result (to predict mortality)
networkOnlyXGBresMortality <- readRDS("processed/xgbResMortalityNetworkOnly_500.rds")
networkOnlyXGBresMortality$AUC
#0.775
#gives confidence that network is right

#let's plot the roc for results on held-out test set for predicting both phenotypes
#skeletal and mortality/aging
networkOnlyXGBres$aug <- networkOnlyXGBres$aug %>% mutate(Phenotype = "Skeletal")
networkOnlyXGBresMortality$aug <- networkOnlyXGBresMortality$aug %>%
  mutate(Phenotype = "Mortality/aging")

networkOnlyPreds <- rbind(networkOnlyXGBres$aug, networkOnlyXGBresMortality$aug)
networkROCdf <- networkOnlyPreds %>% group_by(Phenotype) %>%
  roc_curve(truth = significant, .pred_FALSE)
netAlonecurve <- ggplot(networkROCdf, aes(x=1-specificity, y = sensitivity, 
                                  colour = as.factor(Phenotype)))+
  geom_line(linewidth = 1.5)+
  geom_abline(slope = 1, intercept  = 0, linewidth = 0.5, lty="dashed", alpha = 0.7)+
  theme(panel.border = element_rect(colour = "black", linewidth = 1, fill = "white"),
        aspect.ratio = 1)+
  theme_bw(base_size = 28)+

  scale_colour_manual("Phenotype", values = c("#00A087B2", "#7D5CC6B2"))+
  
  facet_wrap(~"XGB on 500 dimension network")+
  theme(strip.text = element_text(size = 32), legend.position = "bottom")
  
netAlonecurve
ggsave("output/networkAloneROCauc.png", netAlonecurve, height = 12, width = 11)


###read in the results of tissue network
rfNetworkPlusSVnmfFiles <- list.files("processed", full.names = T, pattern = "ResNetworkWithSV_")
rfNetworkPlusSVnmfResList <- lapply(rfNetworkPlusSVnmfFiles, readRDS)
# put them in order
rfNetworkPlusSVnmfResList <- rfNetworkPlusSVnmfResList[order(sapply(rfNetworkPlusSVnmfResList, function(x) x$dim))]
# then get all the Cv metrics
cvRFNetSV <- lapply(rfNetworkPlusSVnmfResList, function(x){
  x$resDf
})
#READ in the nonnetwork for comparison (created in 04_ script )
cvNonetwork <- readRDS("processed/cvAllRf_noNetworkRes.rds")
cvNonetwork <- cvNonetwork %>% mutate(Feature = "gene expression only")
# then plot them together for comparison
cvNetSv <- do.call("rbind", cvRFNetSV) %>% mutate(model = "RF", Feature = "gene expression with network")
networkandNoNetwork <- bind_rows(cvNonetwork, cvNetSv)
networkSVcvBP <- ggplot(networkandNoNetwork, aes(x= as.factor(dim), y = mean, fill=Feature)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#4DBBD5B2", "#7E6148B2"))+
  labs(x="Number of NMF dimensions", y="Cross-validation AUC", fill= "") +
  theme_cowplot(font_size = 26)+
  facet_wrap(~"RF")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4), 
        strip.background = element_rect(color = "black",  linewidth = 0.4),
        legend.position = "bottom")+
  ggtitle("Tuning NMF dimension size with RF")+
  theme(plot.title = element_text(size = 16))
  # ggtitle(expression(atop(italic("Tuning NMF dimension size with RF"))))
#networkSVcvBP

# read in the other plots (created in 04_ script p1 is rf vs xgb and p2 is pca vs nmf)
# so can plot them on one page
p1 <- readRDS("output/p1_RFvsXGB.rds")
p1 <- p1 +
  #ggtitle(expression(atop(italic(""))))
  ggtitle("Comparing RF and XGB")+
  theme(plot.title = element_text(size = 16))
p1
p2 <- readRDS("output/p2_pcaVsNmf.rds")
p2 <- p2 +
  ggtitle("Comparing NMF and PCA")+
  #ggtitle(expression(atop(italic("Comparing NMF and PCA"))))+
  theme(plot.title = element_text(size = 16))
p2
bot <- plot_grid(p1,p2, labels = c("b","c"), label_size = 20)
bot
# allp <- plot_grid(networkSVcvBP / (p1+p2), labels = "auto", label_size = 20)
allp <- plot_grid(networkSVcvBP, bot, ncol = 1, labels = "auto", label_size = 20)
allp
#ggsave("output/allp.png", allp, height = 10, width = 11)
#
#ggsave("output/allp2.png", allp, height = 10, width = 11)
ggsave("output/allp3.png", allp, height = 10, width = 11)

#
cvRfMetrNET_SV <- lapply(rfNetworkPlusSVnmfResList, function(x){
  r <- x$resDf
  maxAuc <- as.data.frame(r[which.max(r$mean),])
})
maxAucNetDims <- do.call("rbind", cvRfMetrNET_SV) %>%
  select(mean, dim, lowerbound, upperbound, Algorithm)
## although the cv metrics are quite better on the training data, they do not reflect as 
#well on the test data:...


# now get all the results for on the held-out test set so we can plot them too as facetted
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
  geom_abline(slope = 1, intercept = 0, linewidth=0.7, lty="dashed", alpha = 0.8)+
  theme(panel.border = element_rect(colour = "black", linewidth = 1.0, fill="white"),
        aspect.ratio = 1, legend.position="none")+
  theme_minimal_grid(font_size = 30)+# theme(legend.position = "none")+
  scale_colour_npg()+
  # facet_wrap(~ paste0(dim, " NMF dimensions"))+
  facet_wrap(~ factor(paste0(dim, " NMF dimensions"), c("50 NMF dimensions", "100 NMF dimensions","150 NMF dimensions",
                                                        "200 NMF dimensions", "500 NMF dimensions")))+
  theme(strip.text = element_text(size=32), legend.position = "none", 
        #this allowed me to increase text size of x and y axes labels!
        text = element_text(size = 48))+
  geom_text(data = dfAUCwNet_text, mapping = aes(x=0.3, y=0.85, label = paste0("AUC = ", round(auc, 3))), 
            size=12)

svWnetROCrf
ggsave("output/SVnetworkROCcurveRF.png", svWnetROCrf, width = 20, height = 12)

# find which dataset the model with the best roc on the held-out test set was 
# trained on

posBestFitNet <- which.max(dfAUCwNet_text$auc)
posBestFitNet
bestDimensionNet <- dfAUCwNet_text[which.max(dfAUCwNet_text$auc),]$dim
bestDimensionNet
# dataForTraining <- sVwithNetworkLabelled
for (n in 1:length(rfNetworkPlusSVnmfResList)){
  if(ncol(dataForTraining[[n]])-1 == bestDimensionNet){#allnmfmatrix is a list containing all the nmf reduced matrices +labels
    bestDfwNet <- as.data.frame(dataForTraining[[n]])
  }
}
dim(bestDfwNet)
bestFitNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$finalFit # what are its hyperparams
bestFitNet

#get the df ready for the corrected repeated kfold cv test
correctRdfBestDimNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$dfForCorrectR
#get CV metrics
bestModMetricsNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$resDf
finalModMetricsNet <- bestModMetricsNet[which.max(bestModMetricsNet$mean),]


# now for actually interpreting the model's most important features

#lets look at the variable importance
bestModVarImportNet <- rfNetworkPlusSVnmfResList[[posBestFitNet]]$importanceDf
rfNetworkPlusSVnmfResList[[posBestFitNet]]$importancePlot
head(bestModVarImportNet)
# now get the first 4 features that are important to the model
bestModVarImportNet <- bestModVarImportNet %>% 
  mutate(sign = case_when(Importance<0 ~"negative", TRUE~"positive"))
head(bestModVarImportNet)
dim(bestModVarImportNet)
top4FeatsNet <- bestModVarImportNet[1:4,]
fourVarsNet <- top4FeatsNet$Variable

#find these vairables and their weightings of genes in the best dataframe (the one we found before
#that gave the best AUC on the held-out test set)
top4FeatsModeldfNet <- bestDfwNet %>% dplyr::select(all_of(fourVarsNet))
head(top4FeatsModeldfNet)
# change the column names to something more informative
colnames(top4FeatsModeldfNet)<- sub("X", "Feature", colnames(top4FeatsModeldfNet))
head(top4FeatsModeldfNet)


dim(top4FeatsModeldfNet)
# labelsFullDf <- read.table("processed/labelsFullDf.txt", header = T, sep = "\t")
# str(labelsFullDf$significant)
labelsAlone <- dataForTraining[[1]] %>% select(significant)
head(labelsAlone)
labelsAlone$significant<- as.factor(labelsAlone$significant)
all(rownames(top4FeatsModeldfNet)%in% rownames(labelsAlone))

all(rownames(top4FeatsModeldfNet)== rownames(labelsAlone))
#they're in order so can just do cbind
top4FeatsModeldfLabelledNet <- cbind(top4FeatsModeldfNet, labelsAlone)

# df <- top5FeatsModeldfLabelledNet[order(top5FeatsModeldfLabelledNet$Feature76, decreasing = T),]
# ggplot(df[1:15,], aes(x=significant, y = Feature76))+
#   geom_boxplot()
# plot a boxplot of the 4 features looking at if there is much of difference of 
# weightings between genes that are associated and not associated with a skeletal disease
for (m in 1:(ncol(top4FeatsModeldfLabelledNet)-1)) {
  #print(m)
  feature <- colnames(top4FeatsModeldfLabelledNet)[m]
  n <- sprintf("%s_nmfSvWithNet_%s",m,colnames(top4FeatsModeldfLabelledNet)[m])
  p <- ggplot(top4FeatsModeldfLabelledNet)+ 
    aes(y = top4FeatsModeldfLabelledNet[,m], x= significant, fill=significant)+
    geom_boxplot()+
    xlab("") + ylab(paste0("NMF " ,colnames(top4FeatsModeldfLabelledNet)[m]))+
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
# don't have to include every single plot though



# now get the dataframe ready to plot as a heatmap, rank it by the 
# model's most important feature
top4FeatsModeldfNet <- top4FeatsModeldfNet[order(top4FeatsModeldfNet[,1], decreasing = T),]
hGenesSymbs <- read.table("processed/human_coding_genes.txt", sep = "\t", header = T)
hGenesSymbs <- hGenesSymbs %>% 
  filter(hgnc_symbol!=""); rownames(hGenesSymbs) <- hGenesSymbs$ensembl_gene_id
hGenesSymbs$ensembl_gene_id <-NULL
head(hGenesSymbs) 
# so now gonna merge with top feats DF
head(top4FeatsModeldfLabelledNet)
topFeatsMapped <- merge(top4FeatsModeldfNet, hGenesSymbs, by=0); rownames(topFeatsMapped) <- topFeatsMapped$hgnc_symbol; topFeatsMapped$Row.names<-NULL; topFeatsMapped$hgnc_symbol <-NULL
head(topFeatsMapped)
topFeatsMapped <- topFeatsMapped[order(topFeatsMapped[,1], decreasing = T),]
head(topFeatsMapped)
# ?rank
matN <- as.matrix(topFeatsMapped)
head(matN)
nrow <- 15
matN <- apply(-matN, 2, rank)
# italicisee gene names!!!
italicNames <- lapply(
  rownames(matN[1:nrow,]), function(x) bquote(italic(.(x)))
)

head(matN)
# brewer.pal.info
heatmapN <- pheatmap::pheatmap(matN[1:nrow,], border_color = "white",
                               cluster_rows = F, 
                               cluster_cols = F, 
                               show_rownames = T, 
                               fontsize = 20, 
                               color = rev(brewer.pal(8, "Reds")),
                               display_numbers = T,
                               number_format = "%.0f", 
                               number_color = "black",
                               fontsize_number = 22,
                               labels_row = as.expression(italicNames)
)

heatmapN
ggsave("output/heatmap.png", heatmapN, width = 11, height = 6)

####extract each feature as a dataframe and put in a list to get ready for clusterprofiler####

getEachFeature_fctn <- function(dataframeOfFeats){
  single <- list()
  for (j in 1:ncol(dataframeOfFeats)){ 
    print(j)
    single[[j]] <- data.frame(dataframeOfFeats[,j],
                              row.names = rownames(dataframeOfFeats))
    colnames(single[[j]]) <- colnames(dataframeOfFeats)[j]
    #single[[j]] <- single[[j]][order(single[[j]][,1], decreasing = T),] #%>% arrange(desc(colnames(single[[j]])))
    single[[j]] <- single[[j]] %>% arrange(desc(colnames(single[[j]])))
  }
  return(single)
}

eachNetFeat <- getEachFeature_fctn(top4FeatsModeldfNet)
saveRDS(eachNetFeat, "processed/top4NetFeats.rds")
# saveRDS(eachNetFeat, "processed/top5NetFeats.rds") # do the gsea on my laptop cause clusterprofiler 
#not installing




############
####let's use the final model to rank the genes that are not yet labelled (by the IMPC####)

colnamesNet <- paste0("Feature", 1:ncol(bestDfwNet))
unlabelledGenes <- readRDS("processed/unStudiedGenes.rds")
head(unlabelledGenes[1:4,1:5])
colnames(unlabelledGenes) <- sub('V', 'X', colnames(unlabelledGenes))
head(unlabelledGenes[1:4,1:3])
unlabelledGenesPredNet <- augment(bestFitNet, unlabelledGenes)


unlabelledGenesPreddfNet <- as.data.frame(unlabelledGenesPredNet); rownames(unlabelledGenesPreddfNet) <- rownames(unlabelledGenesPredNet)
head(unlabelledGenesPreddfNet[1:3,1:4])
length(which(unlabelledGenesPreddfNet$.pred_class=="TRUE")) # 1059
length(which(unlabelledGenesPreddfNet$.pred_class!= "TRUE")) # 7942


########external validation using MGI and HPO

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
mgiCurve <- mgiCurve %>% mutate(database = "MGI")
head(mgiCurve)
rocaucMGI <- roc_auc(unlabelledGenesWithMGIannot, truth = significant, .pred_FALSE)
rocaucMGI <- rocaucMGI %>% mutate(database = "MGI")
autoplot(mgiCurve)
#
head(mgiCurve)

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
head(hpoCurve)
rocaucHPO <- roc_auc(unlabelledGenesWithHPOannot, significant, .pred_FALSE)
rocaucHPO <- rocaucHPO %>% mutate(database = "HPO")
rocaucHPO
#
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
write.table(as.data.frame(top20genesWsymbsNet[,1:3]), "processed/top20rankedgenes", col.names = T, sep = "\t", quote = F)


mgiHPhen <- read.table("processed/mgiHumanPhenotype.txt", header = T, sep = "\t")
hpoHPhen <- read.table("processed/hpoHumanPhenotype.txt", header = T, sep = "\t")
# iwant to find all the rownames of the top 20 that are also in column 1 of mgi
# s <- which(mgiHPhen[,1] %in% rownames(top20genesWsymbsNet))
# s
#so looking in the top 50 genes (predicted by the model) will look for overlaps in mgi/hpo db
top50genes <- orderedGenesNet[1:50,]
top50genesWsymbsNet <- merge(top50genes, hGenesSymbs, by=0)
rownames(top50genesWsymbsNet) <- top50genesWsymbsNet$hgnc_symbol; top50genesWsymbsNet$Row.names <-NULL; top50genesWsymbsNet$hgnc_symbol<-NULL
top50genesWsymbsNet <- top50genesWsymbsNet[order(top50genesWsymbsNet$.pred_TRUE, decreasing = T),]

#first looking in mgi
overlapMGI <- mgiHPhen[mgiHPhen[,1]%in% rownames(top50genesWsymbsNet),]#there are 3
overlapMGI$humanMarkerSymbol
overlapHPO <- hpoHPhen[hpoHPhen[,2] %in% rownames(top50genesWsymbsNet),]#also 3 in hpo
overlapHPO$geneSymbol

# what is the number of genes overlapping in both databases???????????
commonHPOMGI <- intersect(rownames(hpOnt), rownames(mgiGenes))
length(commonHPOMGI)


## I want to compare the best model that was trained on the gene expression data
# with the best from the network data:
NMFRFResSV_50 <- readRDS("processed/NMFRFResSV_50.rds")
NMFRFResNetworkWithSV_200 <- readRDS("processed/NMFRFResNetworkWithSV_200.rds")

nmfRFsv50cr <- NMFRFResSV_50$dfForCorrectR %>%
  mutate(model = "GE")
nmfRFnet200cr <- NMFRFResNetworkWithSV_200$dfForCorrectR %>%
  mutate(model = "Net")

compareNMFsvNetDF <- rbind(nmfRFsv50cr, nmfRFnet200cr)
compareNetWithSV_test<- repkfold_ttest(compareNMFsvNetDF, n1=80, n2=20, k=5, r=3)
compareNetWithSV_test
