# Supervised machine learning to prioritise skeletal disease-associated genes using transcriptomics data

![image](https://github.com/user-attachments/assets/21a216a8-9b57-4235-a416-ef5fb213609d)

## Overview

It is important to identify yet unstudied skeletal disease-associated genes to aid in the discovery of potential therapeutic targets.   
This study aimed to train random forest and gradient boosted trees to prioritise skeletal disease-associated genes using binary mouse gene-to-phenotype associations.

## Installation

install necessary R packages by running:   
 `Rscript install/install.R`
## Directions

**phenotype_labels.R**   
read in IMPC labels and obtain binary labels for machine learning models

**01_geneExpression.R**   
processing of skeletal tissue-specific gene expression data

**02_exploratoryDataVisualisation.R**  
principal componet analysis and non-negative matrix factorisation of the processed gene expression data

**03_modelTraining.R**  
train random forest and gradient boosted trees to predict skeletal disease-associated genes

**04_modelInterpretation.R**  
identify model's most important features and most informative genes

**05_positiveControlGTEX.R**  
train gradient boosted trees to predict mortality/aging-associated genes (positive control)

**06_gtexPredictPhenotype.R**  
use GTEX gene expression data to predict skeletal disease-associated genes

**07_networkProcessing.R**  
process human protein-protein interaction network

**08_integrateNetwork.R**  
create tissue-specific network from expression data

**09_trainTissueNetwork.R**  
train random forest and gradient boosted trees on network

**10_tissueNetworkModelInterpretation.R**  
identify final model's most important features from tissue-specific network 

