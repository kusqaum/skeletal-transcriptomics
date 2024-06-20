#function to install libraries from bioconductor
installPackage <- function(libName) {
  if(libName %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(libName,ask = FALSE)
  }}

#install BiocManager first
install.packages("BiocManager")

#read the libraries needed
packagesToInstall <- read.delim("install.txt",header=FALSE, 
                                stringsAsFactors = FALSE)

#install all the libraries
sapply(packagesToInstall[,1],installPackage)