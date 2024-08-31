# CODE

## Utility scripts

**01a_mapMouseToHumanEnsembl.R**   
convert mouse Ensembl identifiers to human Ensembl identifiers

**01b_extractHumanCodingGenes.R**   
extract all human coding genes from biomaRt

**01c_extractMouseSymbolsToEnsembl.R**   
read in and filter mouse to mouse homology from mouse genome informatics

**01c_extractMouseSymbolsToEnsembl.R**   
read in mouse to mouse homology from mouse genome informatics

**01d_mapIMPCtoHumanEnsembl.R**   
map IMPC mouse gene Ensembl identifiers to human Ensembl identifiers

**07a_getBiogridnetwork.R**   
read in and filter human protein-protein interactions from bioGRID to get ready for 07_networkProcessing.R script

Use PecanPy to run node2vec on edge list created in 07_networkProcessing.R script
install PecanPy with:
`pip3.9 install pecanpy`
run node2vec with:
`pecanpy --input networkEdgeList.edg --output networkEdgeList.emb --dimension 500`

**processMGIgenes.R**   
process human phenotype ontology genes for validating final model

**processHPOgenes.R**   
process mouse genome informatics genes for validating final model
