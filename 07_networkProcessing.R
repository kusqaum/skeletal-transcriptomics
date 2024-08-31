##
library(igraph)
library(stringr)
# first need biogrid data downloaded if i want to run this script!
# source("code/07a_getBiogridnetwork.R")

# biogridNetwork <- read.table("raw/bioGridNetworkData.txt", header = T)
biogridNetwork <- read.table("raw/allNetworkData.txt", header = T)
# make undirected
graph <- graph_from_data_frame(biogridNetwork, directed = F)


# simplify() to remove duplicate edges
graphS <- simplify(graph)

# get the largest connected component
components <- components(graphS)
groups(components)
# https://stackoverflow.com/questions/64344845/getting-the-biggest-connected-component-in-r-igraph
largestComp <- which.max(components$csize)

vertID <- V(graphS)[components$membership==largestComp]
head(vertID)
class(vertID)
subgraph <- igraph::induced_subgraph(graphS, vertID)
subgraph

# edgeList <- as_edgelist(subgraph, names = T)
edgeList <- as.data.frame(as_edgelist(subgraph, names = T))
dim(edgeList)
write.table(edgeList, "processed/networkEdgeList.edg", row.names = F, col.names = F, sep = "\t")
dim(edgeList)
head(edgeList)

# then use PecanPy in python to run node2vec on edge list



