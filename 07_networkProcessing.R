##
library(igraph)
library(stringr)
# biogridNetwork <- read.table("raw/bioGridNetworkData.txt", header = T)
biogridNetwork <- read.table("raw/allNetworkData.txt", header = T)
graph <- graph_from_data_frame(biogridNetwork, directed = F)
vcount(graph)


graphS <- simplify(graph)

components <- components(graphS)
groups(components)
# https://stackoverflow.com/questions/64344845/getting-the-biggest-connected-component-in-r-igraph
largestComp <- which.max(components$csize)

vertID <- V(graphS)[components$membership==largestComp]
head(vertID)
class(vertID)
subgraph <- igraph::induced_subgraph(graphS, vertID)
subgraph

edgeList <- as_edgelist(subgraph, names = T)
write.table(edgeList, "processed/networkEdgeList.edg", row.names = F, col.names = F, sep = "\t")


# temp <- read.table("processed/networkEdgeList.edg", col.names = F)

