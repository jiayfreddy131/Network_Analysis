###### Building Cytoscape network in R #######
library(vegan)
library(reshape)
library(igraph)
library(fdrtool)
library(ggplot2)
library(Hmisc)
library(data.table)
library(dplyr)
library(rgexf)
library(psych)

## First, import "test" data, each column is a serie of parameter 
## (nutrient, temp etc.) or data (relative abun, cell abun, etc.)

## Create correlation matrix
matrix_dist_test <- rcorr(as.matrix(test), type = "spearman")
matrix_corr_test <- matrix_dist_test$r
matrix_p_test <- matrix_dist_test$P
matrix_p_test <- p.adjust(matrix_p_test, method = "BH")

## Change all unnecessary data to 0
# select for strong correlation (|r| >= 0.9) and significant (p < 0.05)
matrix_corr_test[which(matrix_corr_test>=-0.9 & matrix_corr_test<=0.9)] = 0
matrix_corr_test[which(matrix_p_test>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_test[is.nan(matrix_corr_test)] = 0

## Create network attributes
g_network_test <- graph.adjacency(matrix_corr_test, weighted = T, mode = "upper", diag = F)
e_test <- get.edgelist(g_network_test)
edges_test <- as.data.frame(cbind(e_test,E(g_network_test)$weight)) 
colnames(edges_test) <- c("source", "target", "weight")

nodes1_test <-as.data.frame(unique(edges_test[,1]));colnames(nodes1_test)<-c("id")
nodes2_test <-as.data.frame(unique(edges_test[,2]));colnames(nodes2_test)<-c("id")
nodes3_test=rbind(nodes1_test,nodes2_test)
nodes_test <-as.data.frame(unique(nodes3_test[,1]));colnames(nodes_test)<-c("id")

## delete unwanted interactions (in my cese, MF/MF & OTU/OTU)
# using dply filter to remove MF/MF and OTU/OTU corr
edges_test_filter <- dplyr::filter(edges_test, grepl(';', source))
edges_test_filter_2 <- dplyr::filter(edges_test_filter, !grepl(';', target))

## rebuild igraph with filtered edges
g_test <- graph.data.frame(d=edges_test_filter_2, vertices=nodes_test, directed=T)

library(RCy3)
# launch Cytoscape
cytoscapePing()
createNetworkFromIgraph(g_test, "TEST")