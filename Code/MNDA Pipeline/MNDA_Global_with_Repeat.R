## Multiplex Graph Representation Learning and
## multiplex netwoek differential analysis
## works on global netwoeks 
## Graph Representation Learning
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/Microbiome_final/")
source('~/Desktop/R_Root/Microbiome_final/MNDA_functions.R')

library(igraph)
library(MASS)

## Read Microbiome Data
data = data.frame(read.delim(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " " ))
# data = data.frame(read.delim(file="Data/Resulting_net_from_corr_eliminate_SparCC_6M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
MianGraph1 = graph(Nodelist, directed = FALSE)
IndvGraphWeights1 = data[,-1]
IndvGraphWeights1 = abs(IndvGraphWeights1)
IndvGraphWeights1 = (IndvGraphWeights1 - min(IndvGraphWeights1)) / (max(IndvGraphWeights1) - min(IndvGraphWeights1))

data = data.frame(read.delim(file="Data/Resulting_net_notNULL_MAGMACONF9M.txt",sep = " " ))
# data = data.frame(read.delim(file="Data/Resulting_net_from_corr_eliminate_SparCC_9M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
MianGraph2 = graph(Nodelist, directed = FALSE)
IndvGraphWeights2 = data[,-1]
IndvGraphWeights2 = abs(IndvGraphWeights2)
IndvGraphWeights2 = (IndvGraphWeights2 - min(IndvGraphWeights2)) / (max(IndvGraphWeights2) - min(IndvGraphWeights2))

rm(data, Nodelist)

LuckiMap1 = read.delim("Data/LuckiMap_6M.txt")
LuckiMap2 = read.delim("Data/LuckiMap_9M.txt")
Children =  merge(LuckiMap1, LuckiMap2, by = "Child")

IndvGraphWeights1 = IndvGraphWeights1[,Children$Sample.name.x]
IndvGraphWeights2 = IndvGraphWeights2[,Children$Sample.name.y]

N = length(V(MianGraph1))
M = ncol(IndvGraphWeights1)

individual_no = c()
time_point = c()
XX = c()
YY = c()
XXtest = c()

# Graph weights are calculated as the average ISN weights
w1 = apply(IndvGraphWeights1, 1, mean) 
w2 = apply(IndvGraphWeights2, 1, mean)

Threshold = 0
G1 = simplify(set.edge.attribute(MianGraph1, "weight", index=E(MianGraph1), w1))
G1 = delete_edges(G1, E(G1)[E(G1)$weight < Threshold])

G2 = simplify(set.edge.attribute(MianGraph2, "weight", index=E(MianGraph2), w2))
G2 = delete_edges(G2, E(G2)[E(G2)$weight < Threshold])

LatentSpace_list = list()
Rep = 20
for (rep in 1:Rep) {
  
  ## Random walk algorithm
  RW1 = RandomWalkRestart (G1, Nrep = 100, Nstep = 5, weighted_walk = TRUE)
  RW2 = RandomWalkRestart (G2, Nrep = 100, Nstep = 5, weighted_walk = TRUE)
  
  ## Weights adjacency matrices are considered as inputs, and
  # Random walk node visits are considered as the outputs of the decoder
  A1 = as.matrix(as_adj(G1,  attr = "weight"))
  A2 = as.matrix(as_adj(G2,  attr = "weight"))
  
  X = rbind(A1, A2)
  X = X / (apply(X, 1, sum) + .000000001)
  Y = rbind(RW1$P, RW2$P)

  colnames(X) = paste0("V",as.character(1:ncol(X)))
  colnames(Y) = paste0("V",as.character(1:ncol(Y)))

  ## Train encoder-decoder neurak netwoek (EDNN), and
  # calculate the embedding (latent space) for all the nodes (Xtest = X)
  latentSpace = EDNN(X, Y, Xtest = X, latentSize = 10, epochs = 20, batch_size = 10)
  LatentSpace_list[[rep]] = latentSpace
  
  ## Plot
  # 1) select the two high var dimensions of the embedding
  Dim = order(apply(latentSpace, 2, var), decreasing = TRUE)[1:2]
  latentSpace_1 = latentSpace[1:N,Dim]
  latentSpace_2 = latentSpace[(N+1):(2*N),Dim]
  Max = apply(latentSpace_1, 2, max)
  Min = apply(latentSpace_1, 2, min)
  
  # 2) plot them
  par(new = FALSE)
  plot(latentSpace_1, cex = .2, pch = 20, col = "red", xlim = c(Min[1],Max[1]), ylim = c(Min[2],Max[2]))
  text(latentSpace_1, labels = V(G1), cex = .6, col = "red", xlim = c(Min[1],Max[1]), ylim = c(Min[2],Max[2]))
  par(new=TRUE)
  plot(latentSpace_2, cex = .2, pch = 20, col= "blue", xlim = c(Min[1],Max[1]), ylim = c(Min[2],Max[2]))
  text(latentSpace_2, labels = V(G1), cex = .5, col= "blue", xlim = c(Min[1],Max[1]), ylim = c(Min[2],Max[2]))
  
}

saveRDS(LatentSpace_list, file = "Data/LatentSpace_Global_MAGMA_D10_Rep20_Step5.rds")
