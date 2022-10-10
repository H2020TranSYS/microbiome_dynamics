## Multiplex Graph Representation Learning and
## multiplex netwoek differential analysis
## All ISNs in one embedding (Preveius name: Main_GRL_ISN_v6.R)
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

for (indv in 1:M){
  t1 = Sys.time()
  print(indv/M*100)
  
  Threshold = 0
  a = IndvGraphWeights1[,indv]
  G1 = simplify(set.edge.attribute(MianGraph1, "weight", index=E(MianGraph1), a))
  G1 = delete_edges(G1, E(G1)[E(G1)$weight < Threshold])
  
  a = IndvGraphWeights2[,indv]
  G2 = simplify(set.edge.attribute(MianGraph2, "weight", index=E(MianGraph2), a))
  G2 = delete_edges(G2, E(G2)[E(G2)$weight < Threshold])
  
  
  ## Walker
  RW1 = RandomWalkRestart (G1, Nrep = 100, Nstep = 5, weighted_walk = TRUE)
  RW2 = RandomWalkRestart (G2, Nrep = 100, Nstep = 5, weighted_walk = TRUE)
  
  A1 = as.matrix(as_adj(G1,  attr = "weight"))
  A2 = as.matrix(as_adj(G2,  attr = "weight"))
  
  X = rbind(A1, A2)
  X = X / (apply(X, 1, sum) + .000000001)
  Y = rbind(RW1$P, RW2$P)
  
  individual_no = c(individual_no, rep(indv, nrow(X)))
  time_point = c(time_point, rep(6, nrow(X)/2), rep(9, nrow(X)/2))
  
  XX = rbind(XX, X)
  YY = rbind(YY, Y)
  t2 = Sys.time()
}

# RW_param = list(XX,YY,individual_no,time_point)
# saveRDS(RW_param, file = "Data/RW_1.rds")
List = readRDS("Data/RW.rds")
XX = List[[1]]
YY = List[[2]]
individual_no = List[[3]]
time_point = List[[4]]

N = ncol(XX)
M = max(unique(individual_no))

X = XX
Y = YY
Xtest = X

colnames(X) = paste0("V",as.character(1:ncol(X)))
colnames(Y) = paste0("V",as.character(1:ncol(Y)))
colnames(Xtest) = paste0("V",as.character(1:ncol(X)))

Rep = 30
Dist = list()
for (rep in 1:Rep){
  t1 = Sys.time()
  # latentSpace = EDNN(X, Y, Xtest = X, latentSize = 2, epochs = 10, batch_size = 100)
  latentSpace = EDNN(X ,Y, Xtest = X, latentSize = 10, epochs = 50, batch_size = 5, l2reg = .0001)
  
  plot(latentSpace[,1:2],pch = 20, cex = .1)
  
  dist = matrix(0,M,N)
  for (indv in 1:M){
    latentSpace_1 = latentSpace[individual_no==indv & time_point==6, ]
    latentSpace_2 = latentSpace[individual_no==indv & time_point==9, ]
    for (i in 1:N)
      dist[indv,i] = Distance (latentSpace_1[i,], latentSpace_2[i,], method = "cosine")
  }
  Dist[[rep]] = dist
  t2 = Sys.time()
}

saveRDS(Dist, file = "Data/dist_MAGMA_Rep30_D10_L2r0001.rds")

## A little test for classification
Dist_ranked = 0
for (rep in 1:Rep){
  dist = Dist[[rep]]
  dist_ranked = t(apply(dist, 1, function(x){order(x)/length(x)}))
  Dist_ranked = Dist_ranked + dist_ranked
}
dist = Dist_ranked / Rep
