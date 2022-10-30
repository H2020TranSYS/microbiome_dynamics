## Individual Clustering 
## Fist run "MNDA_ISNs_with_Repeat.R" to have a distances of different runs
## The results is a set of clusterings which then need to run
## "Individual_ConsensusClustering.R" for a consensus clustering.
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/Microbiome_final/")
source('~/Desktop/R_Root/Microbiome_final/MNDA_functions.R')
library(cluster)

Dist = readRDS("Data/dist_MAGMA_Rep50_D10_L2r0001.rds")
Rep = length(Dist)
maxNclust = 5

Clusters = c()
for (i in 1:Rep){
  dist_i = Dist[[i]]
  dists_cor = as.matrix(as.dist(1-abs(cor(t(dist_i)))))
  dists_cor[is.na(dists_cor)] = 0

  sil = rep(0,maxNclust)
  for (k in 2:maxNclust){
    clusters = kmeans(dists_cor,k)$cluster
    si = silhouette(clusters, dists_cor)
    sil[k] = mean(si[,"sil_width"])
  }
  
  k = which.max(sil)
  clusters = kmeans(dists_cor,k)$cluster
  Clusters = cbind(Clusters, clusters)
}

