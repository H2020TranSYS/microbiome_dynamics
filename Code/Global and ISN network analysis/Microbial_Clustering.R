## Microbial clustering
## Fist run "MNDA_Global_with_Repeat.R" to have the latent spaces
## The results is a set of clustering which then need to run
## "Micriobial_ConsensusClustering.R" for a consensus clustering.
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/Microbiome_final/")
source('~/Desktop/R_Root/Microbiome_final/MNDA_functions.R')
library(cluster)

LatentSpace_list = readRDS("Data/LatentSpace_Global_MAGMA_D10_Rep20_Step5.rds")
Rep = length(LatentSpace_list)
N = nrow(LatentSpace_list[[1]])/2

maxNclust = 6
method = "cosine"
#method = "euclidean"

Dists_pairs = matrix(0,Rep,N)
Clusters = c()
for (i in 1:Rep){
  z = LatentSpace_list[[i]]
  
  ## distance measure for all microbes
  if (method == "cosine"){
    z = z / apply(z, 1, function(x){sqrt(sum(x^2))})
    Cos = z %*% t(z)
    dists = as.dist(1-Cos)
  }else
    dists = dist(z, method = method)
  
  ## cosine distance of the corresponding microbes
  latentSpace_1 = z[1:N,]
  latentSpace_2 = z[(N+1):(2*N),]
  for (n in 1:N)
    Dists_pairs[i, n] = Distance(latentSpace_1[n,], latentSpace_2[n,], method = "cosine")

  ## clustering
  sil = rep(0,maxNclust)
  for (k in 2:maxNclust){
    clusters = kmeans(dists,k)$cluster
    
    #Tree = hclust(dists,method="ward.D")
    #clusters = cutree(Tree, k=k)
    
    si = silhouette(clusters, dists)
    sil[k] = mean(si[,"sil_width"])
  }
  
  k = which.max(sil)
  clusters = kmeans(dists,k)$cluster
  Clusters = cbind(Clusters, clusters)
}

aa = Clusters[,1] - Clusters[,2]
sum(aa==0)/length(aa)

## Distances analysis
Dists_pairs[is.na(Dists_pairs)] = 0
Dists_pairs = t(apply(Dists_pairs, 1, function(x){order(x)/length(x)}))
dists = apply(Dists_pairs, 2, mean)
hist(dists,30)
plot(sort(dists), pch = 20, cex = .5)
order(dists)[1:4]
