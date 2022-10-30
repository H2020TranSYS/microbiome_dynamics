## Consensus clustering for individual clustering and 
## Run after "Individual_Clustering.R"
## To see the final results run 
## "Figure_microb_cluster.R" for micribial clustering, and
## "Figure_high_low_microbes.R" for highly and lowly variabel microbes
## By: Behnam Yousefi

## Function definition

HirClustFromAdjMat = function (Adj, k = 2, alpha = 1){
  # Aff = exp(-(1-Adj)^2/(2*alpha^2))
  Aff = Adj
  dists = as.dist(1-Aff)
  Tree = hclust(dists,method="ward.D")
  clusters = cutree(Tree, k=k)
  return(clusters)
}

ConnMat = function (Clusters){
  Nsample = length(Clusters)
  M = matrix(0,Nsample,Nsample)
  for (i in 1:Nsample)
    for (j in 1:Nsample)
      if (Clusters[i]*Clusters[j]>0)
        M[i,j] = ifelse(Clusters[i]==Clusters[j],1,0)
  return(M)
}

IndcMat = function (Clusters){
  Nsample = length(Clusters)
  I = matrix(0,Nsample,Nsample)
  for (i in 1:Nsample)
    for (j in 1:Nsample)
      I[i,j] = ifelse(Clusters[i]*Clusters[j]>0,1,0)
  return(I)
}


## Data

# X = data.frame(c1 = sample(3,100, replace = TRUE), 
#                c2 = sample(4,100, replace = TRUE), 
#                c3 = sample(4,100, replace = TRUE),
#                c4 = sample(2,100, replace = TRUE),
#                c5 = sample(3,100, replace = TRUE))
X = Clusters

Nsample = nrow(X)
Nmethod = ncol(X)

## Co-cluster matrix
M = matrix(0,Nsample,Nsample)
for (cl in 1:Nmethod)
  M = M + ConnMat(X[,cl])
coClusterMatrix = M/Nmethod
pheatmap::pheatmap(coClusterMatrix)

## Parameters
MaxClust = 6
ResampleRatio = .8
MaxIt = 100
CM = list()

## Consensus Clustering
for (nClust in 2:MaxClust) {                          # Loop for K
  ## Conectivity matrix
  M = matrix(0,Nsample,Nsample)                      
  rownames(M) = rownames(X)
  colnames(M) = rownames(X)
  ## Indicator matrix
  I = M
  
  for (i in 1:MaxIt){
    RandInd = sample(Nsample, floor(ResampleRatio*Nsample), replace = FALSE)
    coClusterMatrix_i = coClusterMatrix[RandInd,RandInd]
    
    ## Do clustering
    clusters = HirClustFromAdjMat(coClusterMatrix_i, nClust)
    
    Clusters = rep(0,Nsample)
    names(Clusters) = rownames(X)
    Clusters[RandInd] = clusters
    
    ## Conectivity matrix
    Mi = ConnMat(Clusters)
    M = M + Mi
    
    Ii = IndcMat(Clusters)
    I = I + Ii
  }
  
  CM[[nClust]] = M/I
}

# saveRDS(list(CM, coClusterMatrix),"Data/CM_Indv_MAGMA_D10_R50_Step5_L2r0001.rds")



## CM analysis
# List = readRDS("Data/CM_Indv_MAGMA_D10_R50_Step5_L2r0001.rds")
# CM = List[[1]]
# coClusterMatrix = List[[2]]

MaxClust = length(CM)
Nsample = ncol(CM[[2]])
K = 2:MaxClust

par(new=FALSE)
A = rep(0,length(K))
S_diffA_Windowed = rep(0,length(K))
Nbin = 100
for (k in K){
  M = CM[[k]]
  M0 = M
  diag(M0) = 100
  
  ConsDistr = rep(0,Nbin)
  for (i in 1:Nbin)
    ConsDistr[i] = sum(M0<(i/Nbin))/2
  ConsDistr = ConsDistr/(Nsample*(Nsample-1)/2)
  A[k] = sum(ConsDistr)/Nbin

  plt = ConsDistr
  plot(plt, type='l', col=k, ylim = c(0,1))
  text(10*k, plt[10*k], labels = k, ylim = c(0,.04))
  par(new=TRUE)
}

# Monti
par(new=FALSE)
deltaA = A
for (k in 3:(max(K)-1))
  deltaA[k] = (A[k+1]-A[k])/A[k]
plot(deltaA, type='b', pch = 20)

# Maybe better
par(new=FALSE)
deltaA = A
for (k in 3:(max(K)-1))
  deltaA[k] = (A[k]-A[k-1])/A[k-1]
plot(deltaA, type='b', pch = 20)

## Final clusters
k = 2
cm = CM[[k]]

dists = as.dist(1-cm)
Tree = hclust(dists,method="ward.D")
clusters = cutree(Tree, k=k)
table(clusters)

pheatmap::pheatmap(cm)
