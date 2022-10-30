## Consensus clustering for micribial clustering and find highly and lowly variable microbes
## Run after "Microbial_Clustering.R"
## To see the final results run 
## "Figure_microb_cluster.R" for micribial clustering, and
## "Figure_high_low_microbes.R" for highly and lowly variabel microbes
## By: Behnam Yousefi

## Function definition

HirClustFromAdjMat = function (Adj, k = 2, alpha = 1){
  Aff = exp(-(1-Adj)^2/(2*alpha^2))
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
MaxClust = 2
ResampleRatio = .7
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

# saveRDS(list(CM, coClusterMatrix),"Data/CM_111202.rds")

k = 2
clusters = HirClustFromAdjMat(CM[[k]], k)
table(clusters)

MicClust = cbind(clusters[1:95],clusters[96:(2*95)])

Tab = table(MicClust[,1],MicClust[,2])
print(Tab)

which((MicClust[,1]==1) & (MicClust[,2]==2))
which((MicClust[,1]==2) & (MicClust[,2]==1))
which((MicClust[,1]==2) & (MicClust[,2]==2))

# 
# ## CM analysis
List = readRDS("Data/CM_Global_MAGMA_D10_R50_Step5_maxClust6.rds")
CM = List[[1]]
coClusterMatrix = List[[2]]
 
Score = CC_cluster_count(CM)
 
# 
# MaxClust = length(CM)
# K = 2:MaxClust
# 
# par(new=FALSE)
# A = rep(0,length(K))
# S_diffA_Windowed = rep(0,length(K))
# Nbin = 100
# for (k in K){
#   M = CM[[k]]
#   ## Consensus distribution
#   M0 = M
#   diag(M0) = 100
# 
#   ConsDistr = rep(0,Nbin)
#   for (i in 1:Nbin)
#     ConsDistr[i] = sum(M0<(i/Nbin))/2
#   ConsDistr = ConsDistr/(Nsample*(Nsample-1)/2)
#   A[k] = sum(ConsDistr)/Nbin
# 
#   diffA = ConsDistr[2:Nbin]-ConsDistr[1:Nbin-1]
#   #diffA = 2*(diffA-min(diffA))/(max(diffA)-min(diffA)+.00000000001) - 1
# 
#   l = 3
#   diffA_smooth = rep(0,length(diffA))
#   for (n in (l+1):(length(diffA)-l))
#     diffA_smooth[n] =  mean(diffA[(n-l):(n+l)])
# 
#   diffA = diffA_smooth
# 
#   diffA_Windowed = (diffA*dnorm(0:(Nbin-2),0,10) +
#                       diffA*dnorm(0:(Nbin-2),Nbin,10)) - (diffA*dnorm(0:(Nbin-2),Nbin/2,40))
# 
# 
#   diffA_Windowed[is.infinite(diffA_Windowed)] = NA
#   diffA_Windowed[is.na(diffA_Windowed)] = max(diffA_Windowed, na.rm=TRUE)
# 
#   S_diffA_Windowed[k] = sum(abs(diffA_Windowed))/Nbin
# 
#   S_diffA_Windowed[k] = sum(diffA_smooth*dnorm(0:(Nbin-2),Nbin/2,20))
# 
#   plt = ConsDistr
#   #plt = diffA_smooth
#   #plt = diffA_Windowed
# 
#   plot(plt, type='l', col=k, ylim = c(0,1), lwd = 2,
#        xlab = "consensus matrix index", ylab = "CDF")
#   text(10*k, plt[10*k], labels = k, ylim = c(0,.04))
#   par(new=TRUE)
# }
# 
# par(new=FALSE)
# deltaA = A
# for (k in 3:(max(K)-1))
#   deltaA[k] = (A[k+1]-A[k])/A[k]
# plot(deltaA, type='b', pch = 20)
# 
# ## Final clusters
k = 2
clusters = HirClustFromAdjMat(CM[[k]], k)
table(clusters)

MicClust = cbind(clusters[1:95],clusters[96:(2*95)])

Tab = table(MicClust[,1],MicClust[,2])
print(Tab)
apply(Tab, 2, sum)
apply(Tab, 1, sum)

Mat = coClusterMatrix
diag(Mat) = 0
rownames(Mat) = 1:190
colnames(Mat) = 1:190
Variables = data.frame(Clusters = as.factor(clusters),
                       TimePoint = c(rep("6m",95),rep("9m",95)))
rownames(Mat) = 1:190
col.pal = RColorBrewer::brewer.pal(9, "Reds")

pheatmap::pheatmap(Mat[order(clusters), order(clusters)],
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = Variables[order(clusters),],
                   annotation_row = Variables[order(clusters),],
                   #color = col.pal,
                   show_rownames=F, show_colnames=F,
                   fontsize = 6.5, fontsize_row=6,fontsize_col = 6)

pheatmap::pheatmap(coClusterMatrix)

## Highly and lowly variable microbes
dist_cluster = rep(0,N)
for (i in 1:N)
  dist_cluster[i] = coClusterMatrix[i,N+i]

plot(sort(dist_cluster), cex = .5, pch = 20, xlab = "microbial index", ylab = "co-occurrence similarity")
abline(h=.99, col="red")
abline(h=.2, col="red")

LowVar_index = which(dist_cluster==1)
HighVar_index = order(dist_cluster)[1:10]
MicClust[LowVar,]
MicClust[HighVar,]
