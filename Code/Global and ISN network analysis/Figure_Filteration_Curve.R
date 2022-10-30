## Filteration Curve
SEign = function(graph, N){
  Lsym = as.matrix(laplacian_matrix(graph, normalized = TRUE))
  Lsym[is.na(Lsym)] = 0
  PCA = prcomp(Lsym)
  Lambda = PCA$sdev[ncol(Lsym):1]
  return(sum(Lambda[1:N]))
}

Thr = seq(.01,.2,.01)

isn1_vgn = IndvGraphWeights1[, y==levels(y)[1]]
isn1_csc = IndvGraphWeights1[, y==levels(y)[2]] 

isn2_vgn = IndvGraphWeights2[, y==levels(y)[1]]
isn2_csc = IndvGraphWeights2[, y==levels(y)[2]]


ISN = isn1_vgn
MainGraph = MianGraph1

Nindv = ncol(ISN)
FV = matrix(0,Nindv,length(Thr))

for (indv in 1:Nindv){
  w_i = ISN[,indv]
  G_i = simplify(set.edge.attribute(MainGraph, "weight", index=E(MainGraph), w_i))
  A_i = as.matrix(as_adj(G_i,  attr = "weight"))
  # Adj = A_i[nodes_final,nodes_final]
  Adj = A_i
  
  for (i in 1:length(Thr)){
    Adj_bin = ifelse(Adj>Thr[i],1,0)
    g = graph_from_adjacency_matrix(Adj_bin, mode = "undirected")
    # FV[indv,i] = SEign(g,20)
    FV[indv,i] = (mean(degree(g)))
  }
}
FV1 = FV


ISN = isn2_vgn
MainGraph = MianGraph2

Nindv = ncol(ISN)
FV = matrix(0,Nindv,length(Thr))

for (indv in 1:Nindv){
  w_i = ISN[,indv]
  G_i = simplify(set.edge.attribute(MainGraph, "weight", index=E(MainGraph), w_i))
  A_i = as.matrix(as_adj(G_i,  attr = "weight"))
  # Adj = A_i[nodes_final,nodes_final]
  Adj = A_i
  
  for (i in 1:length(Thr)){
    Adj_bin = ifelse(Adj>Thr[i],1,0)
    g = graph_from_adjacency_matrix(Adj_bin, mode = "undirected")
    # FV[indv,i] = SEign(g,20)
    FV[indv,i] = (mean(degree(g)))
    
  }
}
FV2 = FV

df1 = data.frame(Thr = Thr, MeanFV = apply(FV1, 2, mean), sd = apply(FV1, 2, sd), cl = rep("1",length(Thr)))
df2 = data.frame(Thr = Thr, MeanFV = apply(FV2, 2, mean), sd = apply(FV2, 2, sd), cl = rep("2",length(Thr)))

df = rbind(df1,df2)

pdf(file = "Figures/FiltCurve_nonConsistantDiete.pdf", width = 3.5, height = 3.5)
ggplot(df, aes(x=Thr, y=MeanFV, group=cl, color=cl)) + 
  geom_line(size = .5) +
  geom_point() +
  geom_errorbar(aes(ymin=MeanFV-sd, ymax=MeanFV+sd), width=.005, size = .5,
                position=position_dodge(0.003))+
  labs(title="Persistence Diagrams", x="Threshols", y = "Node Degree")+
  theme(legend.title = element_blank()) + 
  theme_classic()+
  scale_color_manual(values=c("red", "green"))
dev.off()
