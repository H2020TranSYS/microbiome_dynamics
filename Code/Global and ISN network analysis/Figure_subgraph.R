## Check different populations for the delivery type and diete mode
## Run after "Table_imp_Classification_microbes.R"

library(igraph)
library(ggraph)

Tab = readRDS(file = "Data/MERGED_Taxonomy_ASVs.rds")
rownames(Tab) = Tab$Node_list

ImpMic = DeliveryImpMic$Node_list
ImpMic = ClusterImpMic
DeliveryImpMic = Tab[ImpMic,3:9]

LuckiMap1 = read.delim("Data/LuckiMap_6M.txt")
LuckiMap2 = read.delim("Data/LuckiMap_9M.txt")
Children =  merge(LuckiMap1, LuckiMap2, by = "Child")

phenotype = data.frame(delivery_type = Children$delivery_type.x, 
                       DMM = as.character(Children$DMM_clust.x),
                       Diet6 = as.character(Children$Diet_TP.x),
                       Diet9 = as.character(Children$Diet_TP.y))

## y: delivery_type
y = phenotype$delivery_type
y[is.na(y)] = 0
y[y==0] = 1
y = as.factor(y)

## y: persistant diet
y = as.numeric(phenotype$Diet6 == phenotype$Diet9)
y[is.na(y)] = 0
y = as.factor(y)

## Read Microbiome Data
data = data.frame(read.delim(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
MianGraph1 = graph(Nodelist, directed = FALSE)
IndvGraphWeights1 = data[,-1]
IndvGraphWeights1 = abs(IndvGraphWeights1)
IndvGraphWeights1 = (IndvGraphWeights1 - min(IndvGraphWeights1)) / (max(IndvGraphWeights1) - min(IndvGraphWeights1))

data = data.frame(read.delim(file="Data/Resulting_net_notNULL_MAGMACONF9M.txt",sep = " " ))
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

# Graph weights are calculated as the average ISN weights
w1_vgn = apply(IndvGraphWeights1[, y==levels(y)[1]], 1, mean)
w1_csc = apply(IndvGraphWeights1[, y==levels(y)[2]], 1, mean) 

w2_vgn = apply(IndvGraphWeights2[, y==levels(y)[1]], 1, mean) 
w2_csc = apply(IndvGraphWeights2[, y==levels(y)[2]], 1, mean) 

G1_vgn = simplify(set.edge.attribute(MianGraph1, "weight", index=E(MianGraph1), w1_vgn))
G1_csc = simplify(set.edge.attribute(MianGraph1, "weight", index=E(MianGraph1), w1_csc))
G2_vgn = simplify(set.edge.attribute(MianGraph2, "weight", index=E(MianGraph2), w2_vgn))
G2_csc = simplify(set.edge.attribute(MianGraph2, "weight", index=E(MianGraph2), w2_csc))

A1_vgn = as.matrix(as_adj(G1_vgn,  attr = "weight"))
A1_csc = as.matrix(as_adj(G1_csc,  attr = "weight"))
A2_vgn = as.matrix(as_adj(G2_vgn,  attr = "weight"))
A2_csc = as.matrix(as_adj(G2_csc,  attr = "weight"))

N_imp_node = 7
ImpNodes = rownames(DeliveryImpMic)[1:N_imp_node]
# ImpNodes = ClusterImpMic
for (i in 1:N_imp_node){
  ImpNodes = c(ImpNodes, names(neighbors(G1_vgn, ImpNodes[i])))
  ImpNodes = c(ImpNodes, names(neighbors(G1_csc, ImpNodes[i])))
  ImpNodes = c(ImpNodes, names(neighbors(G2_vgn, ImpNodes[i])))
  ImpNodes = c(ImpNodes, names(neighbors(G2_csc, ImpNodes[i])))
}
Nodes = unique(ImpNodes)

A1_vgn_sub = A1_vgn[Nodes,Nodes]
A1_csc_sub = A1_csc[Nodes,Nodes]
A2_vgn_sub = A2_vgn[Nodes,Nodes]
A2_csc_sub = A2_csc[Nodes,Nodes]

# qgraph(A1_vgn)
# qgraph(A1_csc)
# qgraph(A2_vgn)
# qgraph(A2_csc)

dA_vgn = A2_vgn_sub - A1_vgn_sub
dA_csc = A2_csc_sub - A1_csc_sub
dd = abs(dA_vgn - dA_csc)
dd[dd<.029] = 0
# dd[dd<.025] = 0
nodes_final = names(which(apply(dd, 2, sum)>0))

# nodes_final = nodes_final[-c(10,13,16,19,20,21,22,23,24,25,26,27,28)]

dA_vgn = dA_vgn[nodes_final,nodes_final]
dA_csc = dA_csc[nodes_final,nodes_final]

# qgraph::qgraph(dA_vgn)
# qgraph::qgraph(dA_csc)

g1 = graph_from_adjacency_matrix(dA_vgn, weighted = TRUE, mode = "undirected")
g2 = graph_from_adjacency_matrix(dA_csc, weighted = TRUE, mode = "undirected")
dg = graph_from_adjacency_matrix(dA_csc - dA_vgn, weighted = TRUE, mode = "undirected")

NodeColor = rep("slategray2", length(nodes_final))
NodeColor[1:N_imp_node] = "slategrey"
TxtColor = rep("slategrey", length(nodes_final))
TxtColor[1:N_imp_node] = "black"
FontFace = rep("plain", length(nodes_final))
FontFace[1:N_imp_node] = "plain"
FontSize = rep(4, length(nodes_final))
FontSize[1:N_imp_node] = 4

pdf(file = "Figures/Subgraph_Delivery_diff_2.pdf")
g = dg
W = E(g)$weight
Names = Tab[names(V(g)),"Family"]
Names0= sapply(Names, function(x) {strsplit(x,"_")[[1]][3] })
Names = Tab[names(V(g)),"Genus"]
Names1 = sapply(Names, function(x) {strsplit(x,"_")[[1]][3] })
Names = Tab[names(V(g)),"Species"]
Names2 = sapply(Names, function(x) {strsplit(x,"_")[[1]][3] })

Names2[is.na(Names1)] = Names0[is.na(Names1)]
Names1[is.na(Names1)] = "unclass."
Names2[is.na(Names2)] = "sp."
Names = paste(Names1, Names2)


ggraph(g, layout = 'linear', circular = TRUE) +
  geom_edge_arc(aes(color = ifelse(W>0, "hotpink", "seagreen2"), edge_width = abs(W)^3), edge_alpha = .7) +
  geom_node_point(aes(x = x, y = y), size = 5, fill = NodeColor, shape = 21, colour = "black", stroke = .4) +
  geom_node_text(aes(x = x*1.12, y = y*1.12, angle = -((-node_angle(x, y)+90)%%180)+90,
                    hjust= c(rep(0,14), rep(1,14)),
                    # hjust= c(rep(0,11), rep(1,11)),
                    ), fontface = FontFace,
                 label = Names, size = FontSize, color = TxtColor) +
  scale_size_identity() +
  scale_edge_width_continuous(range = c(.0001, 3)) + 
  theme_void() +
  theme(legend.position="none") + 
  expand_limits(x = c(-2.5, 2.5), y = c(-2.5, 2.5))

dev.off()
