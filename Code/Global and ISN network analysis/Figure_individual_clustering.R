## Run after "Individual_ConsensusClustering.R"

source("DMM_analysis.R")
col.pal = RColorBrewer::brewer.pal(9, "Reds")

phenotype = data.frame(delivery_type = Children$delivery_type.x,
                       Diet6 = as.character(Children$Diet_TP.x),
                       Diet9 = as.character(Children$Diet_TP.y),
                       DMM = as.character(Children$DMM_clust.x),
                       DMM6 = as.character(Children$DMM_6_9M_Clusters.x),
                       DMM9 = as.character(Children$DMM_6_9M_Clusters.y))

rownames(phenotype) = Children$Child

phenotype$delivery_type = ifelse(phenotype$delivery_type == "Vaginal", 0,1)
phenotype$Diet = as.numeric(phenotype$Diet6 == phenotype$Diet9)
phenotype$DMM = as.numeric(phenotype$DMM)
phenotype$DMM6 = as.numeric(phenotype$DMM6)
phenotype$DMM9 = as.numeric(phenotype$DMM9)
phenotype$Clusters = clusters

plot(cluster::silhouette(clusters,dists))
plot(Tree, labels = Children$Childists)

rownames(cm) = Children$Child
colnames(cm) = Children$Child

ann_colors = list(
  DMM6 = c("white", "firebrick"),
  DMM9 = c("white", "firebrick"),
  DMM = RColorBrewer::brewer.pal(9, "Set3"),
  Diet = RColorBrewer::brewer.pal(5, "Set2"),
  delivery_type = c("white", "firebrick")
)

pheatmap::pheatmap(cm[Tree$order,Tree$order],
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = phenotype[Tree$order,],
                   annotation_row = phenotype[Tree$order,],
                   annotation_colors = ann_colors,
                   color = col.pal,
                   show_rownames=F, show_colnames=F,
                   fontsize = 6.5, fontsize_row=6,fontsize_col = 6)

my_order = order(phenotype$DMM9)
pheatmap::pheatmap(cm[my_order,my_order],
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = phenotype[my_order,],
                   annotation_row = phenotype[my_order,],
                   annotation_colors = ann_colors,
                   color = col.pal,
                   show_rownames=F, show_colnames=F,
                   fontsize = 6.5, fontsize_row=6,fontsize_col = 6)

## Find associations
phenotype$Diet[is.na(phenotype$Diet)] = 0
pred = prediction(phenotype$Diet, phenotype$Clusters)
performance(pred, measure = "auc")@y.values
T = table(phenotype[,c("Diet", "Clusters")])
chisq.test(T, simulate.p.value = TRUE)

## Compare with DMM
table(clusters)
table(phenotype$Clusters, phenotype$DMM6)
table(phenotype$Clusters, phenotype$DMM9)
table(paste0(phenotype$DMM6, phenotype$DMM9), phenotype$Clusters)

## imp nodes

# Dist = readRDS("Data/dist_MAGMA_Rep50_D10_L2r0001.rds")

Rep_dist = length(Dist)
Qval = c()
y = factor(clusters)
for (iii in 1:Rep_dist){
  dist = Dist[[iii]]
  X_main = dist
  X_main[is.na(X_main)] = 0
  
  FeaturePval = apply(X_main, 2, function(x){
  a = x[y==levels(y)[1]]
  b = x[y==levels(y)[2]]
  return(t.test(a,b)$p.value)
  })

  # hist(FeaturePval, 50)
  
  # Qval = rbind(Qval, order(p.adjust(FeaturePval, method = "BH"), decreasing = TRUE))
  Qval = rbind(Qval, order(FeaturePval))
}

rank_sum_mic = rep(0,95)
for (i in 1:95)
  for (j in 1:Rep_dist)
    rank_sum_mic[i] = rank_sum_mic[i] + which(Qval[j,]==i)

plot(sort(rank_sum_mic), pch = 20, cex = .5)
impp = order(rank_sum_mic)[1:7]
impp

MicrobNames = readRDS("Data/MicrobeList.rds")
ClusterImpMic = MicrobNames[impp]


