## Figure 4b
## Run after "Micriobial_ConsensusClustering.R"

library(ggplot2)
library(ggalluvial)

Tab = readRDS(file = "Data/MERGED_Taxonomy_ASVs.rds")
MicrobNames = readRDS("Data/MicrobeList.rds")

Tab2 = Tab[Tab$Node_list %in% MicrobNames,]
rownames(Tab2) = Tab2$Node_list

DynamicList = c()
Cluster = c()
Kcluster = max(MicClust)
for (i in 1:Kcluster){
  for (j in 1:Kcluster){
    Names = MicrobNames[MicClust[,1] == i & MicClust[,2] == j]
    DynamicList = c(DynamicList, Names)
    Cluster = c(Cluster, rep(paste0(i,"_",j), length(Names)))
  }
}

Tab2 = Tab2[DynamicList, 3:9]
Tab2$Cluster = Cluster

is_alluvia_form(as.data.frame(Tab2), silent = TRUE)
colors = RColorBrewer::brewer.pal(4,"Set2")[c(4,1,3,2)]

# pdf(file = "Figures/Figure4b_noLable.pdf", width = 7, height = 3.5)
ggplot(as.data.frame(Tab2),
       aes(axis1 = Phylum, axis2 = Class, axis3 = Order ,axis4 = Cluster)) +
  geom_alluvium(aes(fill = Cluster), col = NA, width = 1/6, alpha = .5, size = 0) +
  geom_stratum(aes(fill = Cluster), width = 1/6, fill = "seashell", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phylum", "Class", "Order", "Cluster"), expand = c(.05, .05)) +
  scale_fill_manual(values = colors) +
  theme_classic()
# dev.off()
