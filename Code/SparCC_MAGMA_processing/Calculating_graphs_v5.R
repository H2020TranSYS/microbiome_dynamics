
# Libraries ---------------------------------------------------------------


library(Rcpp)
library(rMAGMA)
library(dplyr)
library(igraph)
library(stringi)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# name <- read.table("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/name.txt", quote="\"", comment.char="")
# individuals <- read.table("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/individuals.txt", quote="\"", comment.char="")
# Ind_list = cbind("Individual_ID"= name,"NAME"= paste0(rep("Ind_", nrow(column_Individual)),individuals$V1))
# colnames(Ind_list) = c("Individual_ID", "NAME")
# write.table(Ind_list,file = "mapping_nameIndividual.txt",quote = FALSE, row.names = FALSE, col.names = TRUE)

eliminate = T
# Read Paddy data ---------------------------------------------------------
source("C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Code/Graphical/FUNCTION_EVERY_ORDER_GRAPH_v2_see.R")

## DATAS in various folder 
data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
month6_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/6M/"
month9_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/9M/"
sparcc_data_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/SparCC_data/"

result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/6M"
graph_path =  "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/SparCC_graphs/"
# For printing the same layout as in the MAGMA
LAYOUT_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/Graphical_Layout/" #LAYOUT_PHYLUM_6M_graph.txt"



# name <- read.table("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/name.txt", quote="\"", comment.char="")
# individuals <- read.table("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/individuals.txt", quote="\"", comment.char="")
# Ind_list = cbind("Individual_ID"= name,"NAME"= paste0(rep("Ind_", nrow(column_Individual)),individuals$V1))
# colnames(Ind_list) = c("Individual_ID", "NAME")
# write.table(Ind_list,file = "mapping_nameIndividual.txt",quote = FALSE, row.names = FALSE, col.names = TRUE)

setwd(data_path)
# setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/6M")

w0_n335_OTU_table6M <- read.delim("SparCC_ASV6M.tab")
w0_n335_OTU_table9M <- read.delim("SparCC_ASV9M.tab")
w0_n335_OTU_table = merge(w0_n335_OTU_table6M, w0_n335_OTU_table9M, by = "X.OTU.ID")

rownames(w0_n335_OTU_table) = w0_n335_OTU_table[,1]


## Covariates
w0_n335_OTU_table = w0_n335_OTU_table %>% select(-X.OTU.ID)
LuckiMap_6M <- read.delim("LuckiMap_anonym_long_format_6m.txt")
LuckiMap_9M <- read.delim("LuckiMap_anonym_long_format_9m.txt")

w0_n335_metadata = rbind(LuckiMap_6M, LuckiMap_9M)

Taxonomy_ASVs = read.csv("Taxonomy_ASVs.csv")






# OTU TABLE 6m ------------------------------------------------------------
setwd(month6_path)
otu_table_6m = as.matrix(read.table("OTU_TABLE_MAGMA6M.tsv"))


# OTU TABLE 9m ------------------------------------------------------------
setwd(month9_path)
otu_table_9m = as.matrix(read.table("OTU_TABLE_MAGMA9M.tsv"))

# SparCC 6m associations -----------------------------------------------------
setwd(sparcc_data_path)



median_correlation_modified <- read.delim("GlobNet6M/median_correlation_modified.txt")
rownames(median_correlation_modified) = median_correlation_modified$OTU.ID
continuousresult6m = as.matrix(median_correlation_modified[,-1])
write.table(continuousresult6m, file = paste0("SparCC_continuous_Data_GLOBAL6M", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")


median_correlation_modified <- read.delim("GlobNet9M/median_correlation_modified.txt")
rownames(median_correlation_modified) = median_correlation_modified$OTU.ID
continuousresult9m = as.matrix(median_correlation_modified[,-1])
write.table(continuousresult9m, file = paste0("SparCC_continuous_Data_GLOBAL9M", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")



# GRAPHS ------------------------------------------------------------------
plotting_globalnetwork_grouped = function(otu_table, continuousresult, name, month){

  
  
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Phylum"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = month,
    import_layout = paste0(LAYOUT_path, "LAYOUT_PHYLUM_6M_graph.txt")
  )
  
  # CLAss -------------------------------------------------------------------
  
  
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Class"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = month,
    import_layout = paste0(LAYOUT_path, "LAYOUT_Class_6M_graph.txt"), 
    level_taxa_order = "Class"
    
  )
  
  
  # ORRRDER -----------------------------------------------------------------
  
  
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Order"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = month,
    import_layout = paste0(LAYOUT_path, "LAYOUT_Order_6M_graph.txt"), 
    level_taxa_order = "Order"
    
  )
  
  
  
  
  # FAMILY ------------------------------------------------------------------
  
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Family"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = month,
    import_layout = paste0(LAYOUT_path, "LAYOUT_Family_6M_graph.txt"), 
    level_taxa_order = "Family"
    
  )
  
}
plot_circleplot = function(otu_table, name,month, Taxonomy_ASVs = Taxonomy_ASVs){
  magma_Stool_AllCov <- magma(data = otu_table)
  
  Taxonomy_ASVs[2,]
  ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(otu_table),]
  df = (t(ASV_final))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  
  
  temp_graph = igraph::graph_from_adjacency_matrix(   magma_Stool_AllCov$refit, mode="undirected" )
  
  coords <- layout_in_circle(temp_graph ) 
  
  
  ll2 <- igraph::layout_with_fr(igraph::graph_from_adjacency_matrix(
    magma_Stool_AllCov$refit, mode="undirected"))
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[5,]) ) ) )
  
  # png("6m_graphs_Paddy_Covariates_MAGMA_5thlevel_FAMILY.png", width = 465, height = 225, units='mm', res = 300)
  # # magma_Stool_AllCov$
  # plot(magma_Stool_AllCov,layout = ll2, V.color.factor = df[5,])
  # dev.off()
  # 
  # png(month,"_graphs_Paddy_Covariates_MAGMA_5thlevel_FAMILY_circle_plot.png", width = 465, height = 225, units='mm', res = 300)
  # # magma_Stool_AllCov$
  # plot(magma_Stool_AllCov,layout = coords, V.color.factor = df[5,])
  # dev.off()
  
  png(paste0(month,"_MAGMA_", name, "_5thlevel_FAMILY_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord, V.color.factor = ifelse(stri_replace_all_regex(df[5,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[5,], "[a-z]__", "")))
  dev.off()
  
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[3,]) ) ) )
  
  png(paste0(month,"_MAGMA_", name, "_3rdlevel_CLASS_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord, V.color.factor = ifelse(stri_replace_all_regex(df[3,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[3,], "[a-z]__", "")))
  dev.off()
  
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[2,]) ) ) )
  
  png(paste0(month,"_MAGMA_", name, "_2ndlevel_PHYLUM_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord,V.color.factor = ifelse(stri_replace_all_regex(df[2,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[2,], "[a-z]__", "")))
  dev.off()
  
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[4,]) ) ) )
  
  png(paste0(month,"_MAGMA_", name, "_4thlevel_ORDER_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord, V.color.factor = ifelse(stri_replace_all_regex(df[4,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[4,], "[a-z]__", "")))
  dev.off()
}

checkStrict(plotting_globalnetwork_grouped)
checkStrict(plot_circleplot)


setwd(graph_path)
plot_circleplot(otu_table = otu_table_6m, name ="", Taxonomy_ASVs = Taxonomy_ASVs,month = "6M_")
plotting_globalnetwork_grouped(otu_table = otu_table_6m, continuousresult = continuousresult6m, name ="", month = "6M_")


plot_circleplot(otu_table = otu_table_9m, name ="", Taxonomy_ASVs = Taxonomy_ASVs,month = "9M_")
plotting_globalnetwork_grouped(otu_table = otu_table_9m, continuousresult = continuousresult9m, name ="", month = "9M_")



