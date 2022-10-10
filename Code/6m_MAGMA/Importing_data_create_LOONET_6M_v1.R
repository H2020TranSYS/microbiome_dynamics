
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


data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/6M"
graph_path =  "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/MAGMA_graphs/6M"
LAYOUT_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/Graphical_Layout/" #LAYOUT_PHYLUM_6M_graph.txt"

setwd(data_path)
# setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/6M")

w0_n335_OTU_table6M <- read.delim("SparCC_ASV6M.tab")
w0_n335_OTU_table9M <- read.delim("SparCC_ASV9M.tab")
w0_n335_OTU_table = merge(w0_n335_OTU_table6M, w0_n335_OTU_table9M, by = "X.OTU.ID")

rownames(w0_n335_OTU_table) = w0_n335_OTU_table[,1]
# rownames(w0_n335_OTU_table)


## Covariates
w0_n335_OTU_table = w0_n335_OTU_table %>% select(-X.OTU.ID)
LuckiMap_6M <- read.delim("LuckiMap_anonym_long_format_6m.txt")
LuckiMap_9M <- read.delim("LuckiMap_anonym_long_format_9m.txt")

w0_n335_metadata = rbind(LuckiMap_6M, LuckiMap_9M)

Taxonomy_ASVs = read.csv("Taxonomy_ASVs.csv")




# Pre-selection based on OTUs present in more than 15% and samples --------

setwd(result_path)

t_otu = t(w0_n335_OTU_table)
## to have rows as individual amnd columns as taxa
prevalence       <- colMeans(t_otu>0)
sequencing_depth <- rowSums(t_otu)

icol <- prevalence >0.15 ##the prevalence value decreased to assure more OTUS
irow <- sequencing_depth >500
sum(icol)



otu_table <- t_otu[irow,icol] ## w0_n335_OTU_table[irow,icol]

dim(otu_table) #[1] 155  95



## covariates cutting
# final_covariate = first_covariates_selected[first_covariates_selected$FC.nummer %in% rownames(otu_table),]
# 
# final_covariate$FC.nummer == rownames(otu_table)
# MAGMA -------------------------------------------------------------------
# X = as.matrix(final_covariate %>% select(-c(FC.nummer)))
# X
# class(X)
magma_Stool_NoCov <- magma(data = otu_table)





# OTU TABLE 6M ------------------------------------------------------------
colnames(otu_table)
Children =  merge(LuckiMap_6M, LuckiMap_9M, by = "Child")


# MATCHING ----------------------------------------------------------------


# select paired samples  ----------------------------------------------------------

## ALL THE SAMPLES
otu_table6M = otu_table[rownames(otu_table) %in% LuckiMap_6M$Sample.name,]
sum(rownames(otu_table6M) %in% LuckiMap_6M$Sample.name)
# 81 ; 69 shared



## SHARED SAMPLES 
otu_table6M_69 = otu_table[rownames(otu_table) %in% Children$Sample.name.x,]
sum(rownames(otu_table6M) %in% Children$Sample.name.x)
# 69


magma_Stool_NoCov <- magma(data = otu_table6M)
write.table(otu_table6M, file = "OTU_TABLE_MAGMA6M.tsv",quote = F,sep = "\t",row.names = T, col.names = T)
write.table(otu_table6M_69, file = "OTU_TABLE_MAGMA6M_69.tsv",quote = F,sep = "\t",row.names = T, col.names = T)




plot(magma_Stool_NoCov)
png("Graphs__MAGMA6M.png", width = 465, height = 225, units='mm', res = 300)
plot(magma_Stool_NoCov)
dev.off()



# Example with the data of the library ------------------------------------


## CREATE LOO NETWORKs

setwd(result_path)
# mapping_nameIndividual_SPARCC <- read.csv("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/mapping_nameIndividual_SPARCC.txt", sep="")
calculating_LOO = function(table_otu, dirs, name){
  for (i in 1: nrow(table_otu)){
    # if (i > 5){break}
    ## does not work, interrupt on the 7th
    Individual = rownames(table_otu)[i]
    magma_Stool_AllCov <- magma(data = table_otu[-i,])
    index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda )
    binaryresult = magma_Stool_AllCov$path[[index]]
    colnames(binaryresult) = rownames(binaryresult) = colnames(table_otu)
    continuousresult = magma_Stool_AllCov$opt.icov
    colnames(continuousresult) = rownames(continuousresult) = colnames(table_otu)
    print(rownames(table_otu)[i])
    write.table(continuousresult, file = paste0(dirs, "/MAGMA_continuous_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,
                sep = "\t")
    write.table(binaryresult, file = paste0(dirs, "/MAGMA_binary_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,
                sep = "\t")
    
  }
  
  
  magma_Stool_AllCov <- magma(data = table_otu)
  index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda ) 
  binaryresult = magma_Stool_AllCov$path[[index]]
  colnames(binaryresult) = rownames(binaryresult) = colnames(table_otu)
  continuousresult = magma_Stool_AllCov$opt.icov
  colnames(continuousresult) = rownames(continuousresult) = colnames(table_otu)
  # print(rownames(otu_table_eliminated)[i])
  write.table(continuousresult, file = paste0("MAGMA_continuous_Data_GLOBAL",name, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
  write.table(binaryresult, file = paste0("MAGMA_binary_Data_GLOBAL",name, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
  
  return (list(otu_table = table_otu, continuous_r = continuousresult, binary_r = binaryresult )  )
  
}

dir.create(file.path(result_path, "LooNet"), showWarnings = FALSE)
dir.create(file.path(result_path, "LooNet_69"), showWarnings = FALSE)




list_not_eliminate = calculating_LOO(otu_table6M, dirs = "LooNet", name = "")
list_not_eliminate_69 = calculating_LOO(otu_table6M_69, dirs = "LooNet_69", name = "_69")



# 81 -------------------------------------------------------


continuousresult = read.table("MAGMA_continuous_Data_GLOBAL.tsv")
binaryresult = read.table("MAGMA_binary_Data_GLOBAL.tsv")
otu_table = as.matrix(read.table("OTU_TABLE_MAGMA6M.tsv"))

otu_table == list_not_eliminate$otu_table



# 69  -----------------------------------------------------------



continuousresult = read.table("MAGMA_continuous_Data_GLOBAL_69.tsv")
binaryresult = read.table("MAGMA_binary_Data_GLOBAL_69.tsv")
otu_table = as.matrix(read.table("OTU_TABLE_MAGMA6M_69.tsv"))




# PHYLUMS -----------------------------------------------------------------

continuousresult = read.table("MAGMA_continuous_Data_GLOBAL.tsv")
binaryresult = read.table("MAGMA_binary_Data_GLOBAL.tsv")
otu_table = as.matrix(read.table("OTU_TABLE_MAGMA6M.tsv"))

setwd(graph_path)

## 

plotting_globalnetwork_grouped = function(otu_table, binaryresult, continuousresult, name){
  plot_graphical_phylum(
    OTU_RESULT = binaryresult,
    label_name_final = paste0("MAGMA_", name,"_Phylum"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "BINARY",
    month = "6M_",
    import_layout =  paste0(LAYOUT_path, "LAYOUT_PHYLUM_6M_graph.txt") #LAYOUT_PHYLUM_6M_graph.txt"
  )# NA #"../6M/LAYOUT_PHYLUM_6M_graph.txt"
  # )
  
  
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Phylum"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_PHYLUM_6M_graph.txt")
  )
  
  # CLAss -------------------------------------------------------------------
  
  
  # level_dataset = data.frame(name = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), rank = seq(1,7))
  plot_graphical_phylum(
    OTU_RESULT = binaryresult,
    label_name_final = paste0("MAGMA_", name,"_Class"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "BINARY",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_Class_6M_graph.txt"), 
    level_taxa_order = "Class"
    
  )
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Class"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_Class_6M_graph.txt"), 
    level_taxa_order = "Class"
    
  )
  
  
  # ORRRDER -----------------------------------------------------------------
  
  
  
  # level_dataset = data.frame(name = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), rank = seq(1,7))
  plot_graphical_phylum(
    OTU_RESULT = binaryresult,
    label_name_final = paste0("MAGMA_", name,"_Order"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "BINARY",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_Order_6M_graph.txt"), 
    level_taxa_order = "Order"
    
  )
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Order"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_Order_6M_graph.txt"), 
    level_taxa_order = "Order"
    
  )
  
  
  
  
  # FAMILY ------------------------------------------------------------------
  
  plot_graphical_phylum(
    OTU_RESULT = binaryresult,
    label_name_final = paste0("MAGMA_", name,"_Family"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "BINARY",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_Family_6M_graph.txt"), 
    level_taxa_order = "Family"
    
  )
  plot_graphical_phylum(
    OTU_RESULT = continuousresult,
    label_name_final = paste0("MAGMA_", name,"_Family"),
    otu_table = otu_table,
    Taxonomy_ASVs = Taxonomy_ASVs,
    continuous_or_binary = "CONTINUOUS",
    month = "6M_",
    import_layout = paste0(LAYOUT_path, "LAYOUT_Family_6M_graph.txt"), 
    level_taxa_order = "Family"
    
  )
  
}
plot_circleplot = function(otu_table, name,Taxonomy_ASVs = Taxonomy_ASVs){
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
  # png("6m_graphs_Paddy_Covariates_MAGMA_5thlevel_FAMILY_circle_plot.png", width = 465, height = 225, units='mm', res = 300)
  # # magma_Stool_AllCov$
  # plot(magma_Stool_AllCov,layout = coords, V.color.factor = df[5,])
  # dev.off()
  
  png(paste0("6M_MAGMA_", name, "_5thlevel_FAMILY_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord, V.color.factor = ifelse(stri_replace_all_regex(df[5,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[5,], "[a-z]__", "")))
  dev.off()
  
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[3,]) ) ) )
  
  png(paste0("6M_MAGMA_", name, "_3rdlevel_CLASS_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord, V.color.factor = ifelse(stri_replace_all_regex(df[3,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[3,], "[a-z]__", "")))
  dev.off()
  
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[2,]) ) ) )
  
  png(paste0("6M_MAGMA_", name, "_2ndlevel_PHYLUM_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord,V.color.factor = ifelse(stri_replace_all_regex(df[2,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[2,], "[a-z]__", "")))
  dev.off()
  
  
  coords_ord <- layout_in_circle(temp_graph, order = order(as.numeric(factor(df[4,]) ) ) )
  
  png(paste0("6M_MAGMA_", name, "_4thlevel_ORDER_circle_plot_order.png"), width = 465, height = 225, units='mm', res = 300)
  # magma_Stool_AllCov$
  plot(magma_Stool_AllCov,layout = coords_ord, V.color.factor = ifelse(stri_replace_all_regex(df[4,], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(df[4,], "[a-z]__", "")))
  dev.off()
}

checkStrict(plotting_globalnetwork_grouped)
checkStrict(plot_circleplot)


# second graph 81  ------------------------------------------------------------


setwd(result_path)
continuousresult = read.table("MAGMA_continuous_Data_GLOBAL.tsv")
binaryresult = read.table("MAGMA_binary_Data_GLOBAL.tsv")
otu_table = as.matrix(read.table("OTU_TABLE_MAGMA6M.tsv"))
setwd(graph_path)

plot_circleplot(otu_table = otu_table, name ="noteliminate", Taxonomy_ASVs = Taxonomy_ASVs)
plotting_globalnetwork_grouped(otu_table = otu_table, binaryresult = binaryresult, continuousresult = continuousresult, name ="")


# 69 -----------------------------------------------------------


setwd(result_path)

continuousresult = read.table("MAGMA_continuous_Data_GLOBAL_69.tsv")
binaryresult = read.table("MAGMA_binary_Data_GLOBAL_69.tsv")
otu_table = as.matrix(read.table("OTU_TABLE_MAGMA6M_69.tsv"))
setwd(graph_path)

plot_circleplot(otu_table = otu_table, name ="69",Taxonomy_ASVs = Taxonomy_ASVs)
plotting_globalnetwork_grouped(otu_table = otu_table, binaryresult = binaryresult, continuousresult = continuousresult, name ="69")


