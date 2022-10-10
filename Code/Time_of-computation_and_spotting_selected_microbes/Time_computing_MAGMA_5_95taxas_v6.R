
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
plot(magma_Stool_NoCov)


# OTU TABLE 6M ------------------------------------------------------------
colnames(otu_table)
Children =  merge(LuckiMap_6M, LuckiMap_9M, by = "Child")


# MATCHING ----------------------------------------------------------------


# select paired samples  ----------------------------------------------------------

## ALL THE SAMPLES,
otu_table6M = otu_table[rownames(otu_table) %in% LuckiMap_6M$Sample.name,]
sum(rownames(otu_table6M) %in% LuckiMap_6M$Sample.name)
# 81 ; 69 shared


magma_Stool_NoCov <- magma(data = otu_table6M)


# Example with the data of the library ------------------------------------


## CREATE LOO NETWORKs

# mapping_nameIndividual_SPARCC <- read.csv("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/mapping_nameIndividual_SPARCC.txt", sep="")
calculating_LOO_time = function(table_otu, dirs, name){
  start.time <- Sys.time()
  
  for (i in 1: 5){
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
    # write.table(continuousresult, file = paste0(dirs, "/MAGMA_continuous_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,
    #             sep = "\t")
    # write.table(binaryresult, file = paste0(dirs, "/MAGMA_binary_Data",Individual, ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,
    #             sep = "\t")
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  return(time.taken)
}


lists = calculating_LOO_time(otu_table6M, dirs = "LooNet", name = "")
dim(otu_table6M)
# > dim(otu_table6M)
# [1] 81 95
# > lists
# Time difference of 55.14954 secs

