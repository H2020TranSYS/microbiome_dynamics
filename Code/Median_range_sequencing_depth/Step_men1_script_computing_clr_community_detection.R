
# Libraries ---------------------------------------------------------------


library(Rcpp)
library(rMAGMA)
library(dplyr)
library(igraph)
library(stringi)
library(compositions)

library(spoutlier)
library(randcorr)
library(dbscan)
library(pROC)
library(randcorr)
library(DDoutlier)
library(WGCNA)
library(MASS)
library("lionessR")
library("SummarizedExperiment")
# install.packages("lsa")
library(lsa)

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
source("C:/Users/fmelo/Desktop/Backup_Federico/Edge_mod_outlier/LucKi_cohort_analysis/Code/Calculation_function.R")
data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
result_path = "C:/Users/fmelo/Desktop/Backup_Federico/Edge_mod_outlier/LucKi_cohort_analysis/Results/Average_sequencing_depth/"
# graph_path =  "C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Dynamics_in_Microbiome/Graphs/MAGMA_graphs/6M"
# LAYOUT_path = "C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Dynamics_in_Microbiome/Data/Graphical_Layout/" #LAYOUT_PHYLUM_6M_graph.txt"

setwd(data_path)
# setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/6M")

w0_n335_OTU_table6M <- read.delim("SparCC_ASV6M.tab")
w0_n335_OTU_table9M <- read.delim("SparCC_ASV9M.tab")
w0_n335_OTU_table = merge(w0_n335_OTU_table6M, w0_n335_OTU_table9M, by = "X.OTU.ID")

rownames(w0_n335_OTU_table) = w0_n335_OTU_table[,1]
# rownames(w0_n335_OTU_table)


## Covariates
w0_n335_OTU_table = w0_n335_OTU_table %>% dplyr::select(-X.OTU.ID)
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
## From https://www.biorxiv.org/content/10.1101/538579v1.full.pdf
# The structure of real data, we relied on the 16S data of the microbiome of Puerto Rico honey bees
# obtained by MG Dominguez-Bello [58, study ID 1064]. We filtered the data, keeping the 80 OTUs with a
# prevalence greater than 15% and 286 samples with a sequencing depth greater than 100 reads

# https://www.nature.com/articles/s41467-022-28034-z --> on 10% prevalence
# How and whether to conduct independent filtering of data prior to conducting DA tests
# are other important open questions regarding microbiome data analysis7. 
# Although statistical arguments regarding the validity of independent filtering are beyond the scope of this work,
# intuitively it is reasonable to exclude features found in only a small number of samples (regardless of which groups those samples are in).
# The basic reason for this is that otherwise the burden of multiple-test correction becomes so great as to nearly 
# prohibit identifying any differentially abundant features. Despite this drawback, many tools identified large 
# numbers of significant ASVs in the unfiltered data. However, these significant ASVs tended to be more tool-specific
# in the unfiltered data and there was much more variation in the percentage of significant ASVs across tools. 
# Accordingly, we would suggest performing prevalence filtering (e.g., at 10%) of features prior to DA testing, 
# although we acknowledge that more work is needed to estimate an optimal cut-off rather than just arbitrarily selecting one4.

# Before we jump to our analyses, we may want to perform prevalence filtering. Nearing et al. (2021) found that 
# applying a 10% threshold for the prevalence of the taxa generally resulted in more robust results. 
# Some tools have builtin arguments for that. By applying the threshold to our input data, 
# we can make sure it is applied for all tools. Below we show how to do this in mia:


## FROM this, it seems that https://microbiome.github.io/OMA/differential-abundance.html
## clr is after fitlering
irow <- sequencing_depth >500
sum(icol)


df_median_iqr = cbind(median(sequencing_depth),iqr(sequencing_depth))
ind6m = names(sequencing_depth) %in% LuckiMap_6M$Sample.name
df_median_iqr = rbind(df_median_iqr, c(median(sequencing_depth[ind6m]),iqr(sequencing_depth[ind6m])))
df_median_iqr = rbind(df_median_iqr, c(median(sequencing_depth[!ind6m]),iqr(sequencing_depth[!ind6m])))

colnames(df_median_iqr) = c("Median", "IQR")
rownames(df_median_iqr) = c("Global155", "6M", "9M")

min(sequencing_depth)
max(sequencing_depth)
# > min(sequencing_depth)
# [1] 11123
# > max(sequencing_depth)
# [1] 105921
setwd(result_path)
saveRDS(df_median_iqr, file = "Average_sequencing_depth.rds")
# readRDS("Average_sequencing_depth.rds")

