library(dplyr)
library(ggplot2)
library("xlsx")
library(ggpubr)

data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Results/Differential_abundances"
graph_path =  "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/Differential_abundances"


setwd(data_path)

computing_appearing_rel = function(otu_table, Taxonomy_ASVs,  level_taxa_order){
  
  level_dataset = data.frame(name = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), rank = seq(1,7))
  
  otu_table_rel<-as.data.frame(t(apply(otu_table,1, function(x) x/sum(x))))
  rank_current_taxa = level_dataset[level_taxa_order == level_dataset$name, "rank"]
  summed_col = colMeans(otu_table_rel)
  
  ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(otu_table_rel),]
  df = (t(ASV_final))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  level_taxa = unlist(unique(ASV_final %>% dplyr::select (level_taxa_order)))
  
  Taxas_name = colnames(otu_table_rel)
  
  phy_sum = matrix(data = 0, nrow = 1, ncol = length(level_taxa))
  colnames(phy_sum) = level_taxa
  rownames(phy_sum) = c("count")
  phy_cont = phy_sum
  for (i in 1:ncol(otu_table)){
    phy_sum[1,df[rank_current_taxa,Taxas_name[i]] ] = phy_sum[1,df[rank_current_taxa,Taxas_name[i]] ] + summed_col[Taxas_name[i]]
    phy_cont[1,df[rank_current_taxa,Taxas_name[i]] ] = phy_cont[1,df[rank_current_taxa,Taxas_name[i]] ] +1 
    
  }
  taxas_per_phylum = table(df[rank_current_taxa,])
  ordered_taxa_per_phylum = taxas_per_phylum[colnames(phy_cont)]
  return(phy_sum) 
}


typedata = "69"
## TO BE EXTENDED
if (typedata == "81_74"){

OTU_TABLE_MAGMA6M <- read.delim("MAGMA_data/6M/OTU_TABLE_MAGMA6M.tsv")
# 81 x 95
OTU_TABLE_MAGMA9M <- read.delim("MAGMA_data/9M/OTU_TABLE_MAGMA9M.tsv")
# 74 x 95

} else if (typedata == "69") {
  OTU_TABLE_MAGMA6M <- read.delim("MAGMA_data/6M/OTU_TABLE_MAGMA6M_69.tsv")
  # 69 x 95
  OTU_TABLE_MAGMA9M <- read.delim("MAGMA_data/9M/OTU_TABLE_MAGMA9M_69.tsv")
  # 69 x 95
  
}

Taxonomy_ASVs <- read.csv("Taxonomy_ASVs.csv")




otu_table = OTU_TABLE_MAGMA6M

level_taxa_order = "Phylum"

ooo = computing_appearing_rel(OTU_TABLE_MAGMA6M, Taxonomy_ASVs, level_taxa_order)
ooow = computing_appearing_rel(OTU_TABLE_MAGMA9M, Taxonomy_ASVs, level_taxa_order)

dfdata = data.frame("condition" = gsub("p__","",colnames(ooo)), value = as.numeric(ooo ), "specie" = "6Months")
dfdataw = data.frame("condition" = gsub("p__","",colnames(ooow)), value = as.numeric(ooow ),"specie" = "9Months")


setwd(result_path)

Merged_dataset = merge(dfdata, dfdataw, by = "condition")
colnames(Merged_dataset) = c("Condition", "Absolute_abundance.x", "Time.x", "Absolute_abundance.y", "Time.y")
Merged_dataset$Relative_abudance.x = round(Merged_dataset$Absolute_abundance.x / sum(Merged_dataset$Absolute_abundance.x) * 100, 3)
Merged_dataset$Relative_abudance.y = round(Merged_dataset$Absolute_abundance.y / sum(Merged_dataset$Absolute_abundance.y) * 100, 3)
Merged_dataset %>% dplyr::select(c(Condition, Absolute_abundance.x, Relative_abudance.x, Time.x, Absolute_abundance.y, Relative_abudance.y, Time.y   ))  %>%
 write.xlsx(.,file = paste0("Frac_Abundances_BY_IND_", typedata, ".xlsx"),
            sheetName = level_taxa_order, append = FALSE)


dfdata = rbind(dfdata, dfdataw)

# library


# create a dataset
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# data <- data.frame(specie,condition,value)

# Stacked + percent

setwd(graph_path)
g1 =ggplot(dfdata, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")+
  labs(y= "Fractional abundances", x = "Timepoints",fill='') +
  ggtitle(level_taxa_order)+
  theme_classic(base_size = 25, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=25))

png(paste0(typedata,"_IND_FRAC_Phylum_difference.png"), width = 300, height = 300, units='mm', res = 300)
g1
dev.off()



# Class ------------------------------------------------------------------


level_taxa_order =  "Class"# , "Order", "Family"

ooo = computing_appearing_rel(OTU_TABLE_MAGMA6M, Taxonomy_ASVs, level_taxa_order)
ooow = computing_appearing_rel(OTU_TABLE_MAGMA9M, Taxonomy_ASVs, level_taxa_order)

dfdata = data.frame("condition" = gsub("c__","",colnames(ooo)), value = as.numeric(ooo ), "specie" = "6Months")
dfdataw = data.frame("condition" = gsub("c__","",colnames(ooow)), value = as.numeric(ooow ),"specie" = "9Months")

setwd(result_path)

Merged_dataset = merge(dfdata, dfdataw, by = "condition")
colnames(Merged_dataset) = c("Condition", "Absolute_abundance.x", "Time.x", "Absolute_abundance.y", "Time.y")
Merged_dataset$Relative_abudance.x = round(Merged_dataset$Absolute_abundance.x / sum(Merged_dataset$Absolute_abundance.x) * 100, 3)
Merged_dataset$Relative_abudance.y = round(Merged_dataset$Absolute_abundance.y / sum(Merged_dataset$Absolute_abundance.y) * 100, 3)
Merged_dataset %>% dplyr::select(c(Condition, Absolute_abundance.x, Relative_abudance.x, Time.x, Absolute_abundance.y, Relative_abudance.y, Time.y   ))  %>%
  write.xlsx(.,file = paste0("Frac_Abundances_BY_IND_", typedata, ".xlsx"),
             sheetName = level_taxa_order, append = TRUE)


dfdata = rbind(dfdata, dfdataw)

g2 = ggplot(dfdata, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")+
  labs(y= "Fractional abundances", x = "Timepoints",fill='') +
  ggtitle(level_taxa_order)+
  theme_classic(base_size = 25, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=25))

setwd(graph_path)

png(paste0(typedata,"_IND_FRAC_Class_difference.png"), width = 300, height = 300, units='mm', res = 300)
g2
dev.off()



# Order -------------------------------------------------------------------



level_taxa_order =   "Order"#, "Family"


ooo = computing_appearing_rel(OTU_TABLE_MAGMA6M, Taxonomy_ASVs, level_taxa_order)
ooow = computing_appearing_rel(OTU_TABLE_MAGMA9M, Taxonomy_ASVs, level_taxa_order)

dfdata = data.frame("condition" = gsub("o__","",colnames(ooo)), value = as.numeric(ooo ), "specie" = "6Months")
dfdataw = data.frame("condition" = gsub("o__","",colnames(ooow)), value = as.numeric(ooow ),"specie" = "9Months")
dfdata[dfdata == ""] = "Unknown"
dfdataw[dfdataw == ""] = "Unknown"

setwd(result_path)

Merged_dataset = merge(dfdata, dfdataw, by = "condition")
colnames(Merged_dataset) = c("Condition", "Absolute_abundance.x", "Time.x", "Absolute_abundance.y", "Time.y")
Merged_dataset$Relative_abudance.x = round(Merged_dataset$Absolute_abundance.x / sum(Merged_dataset$Absolute_abundance.x) * 100, 3)
Merged_dataset$Relative_abudance.y = round(Merged_dataset$Absolute_abundance.y / sum(Merged_dataset$Absolute_abundance.y) * 100, 3)
Merged_dataset %>% dplyr::select(c(Condition, Absolute_abundance.x, Relative_abudance.x, Time.x, Absolute_abundance.y, Relative_abudance.y, Time.y   ))  %>%
  write.xlsx(.,file = paste0("Frac_Abundances_BY_IND_", typedata, ".xlsx"),
             sheetName = level_taxa_order, append = TRUE)



dfdata = rbind(dfdata, dfdataw)
g3 = ggplot(dfdata, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")+
  labs(y= "Fractional abundances", x = "Timepoints",fill='') +
  ggtitle(level_taxa_order)+
  theme_classic(base_size = 25, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=25))

setwd(graph_path)
png(paste0(typedata,"_IND_FRAC_Order_difference.png"), width = 300, height = 300, units='mm', res = 300)
g3
dev.off()



# Family -------------------------------------------------------------------



level_taxa_order =    "Family"


ooo = computing_appearing_rel(OTU_TABLE_MAGMA6M, Taxonomy_ASVs, level_taxa_order)
ooow = computing_appearing_rel(OTU_TABLE_MAGMA9M, Taxonomy_ASVs, level_taxa_order)

dfdata = data.frame("condition" = gsub("f__","",colnames(ooo)), value = as.numeric(ooo ), "specie" = "6Months")
dfdataw = data.frame("condition" = gsub("f__","",colnames(ooow)), value = as.numeric(ooow ),"specie" = "9Months")
dfdata[dfdata == ""] = "Unknown"
dfdataw[dfdataw == ""] = "Unknown"

setwd(result_path)

Merged_dataset = merge(dfdata, dfdataw, by = "condition")
colnames(Merged_dataset) = c("Condition", "Absolute_abundance.x", "Time.x", "Absolute_abundance.y", "Time.y")
Merged_dataset$Relative_abudance.x = round(Merged_dataset$Absolute_abundance.x / sum(Merged_dataset$Absolute_abundance.x) * 100, 3)
Merged_dataset$Relative_abudance.y = round(Merged_dataset$Absolute_abundance.y / sum(Merged_dataset$Absolute_abundance.y) * 100, 3)
Merged_dataset %>% dplyr::select(c(Condition, Absolute_abundance.x, Relative_abudance.x, Time.x, Absolute_abundance.y, Relative_abudance.y, Time.y   ))  %>%
  write.xlsx(.,file = paste0("Frac_Abundances_BY_IND_", typedata, ".xlsx"),
             sheetName = level_taxa_order, append = TRUE)


dfdata = rbind(dfdata, dfdataw)
g4 = ggplot(dfdata, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")+
  labs(y= "Fractional abundances", x = "Timepoints",fill='') +
  ggtitle(level_taxa_order)+
  theme_classic(base_size = 25, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=15))


setwd(graph_path)
png(paste0(typedata,"_IND_FRAC_Family_difference.png"), width = 300, height = 300, units='mm', res = 300)
g4
dev.off()

png(paste0(typedata,"_IND_FRAC_Family_difference_.png"), width = 350, height = 300, units='mm', res = 300)
g4
dev.off()


png(paste0(typedata,"_IND_FRAC_glob_difference.png"), width = 550, height = 350, units='mm', res = 300)
ggarrange(g1, g2,g3,g4, ncol = 2, nrow = 2)
dev.off()

