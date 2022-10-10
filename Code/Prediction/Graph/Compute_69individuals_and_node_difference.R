library(compositions)
## library for the clr trasformation

setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/OTUS_node_analysis_difference")

## Importing base data for no trasf and filtering + clr


OTU_TABLE_MAGMA6M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/6M/OTU_TABLE_MAGMA6M.tsv")
OTU_TABLE_MAGMA9M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/9M/OTU_TABLE_MAGMA9M.tsv")

OTU_TABLE_MAGMA6M ==OTU_TABLE_MAGMA9M
## we want them to have different dimension



# TRASFORMATION -----------------------------------------------------------



class(OTU_TABLE_MAGMA6M)

## trasformed table
OTU_TABLE_MAGMA_CLR_6M = data.frame( clr( OTU_TABLE_MAGMA6M ))
OTU_TABLE_MAGMA_CLR_9M = data.frame( clr( OTU_TABLE_MAGMA9M ))
## we trasfromed back to a data frame, it was a weird clr object

## compositional trasformation, if I sum over an individual, the sum = 0 
apply(OTU_TABLE_MAGMA_CLR_9M,1,sum) ; apply(OTU_TABLE_MAGMA_CLR_6M,1,sum)


## maybe better to correct BEFORE the filtering?

## MAPPING FILES

LuckiMap_9M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Correlations_Lucki_for_Kristel/9 Months/LuckiMap_9M.txt")
# C   View(LuckiMap_9M)
LuckiMap_6M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Correlations_Lucki_for_Kristel/6 Months/LuckiMap_6M.txt")
# >   View(LuckiMap_6M)


Children =  merge(LuckiMap_6M, LuckiMap_9M, by = "Child")

## double road: the filtered NOT trasformed and the CLR ones 

OTU_TABLE_FILT_6M = OTU_TABLE_MAGMA6M[rownames(OTU_TABLE_MAGMA6M) %in% Children$Sample.name.x,]
OTU_TABLE_FILT_9M = OTU_TABLE_MAGMA9M[rownames(OTU_TABLE_MAGMA9M) %in% Children$Sample.name.y,]
## filtering on the observtion in the 6m or 9m

OTU_TABLE_FILT_CLR_6M = OTU_TABLE_MAGMA_CLR_6M[rownames(OTU_TABLE_MAGMA_CLR_6M) %in% Children$Sample.name.x,]
OTU_TABLE_FILT_CLR_9M = OTU_TABLE_MAGMA_CLR_9M[rownames(OTU_TABLE_MAGMA_CLR_9M) %in% Children$Sample.name.y,]


colnames(OTU_TABLE_FILT_6M) ==  colnames(OTU_TABLE_FILT_9M) 
# sanity check 

Children$Child[match(Children$Sample.name.x,rownames(OTU_TABLE_FILT_6M))]== Children$Child[match(Children$Sample.name.y,rownames(OTU_TABLE_FILT_9M))]
## they are not in the same order, we have to order them 


OT_ORD_6M = OTU_TABLE_FILT_6M[match(Children$Sample.name.x,rownames(OTU_TABLE_FILT_6M)),]
OT_ORD_9M = OTU_TABLE_FILT_9M[match(Children$Sample.name.y,rownames(OTU_TABLE_FILT_9M)),]
# ordering it with the same order as in cildren 

OT_ORD_CLR_6M = OTU_TABLE_FILT_CLR_6M[match(Children$Sample.name.x,rownames(OTU_TABLE_FILT_CLR_6M)),]
OT_ORD_CLR_9M = OTU_TABLE_FILT_CLR_9M[match(Children$Sample.name.y,rownames(OTU_TABLE_FILT_CLR_9M)),]


rownames(OT_ORD_6M)[1:5]; rownames(OT_ORD_9M)[1:5]
# Children[1:5,c("Sample.name.x", "Sample.name.y", "Child")]
# match(rownames(OT_ORD_6M),Children$Sample.name.x) == match(rownames(OT_ORD_9M),Children$Sample.name.y)

rownames(OT_ORD_6M) <- Children$Child[match(rownames(OT_ORD_6M), Children$Sample.name.x)]
rownames(OT_ORD_9M) <- Children$Child[match(rownames(OT_ORD_9M), Children$Sample.name.y)]
rownames(OT_ORD_6M) == rownames(OT_ORD_9M)


rownames(OT_ORD_CLR_6M) <- Children$Child[match(rownames(OT_ORD_CLR_6M), Children$Sample.name.x)]
rownames(OT_ORD_CLR_9M) <- Children$Child[match(rownames(OT_ORD_CLR_9M), Children$Sample.name.y)]
rownames(OT_ORD_CLR_6M) == rownames(OT_ORD_CLR_9M)



diff_OTUS = OT_ORD_9M - OT_ORD_6M

diff_OTUS_CLR = OT_ORD_CLR_9M - OT_ORD_CLR_6M



setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/OTUS_node_analysis_difference")

write.table(OT_ORD_6M, file = paste0("OTU_TABLE_69_SUB_6M", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
write.table(OT_ORD_9M, file = paste0("OTU_TABLE_69_SUB_9M", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
write.table(diff_OTUS, file = paste0("OTU_TABLE_69_DIFFERENCE", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")

write.table(OT_ORD_CLR_6M, file = paste0("OTU_TABLE_69_SUB_CLR_6M", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
write.table(OT_ORD_CLR_9M, file = paste0("OTU_TABLE_69_SUB_CLR_9M", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")
write.table(diff_OTUS_CLR, file = paste0("OTU_TABLE_69_CLR_DIFFERENCE", ".tsv" ), quote = FALSE, row.names = TRUE, col.names = TRUE,sep = "\t")













# CREATE CLR FILTERED STARTING FILE ---------------------------------------


# Libraries ---------------------------------------------------------------


library(Rcpp)
library(rMAGMA)
library(dplyr)
library(igraph)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# name <- read.table("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/name.txt", quote="\"", comment.char="")
# individuals <- read.table("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/Data_for_SPARCC/individuals.txt", quote="\"", comment.char="")
# Ind_list = cbind("Individual_ID"= name,"NAME"= paste0(rep("Ind_", nrow(column_Individual)),individuals$V1))
# colnames(Ind_list) = c("Individual_ID", "NAME")
# write.table(Ind_list,file = "mapping_nameIndividual.txt",quote = FALSE, row.names = FALSE, col.names = TRUE)


# Read Paddy data ---------------------------------------------------------


setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/OTUS_node_analysis_difference/")

w0_n335_OTU_table6M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/SparCC_ASV6M.tab")
w0_n335_OTU_table9M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/SparCC_ASV9M.tab")
min(w0_n335_OTU_table9M[,-1])
max(w0_n335_OTU_table6M[,-1])

w0_n335_OTU_table = merge(w0_n335_OTU_table6M, w0_n335_OTU_table9M, by = "X.OTU.ID")
# w0_n335_OTU_table %>% select(-c("")) = merge(w0_n335_OTU_table6M, w0_n335_OTU_table9M, by = "X.OTU.ID")

rownames(w0_n335_OTU_table) = w0_n335_OTU_table[,1]; w0_n335_OTU_table =w0_n335_OTU_table %>% select(-X.OTU.ID)
rownames(w0_n335_OTU_table6M) = w0_n335_OTU_table6M[,1]; w0_n335_OTU_table6M = w0_n335_OTU_table6M %>% select(-X.OTU.ID)
rownames(w0_n335_OTU_table9M) = w0_n335_OTU_table9M[,1]; w0_n335_OTU_table9M = w0_n335_OTU_table9M %>% select(-X.OTU.ID)

rownames(w0_n335_OTU_table9M) == rownames(w0_n335_OTU_table)
matched_merged_6m = match(rownames(w0_n335_OTU_table),rownames(w0_n335_OTU_table6M))
w0_n335_OTU_table6M = w0_n335_OTU_table6M[matched_merged_6m,]
matched_merged_9m = match(rownames(w0_n335_OTU_table),rownames(w0_n335_OTU_table9M))
w0_n335_OTU_table9M = w0_n335_OTU_table9M[matched_merged_6m,]


rownames(w0_n335_OTU_table6M) == rownames(w0_n335_OTU_table)

colnames(w0_n335_OTU_table) %in% colnames(w0_n335_OTU_table9M)

# reschaled cause in merge the order is changed
# w0_n335_OTU_table = w0_n335_OTU_table[matched_merged,]
rownames(w0_n335_OTU_table9M) == rownames(w0_n335_OTU_table)



t_otu = t(w0_n335_OTU_table)
t_otu6M = t(w0_n335_OTU_table6M)
t_otu9M = t(w0_n335_OTU_table9M)

max(t_otu9M)

## to have rows as individual amnd columns as taxa
prevalence       <- colMeans(t_otu>0)
sequencing_depth <- rowSums(t_otu)

icol <- prevalence >0.15 ##the prevalence value decreased to assure more OTUS
irow <- sequencing_depth >500
sum(icol)


sum(irow) == length(irow) ## means that all the individuals are selected

otu_table <- t_otu[irow,icol] ## w0_n335_OTU_table[irow,icol]
otu_table6M <- t_otu6M[,icol] ## w0_n335_OTU_table[irow,icol]
otu_table9M <- t_otu9M[,icol] ## w0_n335_OTU_table[irow,icol]


rownames(w0_n335_OTU_table9M) == rownames(w0_n335_OTU_table)
dim(otu_table) #[1] 155  95
# write.table(otu_table, file = "OTU_TABLE_MAGMA.tsv",quote = F,sep = "\t",row.names = T, col.names = T)
otu_table[1:5,1:5]

# OTU TABLE 9M ------------------------------------------------------------
colnames(otu_table)

w0_OTU_TABLE_MAGMA_CLR_6M = data.frame( clr( t_otu6M ))
w0_OTU_TABLE_MAGMA_CLR_9M = data.frame( clr( t_otu9M ))
max(w0_OTU_TABLE_MAGMA_CLR_9M)
min(w0_OTU_TABLE_MAGMA_CLR_9M)

filt_w0_OTU_TABLE_MAGMA_CLR_6M = w0_OTU_TABLE_MAGMA_CLR_6M[,icol]
filt_w0_OTU_TABLE_MAGMA_CLR_9M = w0_OTU_TABLE_MAGMA_CLR_9M[,icol]
min(filt_w0_OTU_TABLE_MAGMA_CLR_6M)
max(filt_w0_OTU_TABLE_MAGMA_CLR_9M)
min(filt_w0_OTU_TABLE_MAGMA_CLR_9M)

dim(filt_w0_OTU_TABLE_MAGMA_CLR_6M)
dim(otu_table)


colnames(filt_w0_OTU_TABLE_MAGMA_CLR_9M) == colnames(otu_table)
write.table(filt_w0_OTU_TABLE_MAGMA_CLR_9M, file = "OTU_TABLE_MAGMA_CLR_BEFORE_FILT_9M_74obs.tsv",quote = F,sep = "\t",row.names = T, col.names = T)
write.table(filt_w0_OTU_TABLE_MAGMA_CLR_6M, file = "OTU_TABLE_MAGMA_CLR_BEFORE_FILT_6M_81obs.tsv",quote = F,sep = "\t",row.names = T, col.names = T)





OTU_TABLE_FILT_CLR_6M_before_filt = filt_w0_OTU_TABLE_MAGMA_CLR_6M[rownames(filt_w0_OTU_TABLE_MAGMA_CLR_6M) %in% Children$Sample.name.x,]
OTU_TABLE_FILT_CLR_9M_before_filt = filt_w0_OTU_TABLE_MAGMA_CLR_9M[rownames(filt_w0_OTU_TABLE_MAGMA_CLR_9M) %in% Children$Sample.name.y,]
OT_ORD_CLR_6M_before_filt = OTU_TABLE_FILT_CLR_6M_before_filt[match(Children$Sample.name.x,rownames(OTU_TABLE_FILT_CLR_6M_before_filt)),]
OT_ORD_CLR_9M_before_filt = OTU_TABLE_FILT_CLR_9M_before_filt[match(Children$Sample.name.y,rownames(OTU_TABLE_FILT_CLR_9M_before_filt)),]

rownames(OT_ORD_CLR_6M_before_filt) <- Children$Child[match(rownames(OT_ORD_CLR_6M_before_filt), Children$Sample.name.x)]
rownames(OT_ORD_CLR_9M_before_filt) <- Children$Child[match(rownames(OT_ORD_CLR_9M_before_filt), Children$Sample.name.y)]
rownames(OT_ORD_CLR_6M_before_filt) == rownames(OT_ORD_CLR_9M_before_filt)


diff_OTU_CLR_before_filt = OT_ORD_CLR_9M_before_filt - OT_ORD_CLR_6M_before_filt


dim(OT_ORD_CLR_9M_before_filt)
dim(otu_table)

colnames(OT_ORD_CLR_9M_before_filt) == colnames(otu_table)
write.table(OT_ORD_CLR_9M_before_filt, file = "OTU_TABLE_MAGMA_CLR_BEFORE_FILT_9M_69obs.tsv",quote = F,sep = "\t",row.names = T, col.names = T)
write.table(OT_ORD_CLR_6M_before_filt, file = "OTU_TABLE_MAGMA_CLR_BEFORE_FILT_6M_69obs.tsv",quote = F,sep = "\t",row.names = T, col.names = T)
write.table(diff_OTU_CLR_before_filt, file = "DIFF_OTU_TABLE_MAGMA_CLR_BEFORE_FILT_69obs.tsv",quote = F,sep = "\t",row.names = T, col.names = T)



















