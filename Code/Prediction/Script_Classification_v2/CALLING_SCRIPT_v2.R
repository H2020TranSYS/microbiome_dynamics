

# Parameters --------------------------------------------------------------

# eliminate = T
# diagnosi = c("CD", "UC") ## diagnosi = c("CD","UC")
# terapia  =c("VDZ")# c("VDZ") # terapia = c("TNF", "UST", "VDZ")
# setwd("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/SPARCC_and_MAGMA/MAGMA_result_CORRECT_CONF_preprocessing_per_group//")


# setted_wd = args$setted_wd
# reading_file = args$reading_file
# rankfeature = args$rankfeature


# rm(args)
rm(list = ls())
args = list()

library(dplyr)
library(igraph)
library(ROCR)
library(e1071)

setwd("/massstorage/URT/GEN/BIO3/Federico/LucKi_project/MAGMA_CLASSIFICATION")
files <- list.files(path="Data/.", pattern="*", 
                    full.names=TRUE, recursive=FALSE)

files_m = gsub("Data/./","", files)
# "LuckiMap_9M.txt" %in% files_m
to_rem = c(match("LuckiMap_9M.txt", files_m), match("LuckiMap_6M.txt", files_m) )
final_file = files_m[-to_rem]

final_file[3]

curr_wd = getwd()

ex_data = read.delim(file="Data/OTU_TABLE_69_6M.tsv" )
taxa_name = colnames(ex_data)
# Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])

# substr(a,(nchar(a)-3), 100000)
# final_file2 = final_file[grepl("EDGEW", final_file, fixed = TRUE)]

for (file in final_file){
  setwd(curr_wd)
  extens = substr(file,(nchar(file)-3), 100000)
  type_data = ifelse(extens == ".txt","Edge",ifelse(extens == ".rds", "EDNN", "NODE"))
  # if(type_data != "Edge"){next;}
if (grepl("^dist", file)){
print(file)
next; 
}
  rm(args)
  args = list()
  args$setted_wd = curr_wd
  args$type_data = type_data
  args$reading_file = file
  args$rep = 2
  args$colnm = taxa_name
  args$rep = 100
  args$rankfeature = "internal"
  args$ordering = F
  source ("Script_Classification_v2/script_Delivery_type.R")
  setwd(curr_wd)
  source ("Script_Classification_v2/script_PERSISTENT_vs_NONPERSISTENT_DIET_6m_9m.R")
  setwd(curr_wd)

  # args$k_csection =  5 
  # args$k_vaginal = 29
  
  # source ("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/CLASSIFICATION/Script_Classification_v2/script_SVM_LKO.R")
  # source ("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/CLASSIFICATION/Script_Classification_v2/script_DIET_4_outcome_transitions.R")
  
  print("V2")
  
  args$ordering = T
  if (args$type_data != "EDNN"){next; }
  source ("Script_Classification_v2/script_Delivery_type.R")
  setwd(curr_wd)


  source ("Script_Classification_v2/script_PERSISTENT_vs_NONPERSISTENT_DIET_6m_9m.R")
  setwd(curr_wd)


  
}
