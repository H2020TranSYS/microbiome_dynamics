# # THE SAME BUT WITH CORRELATION -------------------------------------------
# library(data.table)
# 
# 
# w0_n335_OTU_table <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/w0_n335_OTU_table.txt")
# rownames(w0_n335_OTU_table)
# rownames(w0_n335_OTU_table) = w0_n335_OTU_table[,1]
# 
# 
# ## Covariates
# w0_n335_OTU_table = w0_n335_OTU_table %>% select(-SampleList_Metadata_Prediction)
# w0_n335_metadata <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Paddy_work_microbiome/w0_n335_metadata.txt")
# first_covariates_selected = w0_n335_metadata %>% select(c(Age, Gender, BMI,FC.nummer))
# sum(is.na(first_covariates_selected)) ; head(first_covariates_selected)
# first_covariates_selected$Gender = ifelse(first_covariates_selected$Gender == "F",1,0)
# first_covariates_selected$Age = as.numeric(gsub(",",".",first_covariates_selected$Age))
# first_covariates_selected$BMI = as.numeric(gsub(",",".",first_covariates_selected$BMI))
# 
rm(list = ls())
library(data.table)

# before it was on the non correct ones 




# 1 -----------------------------------------------------------------------

sixtynine = "_69" # vs ""

# 2 -----------------------------------------------------------------------

sixtynine = "" # vs ""


data_path = paste0("C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/6M/LooNet", 
                   sixtynine, "/")
result_path = paste0("C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/6M/ISN", 
                     sixtynine, "/")
setwd(data_path)
getwd()


# FILES matching ----------------------------------------------------------

files <- list.files(path=".", pattern="MAGMA_cont*", 
                    full.names=TRUE, recursive=FALSE)
names_file = gsub("./", "",files)
names_file = gsub(".tsv","",names_file)
names_file = gsub("MAGMA_continuous_Data","", names_file)


files_name_matching = cbind(files, names_file)

global_net = read.table(file= paste0("../MAGMA_continuous_Data_GLOBAL", sixtynine,".tsv"), sep = '\t', header = TRUE, row.names = 1)
Node_list = paste0("Node", seq(1,length(rownames(global_net))))
Sequences_nodes = cbind(Node_list, rownames(global_net))
sum(colnames(global_net) != Sequences_nodes[,2])

# CHECK
try(if(sum(rownames(global_net) != Sequences_nodes[,2]) > 0 ) stop("NOT MATCHINGs"))

colnames(global_net) = rownames(global_net) = Node_list


# FLATTEN THE GLOBAL net
global_net_vect <- c(as.matrix(global_net))
lionessOutput <- matrix(NA, nrow(global_net) * ncol(global_net), length(files) + 
                          2)
samples = names_file
colnames(lionessOutput) <- c("reg", "tar", samples)
lionessOutput[, 1] <- rep(row.names(global_net), ncol(global_net))
lionessOutput[, 2] <- rep(colnames(global_net), each = nrow(global_net))
lionessOutput <- as.data.frame(lionessOutput, stringsAsFactors = FALSE)




# Calculation ISNs and average strenght -----------------------------------
start.time <- Sys.time()


for (i in 1: length(files)){
  # for (i in 1: 10){
  recons_net = read.table(file= list.files(path=".", pattern="MAGMA_cont*", 
                                           full.names=TRUE, recursive=FALSE)[i], sep = '\t', header = TRUE, row.names = 1)
  colnames(recons_net) = rownames(recons_net)
  dim(recons_net)
  try(if(sum(rownames(recons_net) != Sequences_nodes[,2]) > 0 ) stop("NOT MATCHINGs"))
  colnames(recons_net) = rownames(recons_net) = Node_list 
  corr_recons_net = c(as.matrix(recons_net)) #  scales a covariance matrix into the corresponding correlation matrix efficiently.
  head(recons_net[1:5,1:5])
  
  
  lionessOutput[, i + 2] <- length(files) * (global_net_vect - corr_recons_net) + corr_recons_net
  print(paste0(i, " = " , str(mean(abs( lionessOutput[, i + 2] )))))
  if (i == 5){end.time <- Sys.time()
}
}

time.taken <- end.time - start.time
time.taken

