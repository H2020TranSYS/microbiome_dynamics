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




sixtynine = "_69" # vs ""


# 2 -----------------------------------------------------------------------

sixtynine = "" # vs ""
# eliminated = ""



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
}


edges <- paste(lionessOutput[, 1], lionessOutput[, 2], sep = "_")
rownames(lionessOutput) = edges
lionessOutput[1:5,1:5]



# Write results -----------------------------------------------------------
dir.create(file.path(result_path), showWarnings = TRUE)


setwd(result_path)
# 
fwrite(lionessOutput[,3:ncol(lionessOutput)],file=paste0("Resulting_net_from_corr_MAGMA",sixtynine, ".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames
# 
colnames(Sequences_nodes) = c("Node_list","NAME")
write.table(Sequences_nodes,file = paste0("../Sequences_nodes",sixtynine,".txt"), quote = FALSE,row.names = F, col.names = T)
# 
write.table(files_name_matching,file = paste0("../files_name_matching",sixtynine,".txt"))



# CROP the resulting net --------------------------------------------------


Resulting_net = lionessOutput[,3:ncol(lionessOutput)]
names_r_net = rownames(Resulting_net)
# colnames(Resulting_net) %in% sample_dataset$sid
dim(Resulting_net)
length(names_r_net)

to_eliminate_r = c(rep(NA, 1000000))
j = 0
for (i in 1:(length(names_r_net)-1))
{
  elemento = names_r_net[i]
  fin = unlist(strsplit(elemento,"_"))
  fin2 = paste0(fin[2],sep = "_", fin[1])
  value = match(fin2,names_r_net, nomatch = 0)
  if (value > i)
  {
    j = j +1
    to_eliminate_r[j] =  value
  }
  if (i %% 100 == 0){ print(i)}
}
to_eliminate_r[1:j]
Resulting_net = Resulting_net[-to_eliminate_r[1:j],]
# Resulting_net[aa,]
# getwd()
fwrite(Resulting_net,file=paste0("Resulting_net_from_corr_MAGMA",sixtynine, ".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames



# Since is full of rows with just 0 --> we create a network with j --------


aa = apply(Resulting_net, 1, function(x) sum(abs(x))) > 0
not_null = Resulting_net[aa,]
dim(not_null)
# [1] 547 153
# OR 
# [1] 546 155
fwrite(not_null,file = paste0("Resulting_net_notNULL_MAGMA",sixtynine, ".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames

not_null



# NO node -----------------------------------------------------------------


name_vect = rownames(not_null) #names_r_net
# 
# 
# # sapply(name_vect, function(x) strsplit(x, split = "_")[[1]] )
# 
cc       <- strsplit(name_vect,'_')
part1    <- unlist(cc)[2*(1:length(name_vect))-1]
part2    <- unlist(cc)[2*(1:length(name_vect))  ]
head(part1) ; head(part2) ; head(name_vect)
not_null_nonode = not_null[part1 != part2,]

fwrite(not_null_nonode,file = paste0("Resulting_net_notNULL_MAGMA_nonode_",sixtynine, ".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames
not_null_nonode


