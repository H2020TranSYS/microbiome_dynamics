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

## ITERATE on those 2

# 1 -----------------------------------------------------------------------


sixtynine = "_69" # vs ""

# 2 -----------------------------------------------------------------------

sixtynine = "" # vs ""


data_path = paste0("C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/9M/LooNet", 
                   sixtynine, "/")
result_path = paste0("C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/9M/ISN", 
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
##446 obs of 69


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

# 
# first <- rep(1:122, ncol(global_net))
# second <- rep(1:122, each = nrow(global_net))
# paste(first, second, sep = "_")
# # TRIALS ------------------------------------------------------------------
# name_vect = names_r_net
# 
# 
# # sapply(name_vect, function(x) strsplit(x, split = "_")[[1]] )
# 
# cc       <- strsplit(name_vect,'_')
# part1    <- unlist(cc)[2*(1:length(name_vect))-1]
# part2    <- unlist(cc)[2*(1:length(name_vect))  ]
# head(part1) ; head(part2) ; head(name_vect)
# 
# 
# 
# 
# # cutting -----------------------------------------------------------------
# 
# Final_net = cbind(first, second,part1,part2, Resulting_net[,2:ncol(Resulting_net)])
# 
# Final_net2 = Final_net[Final_net$first > Final_net$second,]
# Final_net2
# final_namenodes = paste(Final_net2$part1, Final_net2$part2,sep = "_")
# Final_net2 = Final_net2[,3:ncol(Final_net2)]
# First_Ind_spec_net = Final_net2[,1:3]

# getwd()
# fwrite(Final_net2,file="../Resulting_net_from_corr_eliminate.txt",
#        sep = " ", row.names = TRUE ) # keeps the rownames




## BACKUP

# TRIAL ALL THE DEGREE CALCULATION ARE CORRECT ----------------------------




# 
# library("GoFKernel")
# library(igraph)
# 
# average_degree = matrix(NA, ncol = 1144,nrow = ncol(Final_net2)-2)
# average_degree_cumulative =  matrix(NA, ncol = 1144,nrow = ncol(Final_net2)-2)
# between_degree_cumulative =  matrix(NA, ncol = 1144,nrow = ncol(Final_net2)-2)
# Louv_list = list()
# for (i in 3:ncol(Final_net2)){
#   print(i)
#   rete = i - 2 
#   First_Ind_spec_net = subset(Final_net2, select = c(1, 2, i))
#   
#   First_Ind_spec_abs = First_Ind_spec_net[,c(1,2,3)]
#   First_Ind_spec_abs[,3] = abs(First_Ind_spec_abs[,3])
#   # First_Ind_spec_abs$LUC112_9m = abs(First_Ind_spec_abs$LUC112_9m)
#   colnames(First_Ind_spec_abs)[3] = "weight"
#   krack_full <- graph.data.frame(First_Ind_spec_abs, directed = FALSE) 
#   is.directed(krack_full)
#   is.weighted(krack_full)
#   
#   ss = strength(krack_full)
#   average_degree[rete,] = ss
#   colnames(average_degree) = names(ss)
#   
# }
# as.data.frame(average_degree)
# rownames(average_degree) = colnames(Final_net2)[3:ncol(Final_net2)]
# 
# 
# distr_6m = apply(average_degree,1,mean)
# 
# 
# # comp
# LuckiMap_9M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Correlations_Lucki_for_Kristel/9 Months/LuckiMap_9M.txt")
# # C   View(LuckiMap_9M)
# LuckiMap_6M <- read.delim("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Correlations_Lucki_for_Kristel/6 Months/LuckiMap_6M.txt")
# # >   View(LuckiMap_6M)
# 
# Children =  merge(LuckiMap_6M, LuckiMap_9M, by = "Child")
# 
# Children$Sample.name.x = paste(Children$Sample.name.x,"6m",sep = "_") 
# Children$Sample.name.y = paste(Children$Sample.name.y,"9m",sep = "_")
# 
# rr = rownames(average_degree)
# 
# average_filter_degree = average_degree[rownames(average_degree) %in% Children$Sample.name.x  ,]
# dim(average_filter_degree)
# rr
# rownames(average_filter_degree)
# 
# dim(average_degree)
# 
# fwrite(average_degree, file = "../../Result_statistic/6month/average_degree_raw.txt",sep = " ", row.names = TRUE )
# 
# fwrite(average_degree_cumulative, file = "../../Result_statistic/6month/average_degree_cdf.txt",sep = " ", row.names = TRUE )
# 
# fwrite(between_degree_cumulative, file = "../../Result_statistic/6month/abetween_degree_cdf.txt",sep = " ", row.names = TRUE )
# 
# # fwrite(Louv_list, file = "../../Result_statistic/9month/abetween_degree_cdf.txt",sep = " ", row.names = TRUE )
# 
# library(yaml)
# yaml::write_yaml(Louv_list, "../../Result_statistic/6month/listLouv.yaml")
# 
# Louv_list2 = yaml::read_yaml("../../Result_statistic/6month/listLouv.yaml")



