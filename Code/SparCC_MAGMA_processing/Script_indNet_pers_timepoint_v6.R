# THE SAME BUT WITH CORRELATION -------------------------------------------
library(data.table)


data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/SparCC_data/"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/SparCC_data/ISN"




## RUN twice for 6M and 9M

timepoint = "9M"
timepoint = "6M"

setwd(data_path)


# FILES matching ----------------------------------------------------------

files <- list.files(path=paste0("./LooNet",timepoint,"/"), pattern="median_correlation*", 
                    full.names=FALSE, recursive=FALSE)
# names_file = gsub("./", "",files)
names_file = gsub(".tsv","",files)
names_file = gsub("median_correlation_","", names_file)

files_name_matching = cbind(files, names_file)

global_net = read.table(file= paste0("GlobNet",timepoint,"/median_correlation.tsv"), sep = '\t', header = FALSE, row.names = 1)
Node_list = paste0("Node", seq(1,length(rownames(global_net))))


Sequences_nodes = cbind(Node_list, rownames(global_net))
sum(rownames(global_net) != Sequences_nodes[,2])


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
  recons_net = read.table(file= paste0("LooNet",timepoint,"/",list.files(path = paste0("LooNet",timepoint,"/"))[i]), sep = '\t', header = FALSE, row.names = 1)
  head(recons_net[1:5,1:5])
  colnames(recons_net) = rownames(recons_net)
  dim(recons_net)
  try(if(sum(rownames(recons_net) != Sequences_nodes[,2]) > 0 ) stop("NOT MATCHINGs"))
  colnames(recons_net) = rownames(recons_net) = Node_list 
  corr_recons_net = c(as.matrix(recons_net)) #  scales a covariance matrix into the corresponding correlation matrix efficiently.
  
  
  lionessOutput[, i + 2] <- length(files) * (global_net_vect - corr_recons_net) + corr_recons_net
  print(paste0(i, " = " , str(mean(abs( lionessOutput[, i + 2] )))))
}


edges <- paste(lionessOutput[, 1], lionessOutput[, 2], sep = "_")
rownames(lionessOutput) = edges
lionessOutput[1:5,1:5]

# 

dir.create(file.path(result_path), showWarnings = TRUE)
setwd(result_path)

fwrite(lionessOutput[,3:ncol(lionessOutput)],file=paste0("Resulting_net_from_corr_",timepoint,".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames


colnames(Sequences_nodes) = c("Node_list","NAME")
write.table(Sequences_nodes,file = paste0("Sequences_nodes_",timepoint,".txt"), quote = FALSE,row.names = F, col.names = T)
write.table(files_name_matching,file =  paste0("files_name_matching_",timepoint,".txt"))

Resulting_net = lionessOutput[,3:ncol(lionessOutput)]
names_r_net = rownames(Resulting_net)
# colnames(Resulting_net) %in% sample_dataset$sid

length(names_r_net)
# 
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

getwd()
fwrite(Resulting_net,file=paste0("Resulting_net_from_corr_eliminate_",timepoint,".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames



aa = apply(Resulting_net, 1, function(x) sum(abs(x))) > 0
not_null = Resulting_net[aa,]
dim(not_null)
# [1] 547 153
# OR 
# [1] 546 155
fwrite(not_null,file=paste0("Resulting_net_from_corr_eliminate_notNULL_",timepoint,".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames
not_null




# first <- rep(1:122, ncol(global_net))
# second <- rep(1:122, each = nrow(global_net))
# paste(first, second, sep = "_")
# # TRIALS ------------------------------------------------------------------
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

fwrite(not_null_nonode,file=paste0("Resulting_net_from_corr_eliminate_notNULL_nonode_",timepoint,".txt"),
       sep = " ", row.names = TRUE ) # keeps the rownames
not_null_nonode

