# if (!exists("args")) {
#   suppressPackageStartupMessages(library("argparse"))
#   parser <- ArgumentParser()
#   parser$add_argument("-a", "--arg1", type="character", defalt="a",
#                       help="First parameter [default %(defult)s]")
#   parser$add_argument("-b", "--arg2", type="character", defalt="b",
#                       help="Second parameter [default %(defult)s]")
#   args <- parser$parse_args()
# }

print (args)

print(args$setted_wd)
print(args$reading_file)
print(args$rankfeature)





# Parameters --------------------------------------------------------------
Rep = args$rep # TODO mmettre a l'haut and as setting

setted_wd = args$setted_wd
reading_file = args$reading_file
rankfeature = args$rankfeature
type_data = args$type_data
taxas_n = args$colnm
orders = args$ordering
## Classification Results
# rm(list = ls())
setwd(setted_wd)
# setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/CLASSIFICATION/")

library(igraph)
library(ROCR)
library(e1071)

set.seed(123)


LuckiMap1 = read.delim("Data/LuckiMap_6M.txt")
LuckiMap2 = read.delim("Data/LuckiMap_9M.txt")
Children =  merge(LuckiMap1, LuckiMap2, by = "Child")
