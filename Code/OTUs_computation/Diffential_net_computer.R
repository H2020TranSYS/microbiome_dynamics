
# LIBRARIES and wd
library(data.table)
data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Results/OTUS_node_analysis_difference/"
setwd(data_path)



## DATA IMPORTING

Resulting_net_6M <- read.csv("./6M/ISN/Resulting_net_from_corr_MAGMA.txt", row.names=1, sep="")
Resulting_net_9M <- read.csv("./9M/ISN/Resulting_net_from_corr_MAGMA.txt", row.names=1, sep="")

# node check
Sequences_nodes6M <- read.csv("./6M/Sequences_nodes.txt", sep="")
Sequences_nodes9M <- read.csv("./9M/Sequences_nodes.txt", sep="")

Sequences_nodes6M == Sequences_nodes9M

## YES --> ALL of them matching!

# comparison
LuckiMap_9M <- read.delim("../LuckiMap_anonym_long_format_9m.txt")
# C   View(LuckiMap_9M)
LuckiMap_6M <- read.delim("../LuckiMap_anonym_long_format_6m.txt")
# >   View(LuckiMap_6M)

Children =  merge(LuckiMap_6M, LuckiMap_9M, by = "Child")



# MATCHING ----------------------------------------------------------------

name_vect = rownames(Resulting_net_6M)



cc       <- strsplit(name_vect,'_')
part1    <- unlist(cc)[2*(1:length(name_vect))-1]
part2    <- unlist(cc)[2*(1:length(name_vect))  ]
head(part1) ; head(part2) ; head(name_vect)

# select paired samples  ----------------------------------------------------------
Resulting_net_6M_f = Resulting_net_6M[,names(Resulting_net_6M) %in% Children$Sample.name.x ]
Resulting_net_9M_f = Resulting_net_9M[,names(Resulting_net_9M) %in% Children$Sample.name.y ]



names(Resulting_net_6M_f) <- Children$Child[match(names(Resulting_net_6M_f), Children$Sample.name.x)]
names(Resulting_net_9M_f) <- Children$Child[match(names(Resulting_net_9M_f), Children$Sample.name.y)]

Ord_res_net_9M = Resulting_net_9M_f[,names(Resulting_net_6M_f)]
## put the IND in the same order 6M and 9M
names(Ord_res_net_9M) == names(Resulting_net_6M_f)

## Same order


head(Ord_res_net_9M[,1:5]); head(Resulting_net_6M_f[,1:5])
names(Ord_res_net_9M); names(Resulting_net_6M_f)



# Differential network ----------------------------------------------------

setwd(result_path)

differential_net = Ord_res_net_9M - Resulting_net_6M_f

differential_net
fwrite(differential_net,file="Resulting_DIFF_net_from_corr_eliminate_MAGMACONF.txt",
       sep = " ", row.names = TRUE ) # keeps the rownames


# var(c(0,0,0,0))

aa = apply(differential_net, 1, function(x) sum(abs(x))) > 0
not_null = differential_net[aa,]
dim(not_null)
# [1] 430  81

getwd()
fwrite(not_null,file="Resulting_DIFF_net_notNULL_MAGMACONF.txt",
       sep = " ", row.names = TRUE ) # keeps the rownames


fwrite(Children, file= "Children_matching_table6M-9M.txt", sep = " ")




# Net 6m ------------------------------------------------------------------

fwrite(Resulting_net_6M_f,file="Resulting_6M_net_from_corr_eliminate_MAGMACONF.txt",
       sep = " ", row.names = TRUE ) # keeps the rownames


# var(c(0,0,0,0))

aa = apply(Resulting_net_6M_f, 1, function(x) sum(abs(x))) > 0
not_null = Resulting_net_6M_f[aa,]
dim(not_null)
# [1] 430  81

getwd()
fwrite(Resulting_net_6M_f,file="Resulting_6M_net_notNULL_MAGMACONF.txt",
       sep = " ", row.names = TRUE ) # keeps the rownames


# Net 9m ------------------------------------------------------------------

fwrite(Resulting_net_9M_f,file="Resulting_9M_net_from_corr_eliminate_MAGMACONF.txt",
       sep = " ", row.names = TRUE ) # keeps the rownames


# var(c(0,0,0,0))

aa = apply(Resulting_net_9M_f, 1, function(x) sum(abs(x))) > 0
not_null = Resulting_net_9M_f[aa,]
dim(not_null)
# [1] 430  81

getwd()
fwrite(Resulting_net_9M_f,file="Resulting_9M_net_notNULL_MAGMACONF.txt",
       sep = " ", row.names = TRUE ) # keeps the rownames


