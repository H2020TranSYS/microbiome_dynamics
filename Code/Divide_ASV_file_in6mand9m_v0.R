library(readxl)



data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"

setwd(data_path)

SparCC_ASV <- read.delim("SparCC_ASV.tab")

LuckiMap_anonym_long_format <- read_excel("LuckiMap_anonym_long_format.xlsx")


# Divide per month --------------------------------------------------------


LuckiMap_anonym_long_format_6m = LuckiMap_anonym_long_format[LuckiMap_anonym_long_format$Age_months == 6,]
LuckiMap_anonym_long_format_9m = LuckiMap_anonym_long_format[LuckiMap_anonym_long_format$Age_months == 9,]

SparCC_ASV[1:5,1:5]


SparCC_ASV6M = cbind(SparCC_ASV$X.OTU.ID,SparCC_ASV[,colnames(SparCC_ASV) %in% LuckiMap_anonym_long_format_6m$Sample.name ])
SparCC_ASV9M = cbind(SparCC_ASV$X.OTU.ID,SparCC_ASV[,colnames(SparCC_ASV) %in% LuckiMap_anonym_long_format_9m$Sample.name ])

colnames(SparCC_ASV6M)
rownames(SparCC_ASV6M)

setwd(result_path)


# Write the 6 and 9 months ------------------------------------------------


write.table(SparCC_ASV6M, file = "SparCC_ASV6M.tab",quote = F,sep = "\t", row.names = F, col.names = T)
write.table(SparCC_ASV9M, file = "SparCC_ASV9M.tab",quote = F,sep = "\t", row.names = F, col.names = T)
write.table(LuckiMap_anonym_long_format_6m,file = "LuckiMap_anonym_long_format_6m.txt",quote = F,sep = "\t", row.names = F, col.names = T)
write.table(LuckiMap_anonym_long_format_9m,file = "LuckiMap_anonym_long_format_9m.txt",quote = F,sep = "\t", row.names = F, col.names = T)

LuckiMap_anonym_long_format_6m[1:5,1:5]
## n.b. after save we have to add #OTU ID manually 
SparCC_ASV6M[1:5,1:5]


# ## n.b. after save we have to add #OTU ID manually y  -------------------------------------------------------------------------


