# THE SAME BUT WITH CORRELATION -------------------------------------------
library(data.table)



data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/SparCC_data/"


setwd(data_path)
getwd()
OTU_TABLE_MAGMA6M <- read.delim("6M/OTU_TABLE_MAGMA6M.tsv")
OTU_TABLE_MAGMA9M <- read.delim("9M/OTU_TABLE_MAGMA9M.tsv")

Sparcc_OTU_TABLE6M = t(OTU_TABLE_MAGMA6M)
Sparcc_OTU_TABLE9M = t(OTU_TABLE_MAGMA9M)

setwd(result_path)
write.table(x = Sparcc_OTU_TABLE9M, file = "OTU_TABLE_input_sparCC_9M.tsv", sep = "\t",row.names = T, quote = F)
write.table(x = Sparcc_OTU_TABLE6M, file = "OTU_TABLE_input_sparCC_6M.tsv", sep = "\t",row.names = T, quote = F)
