
# Libraries ---------------------------------------------------------------


library(xlsx)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
otu_table_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data/MAGMA_data/6M"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Results/Selected_95_microbes"

setwd(data_path)

Taxonomy_ASVs = read.csv("Taxonomy_ASVs.csv")



# Example with the data of the library ------------------------------------

setwd(otu_table_path)


# 81  -------------------------------------------------------


## Whatever OTU table would be the same
otu_table = as.matrix(read.table("OTU_TABLE_MAGMA6M.tsv"))

Taxonomy_ASVs[2,]
ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(otu_table),]
df = (t(ASV_final))
dim(df)
colnames(df) = df[1,]
df = df[-1,]


rownames(df)


# Solitary Phylums --------------------------------------------------------

setwd(result_path)


df[1:5,1:5]
(TM = df[,df["Phylum",] == "p__TM7"] )
# Kingdom        Phylum         Class         Order        Family         Genus       Species 
# "k__Bacteria"      "p__TM7"    "c__TM7-3"         "o__"         "f__"         "g__"         "s__" 
write(TM, file = "TM7_solitary_phylum_specifications.txt")
(VERRUCO = df[,df["Phylum",] == "p__Verrucomicrobia"] )
write(VERRUCO, file = "Verrucomicrobia_solitary_phylum_specifications.txt")

# Kingdom                   Phylum                    Class                    Order                   Family                    Genus                  Species 
# "k__Bacteria"     "p__Verrucomicrobia"    "c__Verrucomicrobiae"  "o__Verrucomicrobiales" "f__Verrucomicrobiaceae"         "g__Akkermansia"         "s__muciniphila" 




## All the selected 95 microbes
df
class(df)
write.table(df, file = "Complete_annotated_list_95microbes.txt", quote = F, sep = "\t", row.names = T, col.names = T)

write.xlsx2(df, file = "Complete_annotated_list_95microbes.xlsx", sheetName = "Microbes_list" )
