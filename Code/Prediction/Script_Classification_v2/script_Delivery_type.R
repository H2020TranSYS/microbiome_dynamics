# if (!exists("args")) {
#   suppressPackageStartupMessages(library("argparse"))
#   parser <- ArgumentParser()
#   parser$add_argument("-a", "--arg1", type="character", defalt="a",
#                       help="First parameter [default %(defult)s]")
#   parser$add_argument("-b", "--arg2", type="character", defalt="b",
#                       help="Second parameter [default %(defult)s]")
#   args <- parser$parse_args()
# }

print("AAA")
source("Script_Classification_v2/Importing_preprocessing_set_seed_step1.R")


phenotype = data.frame(delivery_type = Children$delivery_type.x, 
                       DMM = as.character(Children$DMM_clust.x),
                       Diet = as.character(Children$Diet_TP.x))

y = phenotype$delivery_type
y[is.na(y)] = 0
# y[y==0] = 1
y = as.factor(ifelse(y == "Vaginal",1,2))

# print("BBB")

source("Script_Classification_v2/Type_of_data_step2.R")
source("Script_Classification_v2/further_correction_step3.R")
source("Script_Classification_v2/Real_SVM_step4.R")

finale_name = paste0(gsub("\\.[a-z]+", "",reading_file), ifelse(orders, "_ORDERED","_NO_ORDERED"))
setwd("SVM")

dir.create(file.path(finale_name), showWarnings = FALSE)

setwd(file.path(finale_name))

plot(MeanAUC)
max(MeanAUC)

FeatureCorrelation = apply(X_main, 2, function(x){abs(cor(x,as.numeric(y)))})
nFeat = which.max(MeanAUC)+1
featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
topFeatures = data.frame(feature_top = featureSet, Taxas = colnames(X_main[,featureSet] ))

write.table(topFeatures,  paste0("TopFeatures_",finale_name,"_", rankfeature,".txt") )


## two value external or internal
write(MeanAUC, paste0("MeanAUC_",finale_name,"_", rankfeature,".txt"),)
write(max(MeanAUC), paste0("MaxAUC_",finale_name,"_", rankfeature,".txt"))
write(sdAUC, paste0("sdAUC_",finale_name,"_", rankfeature,".txt"),)
write(sdAUC[which.max(MeanAUC)], paste0("MAXsdAUC_",finale_name,"_", rankfeature,".txt"))


df = data.frame(max(MeanAUC),sdAUC[which.max(MeanAUC)],paste(finale_name,rankfeature,sep = "_"))
write.table(df,paste0("../table_aucc_",rankfeature,".txt"),append = T, col.names = F, row.names = F, quote = F)

png(filename =  paste0("PlotAUC_",finale_name,"_", rankfeature,".png"),)
plot(MeanAUC)
dev.off()


png(filename =  paste0("histAUC_",reading_file,"_", rankfeature,".png"),)
hist(FeatureCorrelation, 20)
dev.off()




#table(y,y_hat)

#W = t(model$coefs) %*% model$SV

