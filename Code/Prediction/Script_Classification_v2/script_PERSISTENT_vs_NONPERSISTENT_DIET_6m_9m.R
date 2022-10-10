

source("Script_Classification_v2/Importing_preprocessing_set_seed_step1.R")


phenotype = data.frame(delivery_type = Children$delivery_type.x, 
                       DMM = as.character(Children$DMM_clust.x),
                       Diet6 = as.character(Children$Diet_TP.x),
                       Diet9 = as.character(Children$Diet_TP.y)
)
# print("AAAA")
table(phenotype$Diet6,phenotype$Diet9)
# print("BBBB")
y = rep(0,nrow(phenotype))
print(phenotype$Diet6 )
print(phenotype$Diet9 )
phenotype$Diet9 <- factor(phenotype$Diet9, levels=levels(phenotype$Diet6))

y[phenotype$Diet6 != phenotype$Diet9 ] = 1
y[phenotype$Diet6 == phenotype$Diet9 ] = 2

print("CCCC")
source("Script_Classification_v2/Type_of_data_step2.R")
print("DDDD")

source("Script_Classification_v2/further_correction_step3.R")
print("EEEE")

source("Script_Classification_v2/Real_SVM_step4.R")
print("FFFF")



finale_name = paste0(gsub("\\.[a-z]+", "",reading_file), ifelse(orders, "_ORDERED","_NO_ORDERED"))
setwd("DIET_PERSISTENT_vs_NONPERSISTENT")

dir.create(file.path(finale_name), showWarnings = FALSE)

setwd(file.path(finale_name))

plot(MeanAUC)
max(MeanAUC)

FeatureCorrelation = apply(X_main_only_eq, 2, function(x){abs(cor(x,as.numeric(y)))})
nFeat = which.max(MeanAUC)+1
featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
topFeatures = data.frame(feature_top = featureSet, Taxas = colnames(X_main_only_eq[,featureSet] ))

write.table(topFeatures,  paste0("TopFeatures_",finale_name,"_", rankfeature,".txt") )


rankfeature ## two value external or internal
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

print("ENDD")


#table(y,y_hat)

#W = t(model$coefs) %*% model$SV

