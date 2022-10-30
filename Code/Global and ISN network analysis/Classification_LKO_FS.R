## Classification Results
rm(list = ls())
setwd("~/Desktop/R_Root/Microbiome_final/")

library(igraph)
library(ROCR)
library(e1071)
library(randomForest)

LuckiMap1 = read.delim("Data/LuckiMap_6M.txt")
LuckiMap2 = read.delim("Data/LuckiMap_9M.txt")
Children =  merge(LuckiMap1, LuckiMap2, by = "Child")

phenotype = data.frame(delivery_type = Children$delivery_type.x, 
                       DMM = as.character(Children$DMM_clust.x),
                       Diet = as.character(Children$Diet_TP.x))

y = phenotype$delivery_type
y[is.na(y)] = 0
y[y==0] = 1
y = as.factor(y)


X_main = read.delim("Data/OTU_TABLE_69_SUB_CLR_6M.tsv")

X_main = readRDS("Data/dist_MAGMA_2.rds")
# X_main = dist
X_main[is.na(X_main)] = 0
# X_main = t(apply(X_main, 1, function(x){order(x)/length(x)}))

# Feature observation
FeatureCorrelation = apply(X_main, 2, function(x){abs(cor(x,as.numeric(y)))})
FeatureSD = apply(X_main, 2, sd)
FeaturePval = apply(X_main, 2, function(x){
  a = x[y==levels(y)[1]]
  b = x[y==levels(y)[2]]
  t.test(a,b)$p.value
})

hist(FeatureCorrelation, 50)
hist(FeatureSD, 50)
hist(FeaturePval, 50)

MeanAUC = c()
for (nFeat in seq(2,40,1)){
  
  # featureSet = order(FeatureSD, decreasing = TRUE)[1:nFeat]
  # X = X_main[,featureSet] 
  
  X = X_main
  X = scale(X)
  
  ## LOOCV
  N = nrow(X)
  Rep = 2
  kernelParam = 2
  
  Result = list()
  Result[["AUC"]] = c()
  for (rep in 1:Rep){
    yy = c()
    y_hat = c()
    k = 1
    for (i in 1:N){
      sampled_one =sample(x = seq(1:N)[-i],size = k-1,replace = F)
      
      # X_test = X[c(sampled_one,i),]
      y_test = y[c(sampled_one,i)]
      X_train = X[-c(sampled_one,i),]
      y_train = y[-c(sampled_one,i)]
      
      # FeatureCorrelation_internal = apply(X_train, 2, function(x){abs(cor(x,as.numeric(y_train)))})
      FeaturePval_internal = apply(X_main, 2, function(x){
        a = x[y==levels(y)[1]]
        b = x[y==levels(y)[2]]
        t.test(a,b)$p.value
      })
      featureSet = order(FeaturePval_internal, decreasing = FALSE)[1:nFeat]
      X_train = X_train[,featureSet]
      X_test = X[c(sampled_one,i), featureSet]
      
      # Down sampling the majority class
      n1 = table(y_train)[1]
      n2 = table(y_train)[2]
      majorityClass = ifelse(n1>n2,1,2)
      n_resample = abs(n1-n2)
      
      indeces_Cm = which(y_train == levels(y_train)[majorityClass])
      indeces_Cm = sample(indeces_Cm, n_resample, replace = FALSE)
      
      X_train_resampled = X_train[-indeces_Cm,]
      y_train_resampled = y_train[-indeces_Cm]
      
      model = svm(y = as.factor(y_train_resampled), x = data.frame(X_train_resampled), 
                  kernel = "radial", gamma = kernelParam, scale = F)
      
      y_hat = c(y_hat, predict(model, newdata = data.frame(rbind(X_test))))
      
      #model = randomForest(y = as.factor(y_train_resampled),x = X_train_resampled, ntree = 50 ,mtry = 10)
      #y_hat = c(y_hat, predict(model, newdata = rbind(X_test)))
      
      yy = c(yy, y_test)
    }
    
    pred = prediction(y_hat, yy)
    auc = performance(pred, measure = "auc")
    
    Result[["AUC"]] = c(Result[["AUC"]], as.numeric(auc@y.values))
  }
  
  MeanAUC = c(MeanAUC, mean(Result[["AUC"]]))
  sd(Result[["AUC"]])
  
}


plot(MeanAUC)
max(MeanAUC)

#table(y,y_hat)
#W = t(model$coefs) %*% model$SV
