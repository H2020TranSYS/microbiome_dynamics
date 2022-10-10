
# EXTERN Feature ranking --------------------------------------------------

X_main_only_eq = X_main[y!=0,]
vv = round(X_main_only_eq,5)!= 0
nozero = apply(vv, 2, sum)
print("EDGE")


# CORRECTION to have no NA ------------------------------------------------
if (type_data == "Edge"){
  if (rankfeature == "internal"){
    X_main_only_eq = X_main_only_eq[,nozero > 1]
  } else {X_main_only_eq = X_main_only_eq[,nozero > 0]}
}

y = factor(y[y!=0])
# Feature selection
FeatureCorrelation = apply(X_main_only_eq, 2, function(x){abs(cor(x,as.numeric(y)))})
hist(FeatureCorrelation, 20)
featureSet = which(FeatureCorrelation > .2)
