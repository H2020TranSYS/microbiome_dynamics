
MeanAUC = c()
sdAUC = c()
for (nFeat in seq(2,30,1)){

  featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]

  X = X_main_only_eq[,featureSet]
  X = scale(X)

  ## LOOCV
  N = nrow(X_main_only_eq)
  kernelParam = .1
  Result = list()
  Result[["AUC"]] = c()

  for (rep in 1:Rep){
    y_hat_cont = rep(NA, N)
    y_hat = rep(NA, N)

    for (i in 1:N){
      if (rankfeature == "internal"){

        FeatureCorrelation = apply(X_main_only_eq[-i,], 2, function(x){abs(cor(x,as.numeric(y[-i])))})
        featureSet = which(FeatureCorrelation > .2)
        featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
        X = X_main_only_eq[,featureSet]
        X = scale(X)        ## DO we really want to scale for CLR?
      }


      X_test = X[i,]
      y_test = y[i]
      X_train = X[-i,]
      y_train = y[-i]

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

      y_hat_cont[i] = predict(model, newdata = data.frame(rbind(X_test)))
      y_hat[i] = ifelse(y_hat_cont[i]>=1.5,2,1)

    }

    pred = prediction(y_hat, y)
    auc = performance(pred, measure = "auc")

    Result[["AUC"]] = c(Result[["AUC"]], as.numeric(auc@y.values))
  }

  MeanAUC = c(MeanAUC, mean(Result[["AUC"]]))
  sdAUC = c(sdAUC, sd(Result[["AUC"]]))

}
