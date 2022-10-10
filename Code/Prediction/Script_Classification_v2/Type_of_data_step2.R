### If is using IS-edges; or if it is using node-value, enfin if it is using EDNN
if (type_data == "Edge" ){
  data = data.frame(read.delim(file=paste0("Data/",reading_file),sep = " " ))
  data$X = as.character(data$X)
  # print(data$X)
  Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
  # MianGraph = graph(Nodelist, directed = FALSE)
  IndvGraphWeights = data[,-1]
  # IndvGraphWeights = log(abs(IndvGraphWeights)+1)
  d = t(IndvGraphWeights)
  X_main = d
  X_main = data.frame(X_main)
  colnames(X_main) = data$X
  vv = round(X_main,5)!= 0
  nozero = apply(vv, 2, sum)
  print("EDGE")
  if (rankfeature == "internal"){ X_main = X_main[,nozero > 1] }
  
  
} else if(type_data == "EDNN"){
  X_main = readRDS(paste0("Data/",reading_file))
  X_main[is.na(X_main)] = 0
  # X_main = X_main + matrix(runif(length(X_main) ,min = 0, max = .001), nrow(X_main),ncol(X_main))
  if (orders){
    X_main = t(apply(X_main, 1, function(x){order(x)/length(x)}))
  }
  X_main = data.frame(X_main)
  colnames(X_main) = taxas_n
  
  print("EDNN")
} else {
  X_main = read.delim(paste0("Data/",reading_file))
  print("NODE")
  
  
}