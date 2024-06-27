

get_imputed <- function(dat, method = "min") {
  
  dat_imputed <- dat
  
  vec_protein <- which(str_detect(colnames(dat), pattern = "protein"))
  
  if(method == "min") {
    for (i in vec_protein) {
      min_val = min(dat_imputed[,i], na.rm = T)
      dat_imputed[which(is.na(dat_imputed[,i])),i] <- min_val
    }
  }
  
  if(method == "mean") {
    for (i in vec_protein) {
      mean_val = mean(dat_imputed[,i])
      dat_imputed[which(is.na(dat_imputed[,i])),i] <- mean_val
    }
  }
  
  if(method == "knn") {
    require(DMwR2)
    t_dat_imputed <- t(dat_imputed[,vec_protein])
    t_dat_imputed <- knnImputation(data = t_dat_imputed, k = 10, scale = F, meth = "weighAvg")
    dat_imputed[,vec_protein] <- t(t_dat_imputed)
  }
  
  if(method == "seqknn") {
    require(multiUS)
    t_dat_imputed <- t(dat_imputed[,vec_protein])
    t_dat_imputed <- seqKNNimp(data = t_dat_imputed, k = 10)
    dat_imputed[,vec_protein] <- t(t_dat_imputed)
  }
  
  return(dat_imputed)
}
