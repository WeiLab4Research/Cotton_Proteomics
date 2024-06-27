
########## Latent Class Mixed Model Covariates数据整理所需函数 ##########

## 非时依协变量整理
fixed_cov <- function(cov_vec, nsample = 413, times = 8) {
  
  cov_vec_new <- cov_vec[(1:nsample) * times]
  
  return(cov_vec_new)
}

## 时依协变量整理
if_start <- function(vec, t, times) {
  
  if(t == times) { return(0) }
  if(vec[t] < vec[t + 1]) { return(1) } else { return(0) }
  
}

start_age <- function(vec, times, i, age_matrix) {
  for (t in 1:times) {
    if(if_start(vec = vec, t = t, times = times)) { return(age_matrix[t,i]) }
    if(t == times) { return(age_matrix[times,i]) }
  }
}

if_stop <- function(vec, t, times) {
  
  if(t == times) { return(1) }
  if(vec[t + 1] > vec[t]) { return(0) } else { return(1) }
  
}

stop_age <- function(vec, times, i, age_matrix) {
  for (t in 1:times) {
    if(if_stop(vec = vec, t = t, times = times)) { return(age_matrix[t,i]) }
    if(t == times) { return(age_matrix[times,i]) }
  }
}

time_dep_cov <- function(cov_vec, age_vec, cov_name, nsample = 413, times = 8) {
  cov_matrix <- matrix(cov_vec, nrow = times)
  age_matrix <- matrix(age_vec, nrow = times)
  
  for (i in 1:nsample) {
    for (j in 2:times) {
      if(cov_matrix[j,i] < cov_matrix[j-1,i]) {
        cov_matrix[j,i] <- cov_matrix[j-1,i]
      }
    }
  }
  
  cov_baseline <- cov_matrix[1,]
  
  cov_during <- vector()
  for (i in 1:nsample) {
    if(cov_matrix[times,i] == cov_matrix[1,i]) { cov_during[i] = 0 } else {
      cov_during[i] = stop_age(vec = cov_matrix[,i], times = times, i = i, age_matrix = age_matrix) - 
                     start_age(vec = cov_matrix[,i], times = times, i = i, age_matrix = age_matrix)
    }
  }
  
  cov_cumulation <- cov_matrix[times,]
  
  cov_complete <- data.frame(cov_baseline, cov_during, cov_cumulation)
  colnames(cov_complete) <- c(paste0(cov_name,"_baseline"),
                              paste0(cov_name,"_during"),
                              paste0(cov_name,"_cumulation"))
  return(cov_complete)
}



