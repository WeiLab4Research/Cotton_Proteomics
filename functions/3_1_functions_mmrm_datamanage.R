
########## Mixed Model Repeated Measurements 数据整理所需函数 ##########

# 计算NA比例 即蛋白质缺失率
na_ratio <- function(x) {
  return(mean(is.na(x)))
}

# 计算某age时 β_protein + β_age:protein*age 的标准误
SE_age <- function(result, age) {
  se <- sqrt(as.numeric(result[1])^2 + (as.numeric(result[2])*(age-60))^2 + 2*(age-60)*as.numeric(result[3]))
  return(se)
}

# 根据均值及标准误计算P值
P.value <- function(cal_P) {
  tstat <- cal_P[1] / cal_P[2]
  P <- 2 * pt(abs(tstat), df = cal_P[3]-12, lower.tail = FALSE)
  return(P)
}

# MMRM数据整理
get_MMRMData <- function(varvec_repeated, varvec_single, age_central = 60) {
  
  ## 重复测量变量
  repeated_data <- data.frame(pid = rep(cotton_protein$pid, each = 8))
  for (i in varvec_repeated) {
    repeated_data_each <- cotton_protein %>%
      dplyr::select(starts_with(i)) %>%
      dplyr::mutate(pid = cotton_protein$pid) %>%
      pivot_longer(cols = !pid,
                   names_to = "names",
                   values_to = i) %>%
      dplyr::select(starts_with(i))
    
    repeated_data <- repeated_data %>% cbind(repeated_data_each)
  }
  
  ## 单次测量变量
  single_data <- data.frame(pid = rep(cotton_protein$pid, each = 8))
  for (j in varvec_single) {
    single_data_each <- rep(cotton_protein[,j], each = 8)
    
    single_data <- single_data %>% cbind(single_data_each)
  }
  single_data <- single_data[,-1]
  colnames(single_data) <- varvec_single
  
  ## proteins变量
  valid_proteins_imputed <- valid_proteins_imputed %>% apply(2, rep, each = 8)
  
  ## 变量合并
  mmrmdata <- cbind(repeated_data, single_data, valid_proteins_imputed)
  
  mmrmdata$age <- mmrmdata$age - age_central
  
  return(mmrmdata)
}

select_proteins <- function(dat, limit = 0.5) {
  proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
  protein_vec <- proteinlist[which(proteinlist$gene.name != ""),"Protein"]
  valid_proteins <- dat %>% 
    dplyr::select(starts_with("protein")) %>%
    dplyr::select(one_of(protein_vec)) %>%
    select_if(~ na_ratio(.) < limit)
  return(valid_proteins)
}

