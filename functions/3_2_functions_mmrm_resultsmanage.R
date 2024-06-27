
########## Mixed Model Repeated Measurements 结果整理所需函数 ##########

# 根据60y与cotton:age整理其他年龄时蛋白质的效应
protein_effect <- function(target_age, origin_age = 60, results_origin, results_inter, inter = F) {
  
  if(target_age == origin_age) {
    results_target <- results_origin
    colnames(results_target)[1:4] <- c("Protein", "beta", "se", "Pval")
    
    if(inter) {
      results_target <- results_inter
      colnames(results_target)[1:6] <- c("Protein", "beta", "se", "Pval", "cov", "num")
    }
    
    results_target <- results_target[which(results_target$Pval != "NaN"),]
    results_target$beta <- as.numeric(results_target$beta)
    results_target$Pval <- as.numeric(results_target$Pval)
    results_target$FDR <- p.adjust(results_target$Pval, method = "BH")
    results_target$Sig_Bon <- results_target$Pval <= 0.05 / nrow(results_target)
    ## Pval判断显著性（P<0.05 & beta>0 = "Up"; P<0.05 & beta<0 = "Down"; P>=0.05 = "Nosig")
    results_target$change = ifelse(results_target$Pval < 0.05, 
                                   ifelse(results_target$beta > 0 , 'Up', 'Down'),
                                   'Nosig')
    
    return(results_target)
  }
  
  ## 计算方法：
  ## β_y = β_60 + β_age * (y-60)
  ## se_y^2 = se_60^2 + 2 * (y-60) * cov(β_60, β_age) + se_age^2
  results_target_age <- data.frame(Protein = results_origin$Protein, 
                                   beta = results_origin$beta + results_inter$beta * (target_age - 60), 
                                   num = as.numeric(results_inter$num))
  
  results_target_age$se <- cbind(results_origin[,"se"],results_inter[,c("se","cov")]) %>%
    apply(1, SE_age, age = target_age)
  
  results_target_age$Pval <- results_target_age[,c("beta","se","num")] %>% 
    apply(1,P.value)
  
  results_target_age$FDR <- p.adjust(results_target_age$Pval, method = "BH")
  results_target_age$Sig_Bon <- results_target_age$Pval <= 0.05/nrow(results_target_age)
  ## Pval判断显著性（P<0.05 & beta>0 = "Up"; P<0.05 & beta<0 = "Down"; P>=0.05 = "Nosig")
  results_target_age$change = ifelse(results_target_age$Pval < 0.05, 
                                     ifelse(results_target_age$beta > 0 ,'Up','Down'),
                                     'Nosig')
  
  return(results_target_age)
}

## MMRM结果整理
mmrm_results_format <- function(age_vec = c(60,70,80,90),
                                age_inter = TRUE,
                                end_inter = FALSE) {
  
  # 读取蛋白质名称信息
  proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
  colnames(proteinlist)[1] <- "Protein_ID"
  
  # 读取蛋白质注释信息
  annotation_protein <- read.table("original data/report/data/Identification/annotation_allprotein.xls",
                                   sep = "\t",
                                   header = TRUE,
                                   fill = TRUE,
                                   quote = "")
  annotation_protein <- annotation_protein[,-2]
  
  results2write <- results_60y_ordered %>%
    left_join(results_70y_ordered, by = "Protein") %>%
    left_join(results_80y_ordered, by = "Protein") %>%
    left_join(results_90y_ordered, by = "Protein")
  if(age_inter) {
    results2write <- results2write %>%
      left_join(results_age_ordered, by = "Protein")
  }
  if(end_inter) {
    results2write <- results2write %>%
      left_join(results_end_ordered, by = "Protein")
  }
  
  results2write <- results2write %>%
    left_join(proteinlist, by = "Protein") %>%
    left_join(annotation_protein, by = "Protein_ID")
  sig_row <- vector()
  for (i in c(age_vec)) {
    age_results <- get(paste0("results_",i,"_ordered"))
    sig_row_age <- age_results[which(age_results$FDR < 0.05),"Protein"]
    sig_row <- unique(c(sig_row,sig_row_age))
  }
  
  results2write <- results2write[which(results2write$Protein %in% sig_row),]
  
  results2write <- results2write[,-which(str_detect(colnames(results2write), pattern = "cov|num|change|Sig_Bon|gene.name"))]
  colnames(results2write)[2:21] <- c("beta_60y", "se_60y", "Pval_60y", "FDR_60y",
                                     "beta_70y", "se_70y", "Pval_70y", "FDR_70y",
                                     "beta_80y", "se_80y", "Pval_80y", "FDR_80y",
                                     "beta_90y", "se_90y", "Pval_90y", "FDR_90y",
                                     "beta_age", "se_age", "Pval_age", "FDR_age")
  
  return(results2write)
}

