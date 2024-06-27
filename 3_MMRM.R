library(tidyverse)
library(nlme)
library(ggrepel)
library(patchwork)
# library(geepack)
# library(emmeans)
setwd("E:/cotton")
rm(list = ls())
source("codes/functions/1_functions_loadall.R")

########## Inputs ##########

cotton_protein <- read.csv("original data/cotton35_proteins.csv", header = T, row.names = 1)

# load("MMRM_new.RData")

########## Data Reform ##########

valid_proteins <- select_proteins(cotton_protein, limit = 0.5)

valid_proteins_imputed <- get_imputed(valid_proteins, method = "seqknn")

mmrmdata <- get_MMRMData(varvec_repeated = c("fev","fvc","age","pkyrs","end"),
                         varvec_single = c("male81", "fev81", "cesyr16", "height16"),
                         age_central = 60)

########## Mixed Model Repeated Measurements Analysis ##########

# protein主效应(60岁)
results_60y <- data.frame()
# protein:age交互作用
results_age <- data.frame()
# protein:cotton交互作用
results_cotton <- data.frame()

vec_protein <- which(str_detect(colnames(mmrmdata), pattern = "protein"))
for (i in vec_protein) {
  ## 
  temp_mmrmdata <- mmrmdata[!is.na(mmrmdata[,i]),]
  colnames(temp_mmrmdata)[i] <- "protein"
  
  ## GEE模型有误 后续可尝试调整
  # geemodel <- geeglm(FEV1 ~ temp_geedata[,i]*cotton + age0 + sex + cigevr + working + FEV10,
  # family = gaussian,
  # data = temp_geedata,
  # id = interaction(pid, age), 
  # corstr = "ar1")
  
  mmrm_model <- gls(fev ~ protein*age + protein*end + male81 + cesyr16 + pkyrs + height16 + fev81,
                    na.action=na.omit,
                    data=temp_mmrmdata,
                    correlation=nlme::corAR1(form = ~1 | pid),
                    method = "REML")
  # 计算协方差矩阵
  cov_mmrm_model <- vcov(mmrm_model, complete = T)
  
  # 统计结果
  ## 格式：
  ## ——————————————————————————————————————————
  ## Protein       b        se      P      cov(proteins&interactions)
  ## ——————————————————————————————————————————
  ## Protein1      0.55    0.34  1.12E-1  1.23
  ## Protein...     ...    ...    ...     ...
  ## ——————————————————————————————————————————
  results_60y <- rbind(results_60y,
                       c(colnames(mmrmdata)[i],coef(summary(mmrm_model))[2,c(1,2,4)]))
  results_age <- rbind(results_age,
                       c(colnames(mmrmdata)[i],
                         coef(summary(mmrm_model))[10,c(1,2,4)],
                         cov_mmrm_model["protein","protein:age"],
                         nrow(temp_mmrmdata)))
  results_cotton <- rbind(results_cotton,
                          c(colnames(mmrmdata)[i],
                            coef(summary(mmrm_model))[11,c(1,2,4)],
                            cov_mmrm_model["protein","protein:end"]))
  print(colnames(mmrmdata)[i])
}

# 整理60岁protein数据
results_60y <- protein_effect(target_age = 60,
                              origin_age = 60,
                              results_origin = results_60y,
                              results_inter = results_age)

# 整理proten:age交互作用数据
results_age <- protein_effect(target_age = 60,
                              origin_age = 60,
                              results_origin = results_60y,
                              results_inter = results_age,
                              inter = T)

# 整理其他年龄 (70/80/90岁) protein数据
for (i in c(70,80,90)) {
  
  varname <- paste0("results_", i, "y")
  varname_ordered <- paste0("results_", i, "y_ordered")
  assign(varname, 
         protein_effect(target_age = i, 
                        origin_age = 60,
                        results_origin = results_60y,
                        results_inter = results_age))
  ## 排序
  assign(varname_ordered,get(varname)[order(get(varname)$Pval),])
  
}

## 60岁protein数据排序
results_60y_ordered <- results_60y[order(results_60y$Pval),]
## protein:age交互作用排序
results_age_ordered <- results_age[order(results_age$Pval),]

#### 添加gene.name
proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
colnames(proteinlist)[1] <- "Protein_ID"

results_60y_ordered$gene.name <- ""
results_60y_ordered[1:5,"gene.name"] <- apply(results_60y_ordered[1:5,],1,match_gene_name)
# results_60y_ordered[which(results_60y_ordered$Protein == "protein1221"), "gene.name"] <- "HBB"

results_70y_ordered$gene.name <- ""
results_70y_ordered[1:5,"gene.name"] <- apply(results_70y_ordered[1:5,],1,match_gene_name)
# results_70y_ordered[which(results_70y_ordered$Protein == "protein1221"), "gene.name"] <- "HBB"

results_80y_ordered$gene.name <- ""
results_80y_ordered[1:5,"gene.name"] <- apply(results_80y_ordered[1:5,],1,match_gene_name)
# results_80y_ordered[which(results_80y_ordered$Protein == "protein1221"), "gene.name"] <- "HBB"

results_90y_ordered$gene.name <- ""
results_90y_ordered[1:5,"gene.name"] <- apply(results_90y_ordered[1:5,],1,match_gene_name)
# results_90y_ordered[which(results_90y_ordered$Protein == "protein1221"), "gene.name"] <- "HBB"

results_age_ordered$gene.name <- ""
results_age_ordered[1:5,"gene.name"] <- apply(results_age_ordered[1:5,],1,match_gene_name)
# results_age_ordered[which(results_age_ordered$Protein == "protein1221"), "gene.name"] <- "HBB"

########## Data Integration ##########

results_all <- cbind(results_60y, results_70y, results_80y, results_90y, results_age)

results2write <- mmrm_results_format(age_vec = c("60y","70y","80y","90y","age"),
                                     age_inter = T)

########## Volcano Plots ##########

mmrm_60y_plot <- volcano_plot(results_60y_ordered, 
                              titles = "A protein effects at 60y",
                              FDR = F,
                              gene_name = T) +
  geom_hline(yintercept = -log10(min(results_60y_ordered$Pval)) + 0.1, linetype = 2)

mmrm_age_plot <- volcano_plot(results_age_ordered, 
                              inter = T, 
                              titles = "Protein*age interactions",
                              FDR = T,
                              gene_name = F)

for (i in c(70,80,90)) {
  title <- data.frame(y=c(70,80,90), t=c("B","C","D"))
  varname_plot <- paste0("mmrm_", i, "y_plot")
  varname_ordered <- paste0("results_", i, "y_ordered")
  assign(varname_plot,
         volcano_plot(get(varname_ordered), 
                      titles = paste0(title[which(title$y==i),"t"]," protein effects at ",i,"y"),
                      FDR = T,
                       gene_name = T))
}

patch_plot <- mmrm_60y_plot + mmrm_70y_plot + mmrm_80y_plot + mmrm_90y_plot


########## Outputs ##########

save(results_all, file = "results_knn/RData/3_MMRM/results_mmrm.RData")

# 输出60岁protein火山图
jpeg("results_knn/figures/3_MMRM/MMRM_60y_plot.jpg", width = 15, height = 15, units = "cm",res = 500,family="serif")
mmrm_60y_plot
dev.off()

# 输出protein*age交互作用火山图
jpeg("results_knn/figures/3_MMRM/MMRM_age_plot.jpg", width = 15, height = 15, units = "cm",res = 500,family="serif")
mmrm_age_plot
dev.off()

# 输出70岁protein火山图
jpeg("results_knn/figures/3_MMRM/MMRM_70y_plot.jpg", width = 15, height = 15, units = "cm",res = 500,family="serif")
mmrm_70y_plot
dev.off()

# 输出80岁protein火山图
jpeg("results_knn/figures/3_MMRM/MMRM_80y_plot.jpg", width = 15, height = 15, units = "cm",res = 500,family="serif")
mmrm_80y_plot
dev.off()

# 输出90岁protein火山图
jpeg("results_knn/figures/3_MMRM/MMRM_90y_plot.jpg", width = 15, height = 15, units = "cm",res = 500,family="serif")
mmrm_90y_plot
dev.off()

# 输出合并的protein火山图
jpeg("results_knn/figures/3_MMRM/MMRM_patch_plot.jpg", width = 25, height = 25, units = "cm",res = 500,family="serif")
patch_plot
dev.off()


# 保存工作空间
save(list = ls(), file = "results_knn/RData/3_MMRM/RData_mmrm.RData")
