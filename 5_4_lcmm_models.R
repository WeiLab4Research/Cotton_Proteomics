library(tidyverse)
library(ggrepel)
library(patchwork)
library(lcmm)
library(stringr)
library(splines)
library(finalfit)
library(ggpubr)
library(qvcalc)
setwd("E:/cotton")
rm(list = ls())
source("codes/functions/1_functions_loadall.R")

load("results_knn/lcmmdata_fev1_class4_cov.RData")
# load("lcmmdata_fev1_class3.RData")
# load("LCMM.RData")
vec <- read.table("vec_protein.txt") %>% unlist()
lcmmdata_fev1 <- lcmmdata_fev1[,c(1:13,which(colnames(lcmmdata_fev1) %in% vec))]
lcmmdata_fev1 <- get_imputed(lcmmdata_fev1, method = "seqknn")

########## Protein ~ FEV1 模型拟合 ##########

lcmmdata_fev1$class <- relevel(lcmmdata_fev1$class, ref = "1")

results_lcmm <- data.frame()
vec_protein <- which(str_detect(colnames(lcmmdata_fev1), pattern = "protein"))

for (i in vec_protein) {
  ## 
  temp_data <- lcmmdata_fev1[!is.na(lcmmdata_fev1[,i]),]
  colnames(temp_data)[i] <- "protein"
  
  lcmm_model <- lm(protein ~ class + age + end_cumulation + male81 + cesyr + pkyrs_cumulation + fev81 + height,
                 data=temp_data)
  
  results_lcmm <- rbind(results_lcmm,
                      c(colnames(lcmmdata_fev1)[i],
                        coef(summary(lcmm_model))[2,c(1,2,4)],
                        coef(summary(lcmm_model))[3,c(1,2,4)],
                        coef(summary(lcmm_model))[4,c(1,2,4)]))
  
  print(colnames(lcmmdata_fev1)[i])
}

colnames(results_lcmm) <- c("Protein",
                          "beta_class2", "se_class2", "Pval_class2",
                          "beta_class3", "se_class3", "Pval_class3",
                          "beta_class4", "se_class4", "Pval_class4")

get_colnum <- function(dat_vec) {
  dat_vec_origin <- dat_vec[1:10]
  return(which(dat_vec_origin == dat_vec[11]))
}

results_lcmm_sig <- results_lcmm
results_lcmm_sig$Pval_min <- apply(results_lcmm[,c(4,7,10)],1, min)
results_lcmm_sig$Pval_min_col <- apply(results_lcmm_sig,1, get_colnum)
table(results_lcmm_sig$Pval_min_col)
results_lcmm_sig$FDR_min <- p.adjust(results_lcmm_sig$Pval_min, method = "BH")
results_lcmm_sig_ordered <- results_lcmm_sig[order(results_lcmm_sig$Pval_min),]

########## Volcano Plot ##########
results_lcmm_volcano_class2 <- results_lcmm[order(results_lcmm$Pval_class2),]
colnames(results_lcmm_volcano_class2)[c(2,4)] <- c("beta", "Pval")
results_lcmm_volcano_class2$beta <- unlist(as.numeric(results_lcmm_volcano_class2$beta))
results_lcmm_volcano_class2$Pval <- unlist(as.numeric(results_lcmm_volcano_class2$Pval))
results_lcmm_volcano_class2$change <- ifelse(results_lcmm_volcano_class2$Pval < 0.05, 
                                      ifelse(results_lcmm_volcano_class2$beta > 0 ,'Up','Down'),
                                      'Nosig')

proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
colnames(proteinlist)[1] <- "Protein_ID"
results_lcmm_volcano_class2$gene.name <- ""
results_lcmm_volcano_class2[1:6,"gene.name"] <- apply(results_lcmm_volcano_class2[1:6,],1,match_gene_name)
# results_lcmm_volcano[which(results_lcmm_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

lcmm_plot_class2 <- volcano_plot(results_lcmm_volcano_class2,
                                 titles = "A Class2",
                                 x.limits = c(-3,3),
                                 x.breaks = seq(-3,3,1),
                                 y.limits = c(0,4),
                                 y.breaks = seq(0,4,0.5),
                                 gene_name = T)

results_lcmm_volcano_class3 <- results_lcmm[order(results_lcmm$Pval_class3),]
colnames(results_lcmm_volcano_class3)[c(5,7)] <- c("beta", "Pval")
results_lcmm_volcano_class3$beta <- unlist(as.numeric(results_lcmm_volcano_class3$beta))
results_lcmm_volcano_class3$Pval <- unlist(as.numeric(results_lcmm_volcano_class3$Pval))
results_lcmm_volcano_class3$change <- ifelse(results_lcmm_volcano_class3$Pval < 0.05, 
                                             ifelse(results_lcmm_volcano_class3$beta > 0 ,'Up','Down'),
                                             'Nosig')

results_lcmm_volcano_class3$gene.name <- ""
results_lcmm_volcano_class3[1:5,"gene.name"] <- apply(results_lcmm_volcano_class3[1:5,],1,match_gene_name)
# results_lcmm_volcano[which(results_lcmm_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

lcmm_plot_class3 <- volcano_plot(results_lcmm_volcano_class3,
                                 titles = "B Class3",
                                 x.limits = c(-3,3),
                                 x.breaks = seq(-3,3,1),
                                 y.limits = c(0,4),
                                 y.breaks = seq(0,4,0.5),
                                 gene_name = T)

results_lcmm_volcano_class4 <- results_lcmm[order(results_lcmm$Pval_class4),]
colnames(results_lcmm_volcano_class4)[c(8,10)] <- c("beta", "Pval")
results_lcmm_volcano_class4$beta <- unlist(as.numeric(results_lcmm_volcano_class4$beta))
results_lcmm_volcano_class4$Pval <- unlist(as.numeric(results_lcmm_volcano_class4$Pval))
results_lcmm_volcano_class4$change <- ifelse(results_lcmm_volcano_class4$Pval < 0.05, 
                                             ifelse(results_lcmm_volcano_class4$beta > 0 ,'Up','Down'),
                                             'Nosig')

results_lcmm_volcano_class4$gene.name <- ""
results_lcmm_volcano_class4[1:5,"gene.name"] <- apply(results_lcmm_volcano_class4[1:5,],1,match_gene_name)
# results_lcmm_volcano[which(results_lcmm_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

lcmm_plot_class4 <- volcano_plot(results_lcmm_volcano_class4,
                                 titles = "C Class4",
                                 x.limits = c(-3,3),
                                 x.breaks = seq(-3,3,1),
                                 y.limits = c(0,4),
                                 y.breaks = seq(0,4,0.5),
                                 gene_name = T)

lcmm_plot <- lcmm_plot_class2 + lcmm_plot_class3 + lcmm_plot_class4

jpeg("lcmm_plot.jpg", width = 36, height = 12, units = "cm",res = 1000,family="serif")
lcmm_plot
dev.off()

save(list = ls(), file = "LCMM.RData")
save(results_lcmm, file = "results_lcmm.RData")
