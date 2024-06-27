library(dplyr)
library(ACAT)
# library(EmpiricalBrownsMethod)
library(corrplot)
library(ggplot2)
library(stringr)
setwd("E:/cotton")
rm(list = ls())
source("codes/functions/1_functions_loadall.R")

load("results_knn/RData/6_1_Cluster/results_km.RData")
load("results_knn/RData/5_LCMM/results_lcmm.RData")
load("results_knn/RData/6_2_Spline/results_spline.RData")
load("results_knn/RData/3_MMRM/results_mmrm.RData")

vec <- read.table("vec_protein.txt") %>% unlist()
results_all <- results_all[which(results_all$Protein %in% vec),]
results_lcmm <- results_lcmm[which(results_lcmm$Protein %in% vec),]
results_spline <- results_spline[which(results_spline$Protein %in% vec),]
results_km <- results_km[which(results_km$Protein %in% vec),]

results_combined <- data.frame(protein = results_all[,1]) %>%
  mutate(Pval_km = as.numeric(apply(results_km[,c(4,7,10)],1,min))) %>%
  mutate(Pval_lcmm = as.numeric(apply(results_lcmm[,c(4,7,10)],1,min))) %>%
  mutate(Pval_spline = as.numeric(results_spline$Pval_min)) %>%
  mutate(Pval_mmrm = apply(results_all[,c(4,12,20,28,35)],1,min))

row.names(results_combined) <- results_combined[,1]
results_combined <- results_combined[,-1]

ACAT_Pval <- ACAT(t(results_combined))
ACAT_FDR <- p.adjust(ACAT_Pval, method = "BH")

results_ACAT <- data.frame(Protein = results_all$Protein,
                           P_ACAT = ACAT_Pval,
                           FDR_ACAT = ACAT_FDR)

results_ACAT_ordered <- results_ACAT[order(results_ACAT$P_ACAT),]

proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)

results_ACAT_ordered_gene <- left_join(results_ACAT_ordered,proteinlist[,3:4])

write.csv(results_ACAT_ordered_gene,file = "results_ACAT_ordered_gene.csv", row.names = F)
# write.csv(results_ACAT_ordered_uniprot,file = "results_ACAT_ordered_uniprot.csv", row.names = F)

results_combined_tan <- results_combined %>%
  mutate(Pval_km = tan((0.5 - results_combined$Pval_km) * pi)) %>%
  mutate(Pval_lcmm = tan((0.5 - results_combined$Pval_lcmm) * pi)) %>%
  mutate(Pval_spline = tan((0.5 - results_combined$Pval_spline) * pi)) %>%
  mutate(Pval_mmrm = tan((0.5 - results_combined$Pval_mmrm) * pi))

cor_matrix <- cor(results_combined, method = "spearman")

corplot <- corrplot(corr = cor_matrix, type = "full", method = "color", diag = TRUE, addCoef.col = "gray60")

jpeg("cor_plot.jpg", width = 10, height = 10, units = "cm",res = 1000,family="serif")
corrplot(corr = cor_matrix, type = "full", method = "color", diag = TRUE, addCoef.col = "gray60")
dev.off()


## 标记IG和HB蛋白质
results_ACAT_ordered_gene$gene.name.mark <- NA
results_ACAT_ordered_gene$gene.name.mark[which(str_detect(results_ACAT_ordered_gene$gene.name,
                                                          pattern = "^HB|^hb"))] <- "Hemoglobin"
results_ACAT_ordered_gene$gene.name.mark[which(str_detect(results_ACAT_ordered_gene$gene.name,
                                                          pattern = "^IG|^Ig"))] <- "Immunoglobulin"

## 瀑布图
p <- ggplot(data = results_ACAT_ordered_gene[which(results_ACAT_ordered_gene$P_ACAT < 0.05),],
            aes(x = reorder(Protein,P_ACAT), y = (-log10(P_ACAT) - (-log10(9.45e-4))), fill = gene.name.mark, color = gene.name.mark)) +
  geom_bar(stat="identity", width=0.4, position = position_dodge(width=0.2)) +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold",angle=90),
        panel.background = element_blank()) +
  labs(list(title = "Waterfall plot for Pvalues of proteins", x = "Proteins", y = "-log10(Pvalue)"))
p






