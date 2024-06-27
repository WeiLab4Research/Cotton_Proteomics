library(tidyverse)
library(stringr)
library(factoextra)
library(rms)
library(ggplot2)
library(patchwork)
setwd("E:/cotton")
rm(list = ls())
source("codes/functions/1_functions_loadall.R")

cotton_protein <- read.csv("original data/cotton35_proteins.csv", header = T, row.names = 1)

age <- select(cotton_protein, starts_with("age"))
fev1 <- select(cotton_protein, starts_with("fev"))
slope <- vector()
for (i in 1:nrow(age)) {
  age_fev1 <- data.frame(age = unlist(age[i,]),
                         fev1 = unlist(fev1[i,]))
  lm_age_fev1 <- lm(fev1 ~ age, data = age_fev1)
  slope[i] <- coef(lm_age_fev1)["age"]
}
fev1_duration <- cotton_protein$age16 - cotton_protein$age81

clust_data <- data.frame(fev1_baseline = cotton_protein$fev81,
                         fev1_slope = slope)

elbowplot <- fviz_nbclust(clust_data, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)#4

jpeg("Elbow_plot.jpg", width = 15, height = 10, units = "cm",res = 500,family="serif")
elbowplot
dev.off()


km_fev1 <- kmeans(clust_data, 4, nstart = 24)
print(km_fev1)

load("results_knn/lcmmdata_fev1_class4_cov.RData")
vec <- read.table("vec_protein.txt") %>% unlist()
lcmmdata_fev1 <- lcmmdata_fev1[,c(1:13,which(colnames(lcmmdata_fev1) %in% vec))]
lcmmdata_fev1 <- get_imputed(lcmmdata_fev1, method = "seqknn")
lcmmdata_fev1$cluster <- as.factor(km_fev1$cluster)
lcmmdata_fev1$cluster <- relevel(lcmmdata_fev1$cluster, ref = 3)
lcmmdata_fev1$age = lcmmdata_fev1$age + 60

cotton_protein <- read.csv("original data/cotton35_proteins.csv", header = T, row.names = 1)

valid_proteins <- select_proteins(cotton_protein, limit = 0.5)

valid_proteins_imputed <- get_imputed(valid_proteins, method = "seqknn")

mmrmdata <- get_MMRMData(varvec_repeated = c("fev","fvc","age","pkyrs","end"),
                         varvec_single = c("male81", "fev81", "cesyr16", "height16"),
                         age_central = 0)
mmrmdata$clust <- as.factor(rep(km_fev1$cluster,each = 8))
mmrmdata$slope <- rep(slope,each = 8)
cluster_plot <- ggplot() +
  geom_smooth(data = mmrmdata[which(mmrmdata$clust ==1),], aes(x = age,y = fev),
              method = "lm", se = T, formula = y ~ x, color = color_lcmm[4], fill = color_lcmm[4])+
  geom_smooth(data = mmrmdata[which(mmrmdata$clust ==4),], aes(x = age,y = fev),
              method = "lm", se = T, formula = y ~ x, color = color_lcmm[3], fill = color_lcmm[3])+
  geom_smooth(data = mmrmdata[which(mmrmdata$clust ==3),], aes(x = age,y = fev),
              method = "lm", se = T, formula = y ~ x, color = color_lcmm[2], fill = color_lcmm[2])+
  geom_smooth(data = mmrmdata[which(mmrmdata$clust ==2),], aes(x = age,y = fev),
              method = "lm", se = T, formula = y ~ x, color = color_lcmm[1], fill = color_lcmm[1])+
  ylab(bquote(FEV[1]))+
  theme_classic()+
  theme(panel.border = element_rect(fill = NA), 
        axis.line = element_line(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

jpeg("cluster_plot.jpg", width = 16, height = 12, units = "cm",res = 1000,family="serif")
cluster_plot
dev.off()


results_km <- data.frame()
vec_protein <- which(str_detect(colnames(lcmmdata_fev1), pattern = "protein"))
for (i in vec_protein) {
  ## 
  temp_data <- lcmmdata_fev1[!is.na(lcmmdata_fev1[,i]),]
  colnames(temp_data)[i] <- "protein"

  km_model <- lm(protein ~ cluster + age + end_cumulation + male81 + cesyr + pkyrs_cumulation + fev81,
                    data=temp_data)
  
  results_km <- rbind(results_km,
                       c(colnames(lcmmdata_fev1)[i],
                         coef(summary(km_model))[2,c(1,2,4)],
                         coef(summary(km_model))[3,c(1,2,4)],
                         coef(summary(km_model))[4,c(1,2,4)]))
  
  print(colnames(lcmmdata_fev1)[i])
}
colnames(results_km) <- c("Protein",
                          "beta_cluster4", "se_cluster4", "Pval_cluster4",
                          "beta_cluster3", "se_cluster3", "Pval_cluster3",
                          "beta_cluster1", "se_cluster1", "Pval_cluster1")

get_colnum <- function(dat_vec) {
  dat_vec_origin <- dat_vec[1:10]
  return(which(dat_vec_origin == dat_vec[11]))
}

results_km_sig <- results_km
results_km_sig$Pval_min <- apply(results_km[,c(4,7,10)],1, min)
results_km_sig$Pval_min_col <- apply(results_km_sig,1, get_colnum)
table(results_km_sig$Pval_min_col)
results_km_sig$FDR_min <- p.adjust(results_km_sig$Pval_min, method = "BH")
results_km_sig_ordered <- results_km_sig[order(results_km_sig$Pval_min),]


########## Volcano Plot ##########
results_km_volcano_4 <- results_km[order(results_km$Pval_cluster4),]
colnames(results_km_volcano_4)[c(2,4)] <- c("beta", "Pval")
results_km_volcano_4$beta <- unlist(as.numeric(results_km_volcano_4$beta))
results_km_volcano_4$Pval <- unlist(as.numeric(results_km_volcano_4$Pval))
results_km_volcano_4$change <- ifelse(results_km_volcano_4$Pval < 0.05, 
                                      ifelse(results_km_volcano_4$beta > 0 ,'Up','Down'),
                                      'Nosig')
proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
colnames(proteinlist)[1] <- "Protein_ID"
results_km_volcano_4$gene.name <- ""
results_km_volcano_4[1:5,"gene.name"] <- apply(results_km_volcano_4[1:5,],1,match_gene_name)
# results_km_volcano[which(results_km_volcano$protein == "protein1221"), "gene.name"] <- "HBB"


km_plot_4 <- volcano_plot(results_km_volcano_4,
                          titles = "C Cluster4",
                          x.limits = c(-3.2,3.2),
                          x.breaks = seq(-3,3,1),
                          y.limits = c(0,4),
                          y.breaks = seq(0,4,0.5),
                        gene_name = T)

results_km_volcano_3 <- results_km[order(results_km$Pval_cluster3),]
colnames(results_km_volcano_3)[c(5,7)] <- c("beta", "Pval")
results_km_volcano_3$beta <- unlist(as.numeric(results_km_volcano_3$beta))
results_km_volcano_3$Pval <- unlist(as.numeric(results_km_volcano_3$Pval))
results_km_volcano_3$change <- ifelse(results_km_volcano_3$Pval < 0.05, 
                                    ifelse(results_km_volcano_3$beta > 0 ,'Up','Down'),
                                    'Nosig')
results_km_volcano_3$gene.name <- ""
results_km_volcano_3[1:5,"gene.name"] <- apply(results_km_volcano_3[1:5,],1,match_gene_name)
# results_km_volcano[which(results_km_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

km_plot_3 <- volcano_plot(results_km_volcano_3,
                          titles = "B Cluster3",
                          x.limits = c(-3.2,3.2),
                          x.breaks = seq(-3,3,1),
                          y.limits = c(0,4),
                          y.breaks = seq(0,4,0.5),
                          gene_name = T)

results_km_volcano_1 <- results_km[order(results_km$Pval_cluster1),]
results_km_volcano_1 <- results_km_volcano_1[-4,]
colnames(results_km_volcano_1)[c(8,10)] <- c("beta", "Pval")
results_km_volcano_1$beta <- unlist(as.numeric(results_km_volcano_1$beta))
results_km_volcano_1$Pval <- unlist(as.numeric(results_km_volcano_1$Pval))
results_km_volcano_1$change <- ifelse(results_km_volcano_1$Pval < 0.05, 
                                    ifelse(results_km_volcano_1$beta > 0 ,'Up','Down'),
                                    'Nosig')
results_km_volcano_1$gene.name <- ""
results_km_volcano_1[1:5,"gene.name"] <- apply(results_km_volcano_1[1:5,],1,match_gene_name)
# results_km_volcano[which(results_km_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

km_plot_1 <- volcano_plot(results_km_volcano_1,
                          titles = "A Cluster1",
                          x.limits = c(-3.2,3.2),
                          x.breaks = seq(-3,3,1),
                          y.limits = c(0,4),
                          y.breaks = seq(0,4,0.5),
                          gene_name = T)

km_plot <- km_plot_1 + km_plot_3 + km_plot_4

jpeg("km_plot.jpg", width = 36, height = 12, units = "cm",res = 1000,family="serif")
km_plot
dev.off()

save(list = ls(), file = "KMeans.RData")
save(results_km, file = "results_km.RData")

