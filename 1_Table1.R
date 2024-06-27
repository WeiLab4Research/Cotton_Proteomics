rm(list = ls())
library(table1)
library(VIM)
library(tidyverse)
library(Amelia)
library(patchwork)
setwd("E:/cotton")
source("codes/functions/1_functions_loadall.R")

cotton_protein <- read.csv("original data/cotton35_proteomics.csv", header = T, row.names = 1)
cotton_protein <- cotton_protein[,-c(1,2)]
cotton_protein$protein <- T

cotton_no_protein <- read.csv("original data/update/cotton35_no_proteomics.csv", header = T, row.names = 1)
cotton_no_protein <- cotton_no_protein[,-c(1,157,158)]
cotton_no_protein$protein <- F
colnames(cotton_no_protein)[123:130] <- colnames(cotton_protein)[123:130]
cotton_total <- rbind(cotton_protein, cotton_no_protein)

cotton_total$fev81 <- as.numeric(cotton_total$fev81)
# proteins <- read.csv("original data/cotton35_proteins.csv", header = T, row.names = 1)
cotton_total$pfev81 <- cotton_total$fev81 / cotton_total$fvc81 * 100
cotton_total$pfev16 <- cotton_total$fev16 / cotton_total$fvc16 * 100
cotton_total$dfev <- (cotton_total$fev81 - cotton_total$fev16) / (cotton_total$age16 - cotton_total$age81)
cotton_total$pkyrs81[which(cotton_total$pkyrs81 == 0)] <- NA
cotton_total$pkyrs16[which(cotton_total$pkyrs16 == 0)] <- NA

cotton_total$cb16[which(is.na(cotton_total$cb16))] <- 0
cotton_total$chrcof16[which(is.na(cotton_total$chrcof16))] <- 0
cotton_total$byssin16[which(is.na(cotton_total$byssin16))] <- 0
cotton_total$dyspne16[which(is.na(cotton_total$dyspne16))] <- 0

cotton_total$resp16 <- as.character(as.numeric(cotton_total$cb16) | as.numeric(cotton_total$chrcof16) | as.numeric(cotton_total$byssin16) | as.numeric(cotton_total$dyspne16))

cotton_total$male81 <- factor(cotton_total$male81, levels = c(1,0), labels = c("Male","Female"))
cotton_total$smoker81 <- factor(cotton_total$smoker81, levels = c(0,1,2), labels = c("Non-smoker","Current","Former"))
cotton_total$smoker16 <- factor(cotton_total$smoker16, levels = c(0,1,2), labels = c("Non-smoker","Current","Former"))
cotton_total$cotton <- factor(cotton_total$cotton, levels = c(1,0), labels = c("Cotton workers","Silk workers"))

label(cotton_total$male81) <- "Sex"
label(cotton_total$age81) <- "Age"
label(cotton_total$age16) <- "Age"
label(cotton_total$smoker81) <- "Smoking status"
label(cotton_total$smoker16) <- "Smoking status"

follow <- c("81", "86", "92", "96", "01", "06", "11", "16")
for (i in follow) {
  label(cotton_total[,paste0("fev",i)]) <- paste0("FEV-1 ",i)
}
label(cotton_total$pfev81) <- "FEV-1/FVC"
label(cotton_total$pfev16) <- "FEV-1/FVC"
label(cotton_total$dfev) <- "Annual FEV-1 decline"
units(cotton_total$age81) <- "years"
units(cotton_total$age16) <- "years"
for (i in follow) {
  units(cotton_total[,paste0("fev",i)]) <- "ml"
}
units(cotton_total$fev81) <- "ml"
units(cotton_total$fev16) <- "ml"
units(cotton_total$dfev) <- "ml/year"
units(cotton_total$pfev81) <- "%"
units(cotton_total$pfev16) <- "%"

cotton_total[,77:108] <- apply(cotton_total[,77:108], 2, as.character)

table1(~ male81 + age16 + smoker16 + as.numeric(fev81) + fev16 + pfev16 + pkyrs16 + height16 + cesyr16 + end7 + fvc16 + ppfev_gli16 + dfev + cb16 + chrcof16 + byssin16 + dyspne16 + resp16| cotton,
       data = cotton_total[which(cotton_total$protein),], overall = F, topclass = "Rtable1-shade", extra.col = list(`统计量`=statistic,`P-value`=pvalue),
       render.continuous = my.render.cont, render.categorical = my.render.cat)

table1(~ male81 + age81 + smoker81 + fev81 + pfev81 + pkyrs81 + height81 + end0 + fvc81 + ppfev_gli81| cotton,
       data = cotton_total, overall = F, topclass = "Rtable1-shade", extra.col = list(`统计量`=statistic,`P-value`=pvalue),
       render.continuous = my.render.cont, render.categorical = my.render.cat)

summary(cotton_protein$cesyr16)
summary(cotton_protein$age16)
sd(cotton_protein$age16)
#绘制FEV1随访缺失情况
fev <- cotton_protein %>% select(starts_with("fev"))

aggr(fev, prob = F, numbers = F)

matrixplot(fev)

jpeg("results_present/figures/1_Table1/fev1_miss.jpg", width = 20, height = 12, units = "cm",res = 500,family="serif")
aggr(fev, prob = F, numbers = F)
dev.off()

jpeg("results_present/figures/1_Table1/protein_miss.jpg", width = 20, height = 12, units = "cm",res = 500,family="serif")
missmap(proteins[,159:3120])
dev.off()

#计算每个样本的FEV1重复测量次数
num_fev1 <- data.frame(cotton = cotton_protein$cotton)

num_fev1$num <- NA
for (i in 1:nrow(num_fev1)) {
  num_fev1[i,2] <- sum(!is.na(cotton_protein[i,109:116]))
}


## 计算蛋白质基因注释信息比例及缺失率
rm(list = ls())
proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
cotton_protein <- read.csv("original data/cotton35_proteins.csv", header = T, row.names = 1)

mmrmdata <- get_MMRMData(varvec_repeated = c("fev","fvc","age","pkyrs","end"),
                         varvec_single = c("male81", "fev81", "cesyr16", "height16"),
                         protein_limit = 0.5,
                         age_central = 60)

nrow(proteinlist[which(proteinlist$gene.name != ""),]) #2532
nrow(proteinlist[which(proteinlist$gene.name != "" & proteinlist$Protein %in% colnames(mmrmdata)),]) #907





