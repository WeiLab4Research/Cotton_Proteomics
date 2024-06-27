library(dplyr)
library(stringr)
library(factoextra)
library(rms)
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

load("results_knn/lcmmdata_fev1_class4_cov.RData")
vec <- read.table("vec_protein.txt") %>% unlist()
lcmmdata_fev1 <- lcmmdata_fev1[,c(1:13,which(colnames(lcmmdata_fev1) %in% vec))]

lcmmdata_fev1 <- get_imputed(lcmmdata_fev1, method = "seqknn")
lcmmdata_fev1$slope <- slope

dd <- datadist(lcmmdata_fev1)
options(datadist="dd")

results_spline <- data.frame()
vec_protein <- which(str_detect(colnames(lcmmdata_fev1), pattern = "protein"))

for (i in vec_protein) {
  
  temp_lcmmdata_fev1 <- lcmmdata_fev1[!is.na(lcmmdata_fev1[,i]),]
  colnames(temp_lcmmdata_fev1)[i] <- "protein"
  
  spline_ols <- ols(protein ~ rcs(slope,3) + age + end_cumulation + male81 + cesyr + pkyrs_cumulation + fev81 + height,
                    data = temp_lcmmdata_fev1)
  
  an <- anova(spline_ols)
  Fval_linear <- (an[1,2] - an[2,2]) / an[11,3]
  Pval_linear <- 1 - pf(Fval_linear, 1, an[11,1])
  
  results_spline <- results_spline %>% rbind(c(colnames(lcmmdata_fev1)[i], coef(spline_ols)[2], an[1,5], coef(spline_ols)[3], an[2,5], Pval_linear))
  
  print(colnames(lcmmdata_fev1)[i])
}

colnames(results_spline) <- c("Protein", "beta_linear", "Pval_total", "beta_nonlinear", "Pval_nonlinear", "Pval_linear")

results_spline$Pval_min <- apply(results_spline[,5:6],1,min)

results_spline$FDR_min <- p.adjust(results_spline$Pval_min, method = "BH")
results_spline_ordered <- results_spline[order(results_spline$Pval_min),]

########## Volcano Plot ##########
results_spline_volcano <- results_spline[order(results_spline$Pval_linear),]
colnames(results_spline_volcano)[c(2,6)] <- c("beta", "Pval")
results_spline_volcano$beta <- unlist(as.numeric(results_spline_volcano$beta))
results_spline_volcano$Pval <- unlist(as.numeric(results_spline_volcano$Pval))
results_spline_volcano$change <- ifelse(results_spline_volcano$Pval < 0.05, 
                                    ifelse(results_spline_volcano$beta > 0 ,'Up','Down'),
                                    'Nosig')
proteinlist <- read.csv("original data/proteins_list.csv",row.names = 1)
colnames(proteinlist)[1] <- "Protein_ID"
results_spline_volcano$gene.name <- ""
results_spline_volcano[1:5,"gene.name"] <- apply(results_spline_volcano[1:5,],1,match_gene_name)
# results_spline_volcano[which(results_spline_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

spline_plot_linear <- volcano_plot(results_spline_volcano,
                                   titles = "A Spline-linear",
                                   x.limits = c(-0.2,0.2),
                                   x.breaks = seq(-0.2,0.2,0.05),
                                   y.limits = c(0,4),
                                   y.breaks = seq(0,4,0.5),
                                   gene_name = T)

results_spline_volcano_nonlinear <- results_spline[order(results_spline$Pval_nonlinear),]
colnames(results_spline_volcano_nonlinear)[c(4,5)] <- c("beta", "Pval")
results_spline_volcano_nonlinear$beta <- unlist(as.numeric(results_spline_volcano_nonlinear$beta))
results_spline_volcano_nonlinear$Pval <- unlist(as.numeric(results_spline_volcano_nonlinear$Pval))
results_spline_volcano_nonlinear$change <- ifelse(results_spline_volcano_nonlinear$Pval < 0.05, 
                                        ifelse(results_spline_volcano_nonlinear$beta > 0 ,'Up','Down'),
                                        'Nosig')

results_spline_volcano_nonlinear$gene.name <- ""
results_spline_volcano_nonlinear[1:5,"gene.name"] <- apply(results_spline_volcano_nonlinear[1:5,],1,match_gene_name)
# results_spline_volcano[which(results_spline_volcano$protein == "protein1221"), "gene.name"] <- "HBB"

spline_plot_nonlinear <- volcano_plot(results_spline_volcano_nonlinear,
                                      titles = "B Spline-nonlinear",
                                      x.limits = c(-0.2,0.2),
                                      x.breaks = seq(-0.2,0.2,0.05),
                                      y.limits = c(0,4),
                                      y.breaks = seq(0,4,0.5),
                                      gene_name = T)

spline_plot <- spline_plot_linear + spline_plot_nonlinear

jpeg("spline_plot.jpg", width = 24, height = 12, units = "cm",res = 1000,family="serif")
spline_plot
dev.off()

save(list = ls(), file = "Spline.RData")
save(results_spline,file = "results_spline.RData")

# Beta <- Predict(spline_ols, slope)
# plot(Beta,anova=an, pval=T)

