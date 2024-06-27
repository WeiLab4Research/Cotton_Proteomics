library(TwoSampleMR)
library(MendelianRandomization)
library(mr.raps)
library(MRPRESSO)
library(dplyr)
library(bigreadr)
library(patchwork)
setwd("E:/cotton")
rm(list = ls())
source("codes/functions/1_functions_loadall.R")

load("GCST90088811_buildGRCh37.RData")

mr <- mr[-c(2,4,5),]
mr[3:4,] <- mr[1,]
mr[3,c(5,7:9)] <- c("MRPRESSO",mrpresso[[1]][1,c(3,4,6)])
mr[4,c(5,7:9)] <- c("GSMR",gsmr_results$bxy,gsmr_results$bxy_se,gsmr_results$bxy_pval)

mr$b <- unlist(as.numeric(mr$b))
mr$se <- unlist(as.numeric(mr$se))
mr$pval <- unlist(as.numeric(mr$pval))

# mr <- mr[-1,]

hbb_fev1_plot <- mr_scatter_plot(mr, dat, xlab = "HBB",ylab = bquote(FEV[1]))

jpeg("hbb_fev1_plot.jpg", width = 25, height = 20, units = "cm",res = 1000,family="serif")
hbb_fev1_plot
dev.off()
####################
rm(list = ls())
source("codes/functions/1_functions_loadall.R")

load("FEV1_GCST90088811_buildGRCh37.RData")

mr <- mr[-c(2,4,5),]
mr[3:4,] <- mr[1,]
mr[3,c(5,7:9)] <- c("MRPRESSO",mrpresso[[1]][1,c(3,4,6)])
mr[4,c(5,7:9)] <- c("GSMR",gsmr_results$bxy,gsmr_results$bxy_se,gsmr_results$bxy_pval)

mr$b <- unlist(as.numeric(mr$b))
mr$se <- unlist(as.numeric(mr$se))
mr$pval <- unlist(as.numeric(mr$pval))

# mr <- mr[-1,]

fev1_hbb_plot <- mr_scatter_plot(mr, dat, xlab = bquote(FEV[1]), ylab = "HBB")

jpeg("fev1_hbb_plot.jpg", width = 25, height = 20, units = "cm",res = 1000,family="serif")
fev1_hbb_plot
dev.off()


