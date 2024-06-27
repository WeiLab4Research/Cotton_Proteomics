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

load("results_present/cotton_lcmm8_fev1_cov.RData")
# load("cotton_lcmm8_fvc.RData")

# color_lcmm <- c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")


lcmmdata_fev1_class4 <- latentclass_data(dat = lcmmdata,
                                       latentclass = latent_class_list_fev1$m4)

lcmmdata_fev1_class4$cotton_class_df$age = lcmmdata_fev1_class4$cotton_class_df$age + 60
# lcmmdata_fev1_class3 <- latentclass_data(dat = lcmmdata,
                                         # latentclass = latent_class_list_fev1$m3)

# lcmmdata_fvc_class <- latentclass_data(dat = lcmmdata,
                                       # latentclass = latent_class_list_fvc$m3)

predY_plot(latentclass = latent_class_list_fev1$m4,
           y = "fev",
           var.time = "age",
           dat = lcmmdata,
           bty = "l")

# predY_plot(latentclass = latent_class_list_fvc$m3,
           # y = "fvc",
           # var.time = "age",
           # dat = lcmmdata,
           # bty = "l")

class_plot_fev1 <- class_plot(dat = lcmmdata_fev1_class4$cotton_class_df,
                              leng = lcmmdata_fev1_class4$leng,
                              nclass = 4,
                              xvar = "age", yvar = "fev",
                              xlab = "Age", ylab = bquote(FEV[1]),
                              xlim = c(20,90), ylim = c(1000,6000))

# class_plot_fvc <- class_plot(dat = lcmmdata_fvc_class$cotton_class_df,
                             # leng = lcmmdata_fvc_class$leng,
                             # nclass = 3,
                             # xvar = "age", yvar = "fvc",
                             # xlab = "Age", ylab = "FVC",
                             # xlim = c(-40,30), ylim = c(1000,6000))

jpeg("lcmm_class_plot.jpg", width = 32, height = 24, units = "cm",res = 1000,family="serif")
class_plot_fev1[[5]]
dev.off()

jpeg("all_class_plot.jpg", width = 16, height = 12, units = "cm",res = 1000,family="serif")
class_plot_fev1[[6]]
dev.off()

save(lcmmdata_fev1_class4, file = "lcmmdata_fev1_class4_cov.RData")
# save(lcmmdata_fev1_class3, file = "lcmmdata_fev1_class3.RData")
# save(lcmmdata_fvc_class3, file = "lcmmdata_fvc_class3.RData")
