library(dplyr)
library(tidyr)
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
source("R script/functions/1_functions_loadall.R")

########## lcmm——拟合潜变量 ##########

##载入MMRM中使用的整理完成的原始数据
load("results_present/RData/MMRM_new.RData")

## 生成新的LCMM使用的数据格式
lcmmdata <- mmrmdata
#### lcmm()要求pid变量需为numeric型
lcmmdata$pid <- as.numeric(unlist(str_extract_all(lcmmdata$pid, "\\d+")))

latent_class_list_fev1 <- fit_latent_class(lcmmdata = lcmmdata, y = "fev", x = "age",
                                      random = 1,
                                      link = "linear",
                                      subject = "pid",
                                      max_ng = 8,
                                      maxiter = 2000,
                                      nproc = 5)

save(latent_class_list_fev1, lcmmdata, file = "cotton_lcmm8_fev1_cov.RData")

latent_class_list_fvc <- fit_latent_class(lcmmdata = lcmmdata, y = "fvc", x = "age",
                                           random = 1,
                                           link = "linear",
                                           subject = "pid",
                                           max_ng = 8,
                                           maxiter = 2000,
                                           nproc = 5)

save(latent_class_list_fvc, lcmmdata, file = "cotton_lcmm8_fvc.RData")
