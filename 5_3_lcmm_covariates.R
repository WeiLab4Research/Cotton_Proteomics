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
source("R script/functions/1_functions_loadall.R")

load("results_present/lcmmdata_fev1_class4_cov.RData")
# load("lcmmdata_fev1_class3.RData")
# load("lcmmdata_fvc_class3.RData")

dat = lcmmdata_fev1_class4$cotton_class_df
# dat = lcmmdata_fev1_class3$cotton_class_df
pid <- fixed_cov(dat$pid)
class <- fixed_cov(dat$class)
class <- as.factor(class)
male81 <- fixed_cov(dat$male81)
age <- fixed_cov(dat$age)
height <- fixed_cov(dat$height)
cesyr <- fixed_cov(dat$cesyr)
fev81 <- fixed_cov(dat$fev81)


end <- time_dep_cov(cov_vec = dat$end,
                    age_vec = dat$age,
                    cov_name = "end")

pkyrs <- time_dep_cov(cov_vec = dat$pkyrs,
                      age_vec = dat$age,
                      cov_name = "pkyrs")

lcmmdata_fev1 <- cbind(pid, class, male81, age, height, cesyr, fev81, end, pkyrs)

for (i in colnames(dat)[str_detect(colnames(dat), pattern = "protein")]) {
  lcmmdata_fev1 <- cbind(lcmmdata_fev1, protein = fixed_cov(dat[,i]))
}

colnames(lcmmdata_fev1)[str_detect(colnames(lcmmdata_fev1), pattern = "protein")] <- colnames(dat)[str_detect(colnames(dat), pattern = "protein")]

save(lcmmdata_fev1, file = "lcmmdata_fev1_class4_cov.RData")
# save(lcmmdata_fev1, file = "lcmmdata_fev1_class3.RData")
