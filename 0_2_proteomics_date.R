library(ggplot2)
setwd("E:/cotton")

proteomics_date <- read.csv("original data/proteomics_sample_date.csv")
proteomics_date$date <- as.Date(proteomics_date$date, format = "%m/%d/%Y")

ggplot(proteomics_date)+
  geom_histogram(aes(as.Date(date)), bins = 20)+
  scale_x_date(breaks = "6 month")

ggplot(proteomics_date)+
  geom_histogram(aes(age))

summary(proteomics_date$age)
