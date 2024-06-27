library(tidyverse)
library(patchwork)
setwd("E:/cotton")

cotton_protein <- read.csv("cotton35_proteomics.csv", header = T, row.names = 1)
#数据整理
agelist <- cotton_protein %>% select(starts_with("age"))
agelist$pid <- cotton_protein$pid
agelist <- agelist %>% pivot_longer(cols = !pid,
                                    names_to = "agexx",
                                    values_to = "age")

#原始FEV1
fevlist <- cotton_protein %>% select(starts_with("fev"))
fevlist$pid <- cotton_protein$pid
fevlist <- fevlist %>% pivot_longer(cols = !pid,
                                    names_to = "fevxx",
                                    values_to = "fev")

#FEV1预测值百分比
ppfevlist <- cotton_protein %>% select(starts_with("ppfev"))
ppfevlist$pid <- cotton_protein$pid
ppfevlist <- ppfevlist %>% pivot_longer(cols = !pid,
                                        names_to = "ppfevxx",
                                        values_to = "ppfev")

cesdate <- cotton_protein %>% select(starts_with("cesyr"))
cesdate$pid <- cotton_protein$pid
cesdate <- cesdate %>% pivot_longer(cols = !pid,
                                    names_to = "cesyrxx",
                                    values_to = "cesdat")
cesdate$cessation <- cesdate$cesdat > 0

smoke <- cotton_protein %>% select(starts_with("cigevr"))
smoke[is.na(smoke)] <- 0
smoke$pid <- cotton_protein$pid
smoke <- smoke %>% pivot_longer(cols = !pid,
                                names_to = "cigevrxx",
                                values_to = "cigevr")

cotton <- rep(cotton_protein$cotton, each = 8)
fev_traj <- cbind(agelist[,c(1,3)], fevlist$fev, ppfevlist$ppfev, cotton, cesdate$cessation, smoke$cigevr)
colnames(fev_traj) <- c("pid", "age", "FEV1", "PPFEV1", "cotton", "Working Status", "cigevr")
fev_traj$cotton <- factor(fev_traj$cotton, levels = c(0, 1),
                          labels = c("Silk", "Cotton"))
fev_traj$cigevr <- factor(fev_traj$cigevr, levels = c(0, 1),
                          labels = c("Non-smoker", "Smoker"))
fev_traj$`Working Status` <- factor(fev_traj$`Working Status`,levels = c(F,T),
                                    labels = c("Active","Retired"))
#画图
fev1 <- ggplot(data = fev_traj, aes(x=age, y=FEV1, group=pid, color=`Working Status`))+
  geom_line(size=0.3, alpha = 0.5)+
  # stat_summary(fun = mean, geom = "line")+
  facet_grid(cigevr ~ cotton)
fev1

ppfev1 <- ggplot(data = fev_traj, aes(x=age, y=PPFEV1, group=pid, color=`Working Status`))+
  geom_line(size=0.3, alpha = 0.5)+
  # stat_summary(fun = mean, geom = "line")+
  facet_grid(cigevr ~ cotton)
ppfev1

jpeg("results_present/figures/2_FEV1_trajectory/fev1_traj.jpg", width = 20, height = 12, units = "cm",res = 500,family="serif")
fev1
dev.off()

jpeg("results_present/figures/2_FEV1_trajectory/ppfev1_traj.jpg", width = 20, height = 12, units = "cm",res = 500,family="serif")
ppfev1
dev.off()

