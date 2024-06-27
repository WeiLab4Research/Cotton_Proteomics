library(tidyverse)
setwd("E:/cotton")
rm(list = ls())

cotton_protein <- read.csv("original data/cotton35_proteins.csv", header = T, row.names = 1)
vec <- read.table("vec_protein.txt") %>% unlist()
cotton_protein <- cotton_protein[,c(1:158, which(colnames(cotton_protein) %in% vec))]
cotton_protein$diffev <- cotton_protein$fev16 - cotton_protein$fev81

results_simple <- data.frame()
vec_protein <- which(str_detect(colnames(cotton_protein), pattern = "protein"))
for (i in vec_protein) {
  temp_cotton_protein <- cotton_protein[!is.na(cotton_protein[,i]),]
  colnames(temp_cotton_protein)[i] <- "protein"
  
  simplemodel <- lm(protein ~ scale(diffev) + cotton + age16 + male81 + height16 + end6 + pkyrs16 + cesyr16 + fev81,
                    data = temp_cotton_protein)
  
  results_simple <- rbind(results_simple,
                          c(colnames(cotton_protein)[i],coef(summary(simplemodel))[2,c(1,2,4)]))
  print(colnames(cotton_protein)[i])
}
colnames(results_simple) <- c("Protein", "beta", "se", "Pval")
results_simple$beta <- unlist(as.numeric(results_simple$beta))
results_simple$Pval <- unlist(as.numeric(results_simple$Pval))
results_simple$FDR <- p.adjust(results_simple$Pval, method = "BH")
results_simple$change <- ifelse(results_simple$Pval < 0.05, 
                                    ifelse(results_simple$beta > 0 ,'Up','Down'),
                                    'Nosig')

volplot <- ggplot() +
  geom_point(aes(x = beta, y = -log10(Pval), color = change, alpha = 0.4),
             data = results_simple)+
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  labs(x="beta", y="-log10(pvalue)") +
  scale_color_manual(values=c("Nosig" = "#d2dae2",
                              "Up" = "#BC3C29FF",
                              "Down" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1,1),
                     breaks = seq(-1, 1, 0.2),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0,13),
                     breaks = seq(0, 10, 5),
                     expand = c(0.05,0.05)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold",
                                  size = 15),
        axis.text = element_text(size = 12))
volplot

jpeg("results_present/figures/0_4_simple_regression/simple_plot.jpg", width = 12, height = 10, units = "cm", res = 1000, family="serif")
volplot
dev.off()
