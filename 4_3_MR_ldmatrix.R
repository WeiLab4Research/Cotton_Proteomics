library(bigreadr)
library(TwoSampleMR)
library(MendelianRandomization)
library(mr.raps)
library(MRPRESSO)
library(dplyr)
setwd("/public/home/mszhao")

filelist <- unlist(read.csv("proteinlist.csv", header = F))

snp_list <- vector()
for (name in filelist) {
  
  pQTL <- fread2(paste0(name,".tsv"))
  pQTL_clean <- pQTL[which(pQTL$p_value <= 5e-6),]
  pQTL_clean$Phenotype <- name
  pQTL_clean <- format_data(pQTL_clean,
                            type = "exposure",
                            phenotype_col = "Phenotype",
                            snp_col = "variant_id",
                            beta_col = "beta",
                            se_col = "standard_error",
                            eaf_col = "effect_allele_frequency",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            pval_col = "p_value")
  pQTL_clean <- clump_data(pQTL_clean, clump_r2 = 0.001, clump_kb = 10000)
  
  snp_list <- unique(c(snp_list,pQTL_clean$SNP))
}

gwas_clean <- fread2("MR_cotton/MR_extra/gwas_clean.csv")

snp_list <- unique(c(snp_list,gwas_clean$rsid))

write.table(snp_list, file = "snp_list.txt", row.names = F, col.names = F, quote = F)

system("/data1/soft/plink --bfile /data1/user/canju/g1000_eur/g1000_eur --extract snp_list.txt --make-bed --out g1000_mr")
system("/data1/soft/plink --bfile g1000_mr --r square")

ld <- read.table("plink.ld", header = F)

snp_list <- read.table("snp_list.txt")
snp_list <- unlist(snp_list)
rownames(ld) = colnames(ld) = snp_list

write.table(ld, "snp_matrix.txt", quote = F)

write.table(snp_matrix, file = "snp_matrix.txt", quote = F)

