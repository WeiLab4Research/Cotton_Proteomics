library(TwoSampleMR)
library(MendelianRandomization)
library(mr.raps)
library(MRPRESSO)
library(dplyr)
library(bigreadr)
setwd("/public/home/mszhao")

filelist <- unlist(read.csv("proteinlist.csv", header = F))

# gwas <- fread2("/public/home/mszhao/MR_cotton/20150_irnt.gwas.imputed_v3.both_sexes.tsv")
# variants <- fread2("/public/home/mszhao/MR_cotton/variants.tsv")
# gwas <- left_join(gwas,variants[,1:6])
# gwas$Phenotype <- "FEV1"
# rm(variants)

# gwas_clean <- gwas[which(gwas$pval <= 0.00000005),]
# gwas_clean <- format_data(gwas_clean,
# type = "exposure",
# phenotype_col = "Phenotype",
# snp_col = "rsid",
# beta_col = "beta",
# se_col = "se",
# eaf_col = "minor_AF",
# effect_allele_col = "ref",
# other_allele_col = "alt",
# pval_col = "pval")
# gwas_clean <- clump_data(gwas_clean, clump_r2 = 0.001, clump_kb = 10000)

gwas_clean <- fread2("MR_cotton/MR_extra/gwas_clean.csv")

for (name in filelist) {
  
  pQTL <- fread2(paste0(name,".tsv"))
  pQTL_clean <- pQTL[which(pQTL$variant_id %in% gwas_clean$SNP),]
  pQTL_clean$Phenotype <- name
  pQTL_clean <- format_data(pQTL_clean,
                            type = "outcome",
                            phenotype_col = "Phenotype",
                            snp_col = "variant_id",
                            beta_col = "beta",
                            se_col = "standard_error",
                            eaf_col = "effect_allele_frequency",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            pval_col = "p_value")
  
  
  dat <- harmonise_data(
    exposure_dat = gwas_clean,
    outcome_dat = pQTL_clean,
    action=1)
  mr <- mr(dat)
  
  mrpresso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
  
  
  save(list = c("dat", "mr", "mrpresso"), file = paste0("FEV1_",name,".RData"))
}

