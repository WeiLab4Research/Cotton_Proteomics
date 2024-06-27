setwd("E:/cotton/original data/update")

cotton_protein <- read.csv("cotton35_proteomics_updated.csv")
cotton_no_protein <- read.csv("cotton35_no_proteomics.csv")

cotton_protein$dod <- as.Date(cotton_protein$dod, format = "%m/%d/%Y")
cotton_no_protein$dod <- as.Date(cotton_no_protein$dod, format = "%m/%d/%Y")

nrow(cotton_protein[which(cotton_protein$dod > as.Date("1900/1/1")),]) #14
nrow(cotton_no_protein[which(cotton_no_protein$dod > as.Date("1900/1/1")),]) #161

nrow(cotton_no_protein[which(cotton_no_protein$dod > as.Date("1900/1/1") & cotton_no_protein$in16),]) #6


