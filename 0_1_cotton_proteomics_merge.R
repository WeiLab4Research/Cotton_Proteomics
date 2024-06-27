library(tidyverse)
library(factoextra)

#path of files
COVAR_PATH="cotton35_proteomics.csv"
PROTEIN_PATH="original data/report/data/Quantification/dia-proteinSummary-noqc.xls"
PEPTIDE_PATH="original data/report/data/Quantification/dia-peptideSummary-noqc.xls"

#output path
PROTEIN_OUT="cotton35_proteins.csv"
# PROTEIN_LIST_OUT="proteins_list.csv"
PEPTIDE_OUT="cotton35_peptide.csv"
# PEPTIDE_LIST_OUT="peptide_list.csv"

#load covariates data
cotton35_covar=read.csv(COVAR_PATH)

#load protein files
protein0=read.table(file = PROTEIN_PATH, sep = "\t", header=TRUE)
#list of each protein with their protein group and relative variable number in cotton_protein
proteinlist=protein0%>%select(Protein,ProteinGroup)%>%mutate(protein_seq=paste0("protein",row_number()))
protein=as.data.frame(t(protein0%>%select(-Protein,-ProteinGroup)))%>%rename_with(~ tolower(gsub("V", "protein", .x, fixed = TRUE)))
protein=protein%>%mutate(pid=rownames(protein))%>%remove_rownames()%>%relocate(pid)

#merge protein and covariates
cotton_protein=cotton35_covar%>%left_join(protein,by="pid")

#output protein
# write.csv(cotton_protein,PROTEIN_OUT)
# write.csv(proteinlist,PROTEIN_LIST_OUT)

#percent of missing by pid
cotton_protein%>%select(pid,starts_with("protein"))%>%mutate(missingness=rowMeans(is.na(.)))%>%select(pid,missingness)

#first 2 pcs stratified by cotton vs silk
pcs=prcomp(cotton_protein%>%select(starts_with("protein"))%>%mutate_all( ~replace_na(.,mean(., na.rm = TRUE))),center=TRUE,scale=TRUE)
fviz_eig(pcs)
fviz_pca_ind(pcs,
             col.ind=as.factor(cotton_protein$cotton),
             addEllipses = TRUE,
             ellipse.type="confidence",
             repel=TRUE)

#load peptide
peptide0=read.table(file = PEPTIDE_PATH, sep = "\t", header=TRUE)
#list of each peptide with their sequence and protein group as well as their relative variable number in cotton_peptide
peptidelist=peptide0%>%select(PeptideSequence,PrecursorCharge,ProteinGroup)%>%mutate(peptide_seq=paste0("peptide",row_number()))
peptide=as.data.frame(t(peptide0%>%select(-PeptideSequence,-PrecursorCharge,-ProteinGroup)))%>%rename_with(~ tolower(gsub("V", "peptide", .x, fixed = TRUE)))
peptide=peptide%>%mutate(pid=rownames(peptide))%>%remove_rownames()%>%relocate(pid)

#merge peptide and covariates
cotton_peptide=cotton35_covar%>%left_join(peptide,by="pid")

#output peptide
#write.csv(cotton_peptide,PEPTIDE_OUT)
# write.csv(peptidelist,PEPTIDE_LIST_OUT)