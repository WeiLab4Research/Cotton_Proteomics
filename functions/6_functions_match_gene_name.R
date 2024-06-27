match_gene_name <- function(vec) {
  gene_name <- proteinlist[which(proteinlist$Protein == vec["Protein"]),"gene.name"]
  return(gene_name)
}
