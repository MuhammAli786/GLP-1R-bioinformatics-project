library(GEOquery)

gse <- getGEO("GSE144455", GSEMatrix = TRUE)

length(gse)

pheno <- pData(gse[[1]])

colnames(pheno)

table(pheno$title)

saveRDS(
  gse,
  file = "data/GSE144455.rds"
)