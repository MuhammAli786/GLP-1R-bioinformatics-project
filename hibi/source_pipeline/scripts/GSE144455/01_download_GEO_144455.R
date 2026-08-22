# Download GSE144455 GEO dataset and save as RDS
# Input: none (downloads from GEO) -> Output: data/GSE144455.rds

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