# Inspect GSE23317 sample metadata and save raw metadata CSV
# Input: data/GSE23317.rds -> Output: metadata/GSE23317_metadata.csv

library(GEOquery)

gse <- readRDS("data/GSE23317.rds")

eset <- gse[[1]]

pheno <- pData(eset)

colnames(pheno)

head(pheno)

dim(pheno)

write.csv(
  pheno,
  file = "metadata/GSE23317_metadata.csv",
  row.names = TRUE
)

cat("Metadata saved!\n")