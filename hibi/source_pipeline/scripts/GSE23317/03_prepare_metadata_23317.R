# Prepare GSE23317 metadata: filter Sham and HI samples
# Input: data/GSE23317.rds -> Output: (prints filtered sample titles)

library(GEOquery)

gse <- readRDS("data/GSE23317.rds")
eset <- gse[[1]]

pheno <- pData(eset)

pheno$title

keep <- grepl("Sham|HI-", pheno$title)

table(keep)

pheno[keep, "title"]