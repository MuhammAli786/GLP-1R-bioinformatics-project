library(GEOquery)

gse <- readRDS("data/GSE23317.rds")
eset <- gse[[1]]

pheno <- pData(eset)

# Look at sample names
pheno$title

keep <- grepl("Sham|HI-", pheno$title)

table(keep)

pheno[keep, "title"]