############################################################
## Inspect metadata
############################################################

gse <- readRDS("data/GSE23319.rds")

eset <- gse[[1]]

pheno <- pData(eset)

View(pheno)

write.csv(
  pheno,
  "metadata/GSE23319_metadata.csv",
  row.names = FALSE
)

table(pheno$title)

annotation(eset)