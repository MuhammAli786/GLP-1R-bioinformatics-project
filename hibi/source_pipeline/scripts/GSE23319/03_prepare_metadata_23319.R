############################################################
## Prepare expression matrix
############################################################

gse <- readRDS("data/GSE23319.rds")

eset <- gse[[1]]

expr <- exprs(eset)

pheno <- pData(eset)

dim(expr)

head(pheno)