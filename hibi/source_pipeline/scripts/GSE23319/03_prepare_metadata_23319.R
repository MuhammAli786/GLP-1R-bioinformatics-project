# Prepare GSE23319 expression matrix and sample metadata; input: data/GSE23319.rds -> output: expr, pheno objects

gse <- readRDS("data/GSE23319.rds")

eset <- gse[[1]]

expr <- exprs(eset)

pheno <- pData(eset)

dim(expr)

head(pheno)