# Inspect GSE144455 sample metadata and save raw metadata CSV
# Input: data/GSE144455.rds -> Output: metadata/GSE144455_metadata.csv

library(GEOquery)

gse <- getGEO("GSE144455", GSEMatrix = TRUE)

length(gse)

eset <- gse[[1]]

dim(exprs(eset))

annotation(eset)

pheno <- pData(eset)

dim(pheno)
colnames(pheno)

table(pheno$title)

table(pheno$`time point:ch1`)
table(pheno$`hypoxia-ischemia (rv) or none:ch2`)
table(pheno$`treatment:ch1`)
table(pheno$`treatment:ch2`)

write.csv(
  pheno,
  "metadata/GSE144455_metadata.csv",
  row.names = TRUE
)