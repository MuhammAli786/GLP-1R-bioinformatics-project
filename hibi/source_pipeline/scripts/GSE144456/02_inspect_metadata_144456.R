# Inspect GSE144456 sample metadata and save raw metadata CSV
# Input: data/GSE144456.rds -> Output: metadata/GSE144456_metadata_raw.csv

library(GEOquery)

gse <- readRDS("data/GSE144456.rds")

cat("Number of ExpressionSet objects:", length(gse), "\n")

eset <- gse[[1]]
pheno <- pData(eset)

cat("\nPlatform:\n")
print(annotation(eset))

cat("\nExpression dimensions:\n")
print(dim(exprs(eset)))

cat("\nMetadata dimensions:\n")
print(dim(pheno))

cat("\nMetadata columns:\n")
print(colnames(pheno))

cat("\nSample titles:\n")
print(pheno$title)

cat("\nUnique values in each metadata column:\n")

for (column_name in colnames(pheno)) {
  
  values <- unique(as.character(pheno[[column_name]]))
  
  if (length(values) <= 20) {
    cat("\n==============================\n")
    cat(column_name, "\n")
    print(values)
  }
}

write.csv(
  pheno,
  "metadata/GSE144456_metadata_raw.csv",
  row.names = TRUE
)

cat("\nMetadata inspection complete.\n")