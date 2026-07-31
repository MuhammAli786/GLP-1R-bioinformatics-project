############################################################
# GSE144455 - Inspect Metadata
############################################################

library(GEOquery)

# Download GEO dataset
gse <- getGEO("GSE144455", GSEMatrix = TRUE)

# Check number of ExpressionSets
length(gse)

# Use first ExpressionSet
eset <- gse[[1]]

# Expression matrix dimensions
dim(exprs(eset))

# Platform
annotation(eset)

# Sample metadata
pheno <- pData(eset)

# Metadata overview
dim(pheno)
colnames(pheno)

# Sample titles
table(pheno$title)

# Useful metadata
table(pheno$`time point:ch1`)
table(pheno$`hypoxia-ischemia (rv) or none:ch2`)
table(pheno$`treatment:ch1`)
table(pheno$`treatment:ch2`)

# Save metadata for later
write.csv(
  pheno,
  "metadata/GSE144455_metadata.csv",
  row.names = TRUE
)