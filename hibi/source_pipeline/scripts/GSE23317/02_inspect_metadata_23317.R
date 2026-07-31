# ==========================================
# 02_inspect_metadata.R
# Inspect sample metadata
# ==========================================

library(GEOquery)

# Load the downloaded dataset
gse <- readRDS("data/GSE23317.rds")

# If there are multiple platforms, use the first one
eset <- gse[[1]]

# Sample metadata
pheno <- pData(eset)

# View column names
colnames(pheno)

# View the first few rows
head(pheno)

# Dimensions
dim(pheno)

# Save metadata to CSV for easier viewing
write.csv(
  pheno,
  file = "metadata/GSE23317_metadata.csv",
  row.names = TRUE
)

cat("Metadata saved!\n")