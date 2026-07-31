# ==========================================
# 01_download_GEO.R
# Download GEO dataset
# ==========================================

# Install GEOquery if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

# Load package
library(GEOquery)

# Create data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Download GSE23317 supplementary files
getGEOSuppFiles(
  GEO = "GSE23317",
  baseDir = "data"
)

# Download the GEO series matrix
gse <- getGEO(
  "GSE23317",
  GSEMatrix = TRUE
)

# Display basic information
print(gse)

# Save the object for later use
saveRDS(gse, file = "data/GSE23317.rds")

cat("Download complete!\n")