# Download GSE23317 GEO dataset and supplementary files, save as RDS
# Input: none (downloads from GEO) -> Output: data/GSE23317.rds

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

library(GEOquery)

if (!dir.exists("data")) {
  dir.create("data")
}

getGEOSuppFiles(
  GEO = "GSE23317",
  baseDir = "data"
)

gse <- getGEO(
  "GSE23317",
  GSEMatrix = TRUE
)

print(gse)

saveRDS(gse, file = "data/GSE23317.rds")

cat("Download complete!\n")