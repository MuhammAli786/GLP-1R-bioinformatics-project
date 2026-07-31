############################################################
# GSE238220
# Download GEO metadata and supplementary files
############################################################

library(GEOquery)
library(Biobase)

dir.create(
  "data/GSE238220_supplementary",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "metadata",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/GSE238220",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/GSE238220",
  recursive = TRUE,
  showWarnings = FALSE
)

#-----------------------------------------------------------
# Download GEO metadata
#-----------------------------------------------------------

cat("Downloading GSE238220 GEO metadata...\n")

gse <- getGEO(
  "GSE238220",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

cat("Number of ExpressionSet objects:", length(gse), "\n")

eset <- gse[[1]]

cat("Platform:", annotation(eset), "\n")
cat("Number of samples:", ncol(eset), "\n")
cat(
  "Expression dimensions:",
  paste(dim(exprs(eset)), collapse = " x "),
  "\n"
)

saveRDS(
  gse,
  "data/GSE238220_GEO_metadata.rds"
)

#-----------------------------------------------------------
# Save sample metadata
#-----------------------------------------------------------

sample_metadata <- pData(eset)

write.csv(
  sample_metadata,
  "metadata/GSE238220_raw_metadata.csv",
  row.names = FALSE
)

cat(
  "Saved metadata for",
  nrow(sample_metadata),
  "samples.\n"
)

#-----------------------------------------------------------
# Download supplementary files
#-----------------------------------------------------------

cat("\nDownloading supplementary files...\n")

supplementary_files <- getGEOSuppFiles(
  GEO = "GSE238220",
  makeDirectory = FALSE,
  baseDir = "data/GSE238220_supplementary"
)

print(supplementary_files)

cat("\nDownloaded files:\n")

downloaded_files <- list.files(
  "data/GSE238220_supplementary",
  full.names = TRUE,
  recursive = TRUE
)

print(downloaded_files)

writeLines(
  downloaded_files,
  "metadata/GSE238220_downloaded_files.txt"
)

cat("\nScript 01 complete.\n")