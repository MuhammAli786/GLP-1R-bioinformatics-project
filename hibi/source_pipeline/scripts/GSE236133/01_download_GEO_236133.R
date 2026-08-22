# GSE236133 download GEO metadata and supplementary files; input: GEO accession -> output: data/GSE236133_GEO_metadata.rds, metadata/GSE236133_raw_metadata.csv

library(GEOquery)
library(Biobase)

# Create folders
dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create(
  "data/GSE236133_supplementary",
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create("metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Download GEO Series metadata
cat("Downloading GSE236133 metadata...\n")

gse <- getGEO(
  "GSE236133",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

cat("Number of ExpressionSet objects:", length(gse), "\n")

eset <- gse[[1]]

cat("Platform:", annotation(eset), "\n")
cat(
  "Expression matrix dimensions:",
  paste(dim(exprs(eset)), collapse = " x "),
  "\n"
)

# Save GEO object
saveRDS(
  gse,
  file = "data/GSE236133_GEO_metadata.rds"
)

# Save sample metadata
sample_metadata <- pData(eset)

write.csv(
  sample_metadata,
  file = "metadata/GSE236133_raw_metadata.csv",
  row.names = TRUE
)

cat(
  "Saved metadata for",
  nrow(sample_metadata),
  "samples.\n"
)

# Download supplementary files
cat("\nDownloading supplementary files...\n")

supplementary_files <- getGEOSuppFiles(
  GEO = "GSE236133",
  makeDirectory = FALSE,
  baseDir = "data/GSE236133_supplementary"
)

cat("\nSupplementary-file download results:\n")
print(supplementary_files)

# List downloaded files
downloaded_files <- list.files(
  "data/GSE236133_supplementary",
  full.names = TRUE,
  recursive = TRUE
)

cat("\nDownloaded files:\n")
print(downloaded_files)

writeLines(
  downloaded_files,
  con = "metadata/GSE236133_downloaded_files.txt"
)

cat("\nDownload completed.\n")