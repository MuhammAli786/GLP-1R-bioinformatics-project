# GSE238220 inspect sample metadata and read-count file; input: metadata/GSE238220_raw_metadata.csv, data/GSE238220_supplementary/GSE238220_read_counts.tsv.gz -> output: console inspection

library(readr)
library(dplyr)

metadata_file <- "metadata/GSE238220_raw_metadata.csv"

counts_file <- paste0(
  "data/GSE238220_supplementary/",
  "GSE238220_read_counts.tsv.gz"
)

if (!file.exists(metadata_file)) {
  stop("Metadata file missing. Run Script 01 first.")
}

if (!file.exists(counts_file)) {
  stop(
    "Read-count file missing:\n",
    counts_file,
    "\nInspect data/GSE238220_supplementary."
  )
}

# Metadata inspection
metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata columns:\n")
print(colnames(metadata))

cat("\nSample titles:\n")
print(metadata[, c("geo_accession", "title")])

characteristic_columns <- grep(
  "^characteristics",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nCharacteristics columns:\n")
print(characteristic_columns)

for (column in characteristic_columns) {
  
  cat("\n========================================\n")
  cat(column, "\n")
  cat("========================================\n")
  
  print(
    table(
      metadata[[column]],
      useNA = "ifany"
    )
  )
}

# Count-file inspection
counts_preview <- read_tsv(
  counts_file,
  n_max = 6,
  show_col_types = FALSE
)

cat("\nCount-file column names:\n")
print(colnames(counts_preview))

cat("\nCount-file preview:\n")
print(counts_preview)

cat("\nScript 02 complete.\n")