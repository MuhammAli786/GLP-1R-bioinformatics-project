# GSE236133 inspect supplementary files and sample metadata; input: data/GSE236133_supplementary/, metadata/GSE236133_raw_metadata.csv -> output: metadata/GSE236133_supplementary_file_inventory.csv, metadata/GSE236133_metadata_inspection.csv

library(dplyr)

# List supplementary files
supplementary_directory <- "data/GSE236133_supplementary"

if (!dir.exists(supplementary_directory)) {
  stop("Supplementary directory not found. Run Script 01 first.")
}

files <- list.files(
  supplementary_directory,
  full.names = TRUE,
  recursive = TRUE
)

cat("\nSupplementary files:\n")

file_information <- data.frame(
  filename = basename(files),
  full_path = files,
  size_MB = round(file.info(files)$size / 1024^2, 3),
  stringsAsFactors = FALSE
)

print(file_information, row.names = FALSE)

write.csv(
  file_information,
  "metadata/GSE236133_supplementary_file_inventory.csv",
  row.names = FALSE
)

# Inspect sample metadata
metadata_file <- "metadata/GSE236133_raw_metadata.csv"

if (!file.exists(metadata_file)) {
  stop("Metadata file not found. Run Script 01 first.")
}

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata column names:\n")
print(colnames(metadata))

cat("\nSample titles:\n")
print(metadata$title)

cat("\nSource names:\n")
print(metadata$source_name_ch1)

# Inspect characteristics fields
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

# Show one row per sample
candidate_columns <- c(
  "title",
  "geo_accession",
  "source_name_ch1",
  characteristic_columns
)

# Keep only columns that genuinely exist
display_columns <- intersect(
  candidate_columns,
  colnames(metadata)
)

cat("\nColumns selected for sample-level display:\n")
print(display_columns)

if (length(display_columns) == 0) {
  warning("No requested display columns were found.")
} else {
  
  sample_metadata_display <- metadata[
    ,
    display_columns,
    drop = FALSE
  ]
  
  cat("\nSample-level metadata:\n")
  print(
    sample_metadata_display,
    row.names = FALSE
  )
  
  write.csv(
    sample_metadata_display,
    file = "metadata/GSE236133_metadata_inspection.csv",
    row.names = FALSE
  )
}

cat("\nInspection complete.\n")