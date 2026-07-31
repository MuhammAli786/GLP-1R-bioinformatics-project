############################################################
# GSE144456 - Prepare Metadata
# Comparison selected for downstream analysis:
# P5 hypoxia-ischemia versus control at 3 hours
############################################################

library(GEOquery)

#-----------------------------------------------------------
# 1. Load dataset
#-----------------------------------------------------------

gse_file <- "data/GSE144456.rds"

if (!file.exists(gse_file)) {
  stop(
    "GSE144456 RDS file not found: ",
    gse_file,
    "\nRun script 01 first."
  )
}

gse <- readRDS(gse_file)

eset <- gse[[1]]
pheno <- pData(eset)

# Confirm metadata rows match expression-matrix columns
stopifnot(
  identical(
    rownames(pheno),
    colnames(exprs(eset))
  )
)

#-----------------------------------------------------------
# 2. Build clean metadata table
#-----------------------------------------------------------

sample_metadata <- data.frame(
  
  sample_id = rownames(pheno),
  
  title = as.character(
    pheno$title
  ),
  
  condition_ch1 = as.character(
    pheno[["condition:ch1"]]
  ),
  
  condition_ch2 = as.character(
    pheno[["condition:ch2"]]
  ),
  
  tissue_ch1 = as.character(
    pheno[["sublocalization/tissue:ch1"]]
  ),
  
  tissue_ch2 = as.character(
    pheno[["sublocalization/tissue:ch2"]]
  ),
  
  stringsAsFactors = FALSE
)

#-----------------------------------------------------------
# 3. Identify true ischemic-versus-control arrays
#-----------------------------------------------------------

sample_metadata$include_HI_vs_control <-
  tolower(trimws(sample_metadata$condition_ch1)) == "control" &
  tolower(trimws(sample_metadata$condition_ch2)) == "ischemic"

selected_metadata <- sample_metadata[
  sample_metadata$include_HI_vs_control,
  ,
  drop = FALSE
]

cat("\n============================================\n")
cat("All ischemic-versus-control arrays\n")
cat("============================================\n")

print(
  selected_metadata[, c(
    "sample_id",
    "title",
    "condition_ch1",
    "condition_ch2"
  )]
)

cat(
  "\nNumber of ischemic-versus-control arrays:",
  nrow(selected_metadata),
  "\n"
)

if (nrow(selected_metadata) == 0) {
  stop(
    "No ischemic-versus-control arrays were selected. ",
    "Check condition:ch1 and condition:ch2."
  )
}

#-----------------------------------------------------------
# 4. Extract developmental age from title
#-----------------------------------------------------------

selected_metadata$age <- ifelse(
  grepl(
    "day 5|P5\\+",
    selected_metadata$title,
    ignore.case = TRUE
  ),
  "P5",
  ifelse(
    grepl(
      "day 10|P10\\+",
      selected_metadata$title,
      ignore.case = TRUE
    ),
    "P10",
    NA_character_
  )
)

#-----------------------------------------------------------
# 5. Extract post-injury time from title
#-----------------------------------------------------------

selected_metadata$time_point <- sub(
  ".*\\+([0-9]+h).*",
  "\\1",
  selected_metadata$title
)

# Replace failed matches with NA
failed_time_match <-
  !grepl(
    "\\+[0-9]+h",
    selected_metadata$title
  )

selected_metadata$time_point[
  failed_time_match
] <- NA_character_

#-----------------------------------------------------------
# 6. Create clean group labels
#-----------------------------------------------------------

selected_metadata$analysis_group <- paste(
  selected_metadata$age,
  selected_metadata$time_point,
  sep = "_"
)

cat("\n============================================\n")
cat("Arrays by age and time point\n")
cat("============================================\n")

print(
  table(
    selected_metadata$age,
    selected_metadata$time_point,
    useNA = "ifany"
  )
)

#-----------------------------------------------------------
# 7. Select P5 at 3 hours
#-----------------------------------------------------------

selected_P5_3h <- selected_metadata[
  selected_metadata$age == "P5" &
    selected_metadata$time_point == "3h",
  ,
  drop = FALSE
]

cat("\n============================================\n")
cat("Selected P5, 3-hour arrays\n")
cat("============================================\n")

print(
  selected_P5_3h[, c(
    "sample_id",
    "title",
    "condition_ch1",
    "condition_ch2",
    "age",
    "time_point"
  )]
)

cat(
  "\nNumber of selected P5 3-hour arrays:",
  nrow(selected_P5_3h),
  "\n"
)

if (nrow(selected_P5_3h) == 0) {
  stop(
    "No P5 3-hour arrays were selected. ",
    "Inspect the sample titles and extracted metadata."
  )
}

if (nrow(selected_P5_3h) != 3) {
  warning(
    "Expected 3 biological replicate arrays for P5 at 3 hours, ",
    "but selected ",
    nrow(selected_P5_3h),
    "."
  )
}

#-----------------------------------------------------------
# 8. Confirm selected samples exist in expression matrix
#-----------------------------------------------------------

if (!all(
  selected_P5_3h$sample_id %in%
  colnames(exprs(eset))
)) {
  stop(
    "One or more selected samples are missing ",
    "from the expression matrix."
  )
}

#-----------------------------------------------------------
# 9. Save prepared metadata
#-----------------------------------------------------------

write.csv(
  sample_metadata,
  "metadata/GSE144456_metadata_prepared_all_samples.csv",
  row.names = FALSE
)

write.csv(
  selected_metadata,
  "metadata/GSE144456_HI_vs_control_all_times.csv",
  row.names = FALSE
)

write.csv(
  selected_P5_3h,
  "metadata/GSE144456_P5_3h_samples.csv",
  row.names = FALSE
)

saveRDS(
  selected_metadata,
  "metadata/GSE144456_HI_vs_control_all_times.rds"
)

saveRDS(
  selected_P5_3h,
  "metadata/GSE144456_P5_3h_samples.rds"
)

cat("\nMetadata preparation complete.\n")