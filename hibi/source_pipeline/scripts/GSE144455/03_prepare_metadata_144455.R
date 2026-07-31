############################################################
# GSE144455 - Prepare Metadata
# Comparison: HI + PBS versus Naive + PBS at 3 hours
############################################################

library(GEOquery)

#-----------------------------------------------------------
# 1. Load the GEO object
#-----------------------------------------------------------

# Use the saved file if script 01 created it
if (file.exists("data/GSE144455.rds")) {
  
  gse <- readRDS("data/GSE144455.rds")
  
} else if (exists("gse")) {
  
  message("Using the GSE144455 object already in memory.")
  
} else {
  
  message("Saved GEO object not found. Downloading GSE144455.")
  gse <- getGEO("GSE144455", GSEMatrix = TRUE)
  
}

eset <- gse[[1]]
pheno <- pData(eset)

# Confirm sample order matches the expression matrix
stopifnot(
  identical(rownames(pheno), colnames(exprs(eset)))
)

#-----------------------------------------------------------
# 2. Inspect the channel-specific column names
#-----------------------------------------------------------

channel_columns <- grep(
  "time point|delay for observation|hypoxia|treatment|tissue",
  colnames(pheno),
  value = TRUE,
  ignore.case = TRUE
)

print(channel_columns)

#-----------------------------------------------------------
# 3. Helper function
#-----------------------------------------------------------

# Returns a metadata column when it exists.
# Otherwise, fills the column with NA.
get_metadata_column <- function(data, column_name) {
  
  if (column_name %in% colnames(data)) {
    as.character(data[[column_name]])
  } else {
    rep(NA_character_, nrow(data))
  }
  
}

#-----------------------------------------------------------
# 4. Construct a clean metadata table
#-----------------------------------------------------------

sample_metadata <- data.frame(
  
  sample_id = rownames(pheno),
  
  title = as.character(pheno$title),
  
  time_ch1 = get_metadata_column(
    pheno,
    "time point:ch1"
  ),
  
  time_ch2 = get_metadata_column(
    pheno,
    "delay for observation:ch2"
  ),
  
  injury_ch1 = get_metadata_column(
    pheno,
    "hypoxia-ischemia (rv) or none:ch1"
  ),
  
  injury_ch2 = get_metadata_column(
    pheno,
    "hypoxia-ischemia (rv) or none:ch2"
  ),
  
  treatment_ch1 = get_metadata_column(
    pheno,
    "treatment:ch1"
  ),
  
  treatment_ch2 = get_metadata_column(
    pheno,
    "treatment:ch2"
  ),
  
  stringsAsFactors = FALSE
)

# Some channel-1 values may be absent because channel 1 is consistently
# the naive/control reference. Retain the title for verification.
print(sample_metadata)

#-----------------------------------------------------------
# 5. Identify the HI versus naive PBS arrays at 3 hours
#-----------------------------------------------------------

# Titles for this comparison contain "RV PBS 3h".
# We also use channel-specific columns wherever available.

is_3h <- sample_metadata$time_ch2 == "3h" |
  grepl("PBS 3h", sample_metadata$title, ignore.case = TRUE)

is_HI <- sample_metadata$injury_ch2 == "Hypoxia-ischemia (RV)" |
  grepl("Hypoxia-ischemia vs naive", sample_metadata$title,
        ignore.case = TRUE)

is_PBS_ch1 <- sample_metadata$treatment_ch1 == "PBS" |
  grepl("RV PBS", sample_metadata$title, ignore.case = TRUE)

is_PBS_ch2 <- sample_metadata$treatment_ch2 == "PBS" |
  is.na(sample_metadata$treatment_ch2) |
  grepl("RV PBS", sample_metadata$title, ignore.case = TRUE)

sample_metadata$include_HI_vs_naive_PBS_3h <-
  is_3h & is_HI & is_PBS_ch1 & is_PBS_ch2

selected_metadata <- sample_metadata[
  sample_metadata$include_HI_vs_naive_PBS_3h,
  ,
  drop = FALSE
]

#-----------------------------------------------------------
# 6. Verify the selected arrays
#-----------------------------------------------------------

cat("\n====================================================\n")
cat("Selected arrays: HI + PBS versus Naive + PBS, 3 h\n")
cat("====================================================\n")

print(
  selected_metadata[, c(
    "sample_id",
    "title",
    "time_ch1",
    "time_ch2",
    "injury_ch1",
    "injury_ch2",
    "treatment_ch1",
    "treatment_ch2"
  )]
)

cat("\nNumber of selected arrays:",
    nrow(selected_metadata), "\n")

if (nrow(selected_metadata) == 0) {
  stop(
    "No 3-hour HI-versus-naive PBS arrays were selected. ",
    "Inspect the printed metadata column names and values."
  )
}

# The title list shown previously suggests there should be
# approximately three biological replicate arrays.
if (nrow(selected_metadata) < 2) {
  warning(
    "Fewer than two arrays were selected. Check the metadata carefully."
  )
}

#-----------------------------------------------------------
# 7. Save prepared metadata
#-----------------------------------------------------------

write.csv(
  sample_metadata,
  "metadata/GSE144455_metadata_prepared.csv",
  row.names = FALSE
)

write.csv(
  selected_metadata,
  "metadata/GSE144455_HI_vs_naive_PBS_3h_samples.csv",
  row.names = FALSE
)

saveRDS(
  selected_metadata,
  "metadata/GSE144455_HI_vs_naive_PBS_3h_samples.rds"
)

cat("\nMetadata preparation complete.\n")