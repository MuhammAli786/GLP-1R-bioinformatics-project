# Prepare GSE144455 metadata: select HI+PBS vs Naive+PBS at 3h arrays
# Input: data/GSE144455.rds -> Output: metadata/GSE144455_metadata_prepared.csv, metadata/GSE144455_HI_vs_naive_PBS_3h_samples.csv/.rds

library(GEOquery)

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

stopifnot(
  identical(rownames(pheno), colnames(exprs(eset)))
)

channel_columns <- grep(
  "time point|delay for observation|hypoxia|treatment|tissue",
  colnames(pheno),
  value = TRUE,
  ignore.case = TRUE
)

print(channel_columns)

get_metadata_column <- function(data, column_name) {
  
  if (column_name %in% colnames(data)) {
    as.character(data[[column_name]])
  } else {
    rep(NA_character_, nrow(data))
  }
  
}

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

print(sample_metadata)

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

if (nrow(selected_metadata) < 2) {
  warning(
    "Fewer than two arrays were selected. Check the metadata carefully."
  )
}

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