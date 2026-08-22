# GSE238220 prepare metadata and raw read-count matrix; input: metadata/GSE238220_raw_metadata.csv, data/GSE238220_supplementary/GSE238220_read_counts.tsv.gz -> output: data/GSE238220_prepared_data.rds, metadata/GSE238220_prepared_metadata.csv

library(readr)
library(dplyr)
library(stringr)
library(tibble)

# File paths
metadata_file <- "metadata/GSE238220_raw_metadata.csv"

counts_file <- paste0(
  "data/GSE238220_supplementary/",
  "GSE238220_read_counts.tsv.gz"
)

output_file <- "data/GSE238220_prepared_data.rds"

output_metadata_file <- paste0(
  "metadata/",
  "GSE238220_prepared_metadata.csv"
)

# Check input files
if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found:\n",
    metadata_file,
    "\nRun Script 01 first."
  )
}

if (!file.exists(counts_file)) {
  stop(
    "Count file not found:\n",
    counts_file,
    "\nInspect the supplementary-file folder."
  )
}

# Read metadata and counts
metadata_raw <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts_raw <- read_tsv(
  counts_file,
  show_col_types = FALSE
)

cat("\nRaw metadata dimensions:\n")
print(dim(metadata_raw))

cat("\nRaw count-table dimensions:\n")
print(dim(counts_raw))

cat("\nCount-table columns:\n")
print(colnames(counts_raw))

# Identify GEO characteristics columns
characteristic_columns <- grep(
  "^characteristics",
  colnames(metadata_raw),
  value = TRUE,
  ignore.case = TRUE
)

if (length(characteristic_columns) == 0) {
  stop("No GEO characteristics columns were found.")
}

cat("\nCharacteristics columns:\n")
print(characteristic_columns)

# Function to extract one characteristic
get_characteristic <- function(data, pattern) {
  
  apply(
    data[, characteristic_columns, drop = FALSE],
    1,
    function(row_values) {
      
      row_values <- as.character(row_values)
      
      matched_values <- row_values[
        str_detect(
          row_values,
          regex(
            pattern,
            ignore_case = TRUE
          )
        )
      ]
      
      if (length(matched_values) == 0) {
        return(NA_character_)
      }
      
      str_trim(
        str_remove(
          matched_values[1],
          "^[^:]+:\\s*"
        )
      )
    }
  )
}

# Build clean metadata
metadata <- tibble(
  sample_id = metadata_raw$geo_accession,
  title = metadata_raw$title,
  
  tissue_raw = get_characteristic(
    metadata_raw,
    "^tissue:"
  ),
  
  cell_type_raw = get_characteristic(
    metadata_raw,
    "^cell type:"
  ),
  
  pup_id_raw = get_characteristic(
    metadata_raw,
    "^pup id:"
  ),
  
  sex_raw = get_characteristic(
    metadata_raw,
    "^sex:"
  ),
  
  time_raw = get_characteristic(
    metadata_raw,
    "time after hypoxia"
  ),
  
  treatment_raw = get_characteristic(
    metadata_raw,
    "^treatment:"
  )
) %>%
  mutate(
    tissue = case_when(
      str_detect(
        tissue_raw,
        regex("hippocampus", ignore_case = TRUE)
      ) ~ "Hippocampus",
      
      TRUE ~ tissue_raw
    ),
    
    cell_type = case_when(
      str_detect(
        cell_type_raw,
        regex("microglia", ignore_case = TRUE)
      ) ~ "Microglia",
      
      str_detect(
        cell_type_raw,
        regex("macrophage", ignore_case = TRUE)
      ) ~ "Macrophage",
      
      TRUE ~ cell_type_raw
    ),
    
    sex = case_when(
      str_detect(
        sex_raw,
        regex("^female$|^f$", ignore_case = TRUE)
      ) ~ "Female",
      
      str_detect(
        sex_raw,
        regex("^male$|^m$", ignore_case = TRUE)
      ) ~ "Male",
      
      TRUE ~ sex_raw
    ),
    
    time = case_when(
      str_detect(
        time_raw,
        regex("^1 day$|1 day|1d|24 hour", ignore_case = TRUE)
      ) ~ "1d",
      
      str_detect(
        time_raw,
        regex("^3 days$|3 day|3d|72 hour", ignore_case = TRUE)
      ) ~ "3d",
      
      TRUE ~ time_raw
    ),
    
    treatment = case_when(
      str_detect(
        treatment_raw,
        regex(
          "hypoxia|hypoxic.?ischemic|\\bHI\\b",
          ignore_case = TRUE
        )
      ) ~ "HI",
      
      str_detect(
        treatment_raw,
        regex(
          "^none$|control|untreated|normal|naive|sham",
          ignore_case = TRUE
        )
      ) ~ "Control",
      
      TRUE ~ treatment_raw
    ),
    
    pup_id = pup_id_raw
  ) %>%
  select(
    sample_id,
    title,
    tissue,
    cell_type,
    pup_id,
    sex,
    time,
    treatment,
    tissue_raw,
    cell_type_raw,
    pup_id_raw,
    sex_raw,
    time_raw,
    treatment_raw
  )

# Validate parsed metadata
cat("\nPrepared metadata:\n")
print(metadata)

cat("\nCell type x treatment x time:\n")

print(
  with(
    metadata,
    table(
      cell_type,
      treatment,
      time,
      useNA = "ifany"
    )
  )
)

cat("\nSex x treatment x time:\n")

print(
  with(
    metadata,
    table(
      sex,
      treatment,
      time,
      useNA = "ifany"
    )
  )
)

if (anyNA(metadata$sample_id)) {
  stop("One or more sample IDs are missing.")
}

if (anyDuplicated(metadata$sample_id)) {
  stop("Duplicate sample IDs were detected.")
}

if (anyNA(metadata$cell_type)) {
  warning("One or more cell types could not be parsed.")
}

if (anyNA(metadata$treatment)) {
  warning("One or more treatment values could not be parsed.")
}

if (anyNA(metadata$time)) {
  warning("One or more time values could not be parsed.")
}

# Identify gene-ID column in count table
possible_gene_columns <- c(
  "gene_id",
  "gene",
  "Gene",
  "Geneid",
  "gene_name",
  "symbol"
)

gene_column <- intersect(
  possible_gene_columns,
  colnames(counts_raw)
)

if (length(gene_column) == 0) {
  gene_column <- colnames(counts_raw)[1]
} else {
  gene_column <- gene_column[1]
}

cat("\nUsing gene column:\n")
print(gene_column)

# Match count columns to GEO samples
count_columns <- setdiff(
  colnames(counts_raw),
  gene_column
)

normalize_name <- function(x) {
  tolower(
    gsub(
      "[^a-zA-Z0-9]",
      "",
      x
    )
  )
}

normalized_count_columns <- normalize_name(
  count_columns
)

normalized_sample_ids <- normalize_name(
  metadata$sample_id
)

normalized_titles <- normalize_name(
  metadata$title
)

gsm_match <- match(
  normalized_sample_ids,
  normalized_count_columns
)

title_match <- match(
  normalized_titles,
  normalized_count_columns
)

matched_positions <- ifelse(
  !is.na(gsm_match),
  gsm_match,
  title_match
)

metadata$count_column <- count_columns[
  matched_positions
]

cat("\nSample-to-count-column matching:\n")

print(
  metadata %>%
    select(
      sample_id,
      title,
      cell_type,
      treatment,
      time,
      sex,
      count_column
    )
)

if (anyNA(metadata$count_column)) {
  
  cat("\nUnmatched samples:\n")
  
  print(
    metadata %>%
      filter(is.na(count_column)) %>%
      select(
        sample_id,
        title
      )
  )
  
  cat("\nAvailable count columns:\n")
  print(count_columns)
  
  stop(
    "One or more count columns could not be matched ",
    "to the metadata."
  )
}

if (anyDuplicated(metadata$count_column)) {
  stop("Duplicate count-column matches were detected.")
}

# Construct raw count matrix
count_table <- counts_raw %>%
  select(
    all_of(metadata$count_column)
  ) %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  )

count_matrix <- as.matrix(
  count_table
)

gene_ids <- as.character(
  counts_raw[[gene_column]]
)

valid_gene <- !is.na(gene_ids) &
  trimws(gene_ids) != ""

count_matrix <- count_matrix[
  valid_gene,
  ,
  drop = FALSE
]

gene_ids <- trimws(
  gene_ids[valid_gene]
)

rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- metadata$sample_id

# Collapse duplicate genes by sum
if (anyDuplicated(rownames(count_matrix))) {
  
  cat(
    "\nDuplicate gene identifiers detected.",
    "\nCollapsing duplicates by sum.\n"
  )
  
  count_matrix <- rowsum(
    count_matrix,
    group = rownames(count_matrix),
    reorder = FALSE,
    na.rm = TRUE
  )
}

# Validate count matrix
if (anyNA(count_matrix)) {
  stop("Missing count values were detected.")
}

if (any(count_matrix < 0)) {
  stop("Negative count values were detected.")
}

count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

if (!identical(
  colnames(count_matrix),
  metadata$sample_id
)) {
  stop(
    "Count matrix and metadata are not in the same order."
  )
}

cat("\nPrepared count-matrix dimensions:\n")
print(dim(count_matrix))

cat("\nCount-value summary:\n")
print(summary(as.vector(count_matrix)))

# Save prepared metadata
write.csv(
  metadata,
  output_metadata_file,
  row.names = FALSE
)

cat(
  "\nPrepared metadata saved to:\n",
  output_metadata_file,
  "\n"
)

# Save combined prepared object
prepared_data <- list(
  counts = count_matrix,
  metadata = metadata,
  gene_column = gene_column
)

saveRDS(
  prepared_data,
  output_file
)

cat(
  "\nPrepared data saved to:\n",
  output_file,
  "\n"
)

cat("\nScript 03 complete.\n")