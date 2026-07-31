############################################################
# GSE236133
# Prepare metadata and normalized expression matrix
############################################################

library(readr)
library(dplyr)
library(stringr)
library(tibble)

#-----------------------------------------------------------
# 1. File paths
#-----------------------------------------------------------

counts_file <- paste0(
  "data/GSE236133_supplementary/",
  "GSE236133_genelist_full_normalized_counts.txt.gz"
)

output_rds <- "data/GSE236133_prepared_data.rds"

output_metadata <- paste0(
  "metadata/",
  "GSE236133_prepared_metadata.csv"
)

#-----------------------------------------------------------
# 2. Check input file
#-----------------------------------------------------------

if (!file.exists(counts_file)) {
  stop(
    "Normalized-count file not found:\n",
    counts_file,
    "\nRun Script 01 first."
  )
}

#-----------------------------------------------------------
# 3. Read normalized-count table
#-----------------------------------------------------------

cat("\nReading normalized-count file...\n")

counts_raw <- read_tsv(
  counts_file,
  show_col_types = FALSE
)

cat("\nCount-table dimensions:\n")
print(dim(counts_raw))

cat("\nCount-table columns:\n")
print(colnames(counts_raw))

#-----------------------------------------------------------
# 4. Confirm gene ID column
#-----------------------------------------------------------

if (!"gene_id" %in% colnames(counts_raw)) {
  stop("The expected gene_id column was not found.")
}

# All columns except gene_id are expression samples
sample_columns <- setdiff(
  colnames(counts_raw),
  "gene_id"
)

cat("\nNumber of expression samples:\n")
print(length(sample_columns))

if (length(sample_columns) != 24) {
  warning(
    "Expected 24 expression columns, but found ",
    length(sample_columns),
    "."
  )
}

#-----------------------------------------------------------
# 5. Build metadata directly from count-column names
#
# Examples:
# wt_con1_3h
# wt_ips2_6h
# neil1_con1_3h
# neil2_ips2_6h
#-----------------------------------------------------------

metadata <- tibble(
  count_column = sample_columns
) %>%
  mutate(
    genotype_code = str_extract(
      count_column,
      "^(wt|neil1|neil2)"
    ),
    
    region_code = str_extract(
      count_column,
      "(?<=_)(con|ips)(?=[0-9])"
    ),
    
    biological_replicate = str_extract(
      count_column,
      "(?<=(con|ips))[0-9]+"
    ),
    
    time = str_extract(
      count_column,
      "[0-9]+h$"
    ),
    
    genotype = case_when(
      genotype_code == "wt" ~ "WT",
      genotype_code == "neil1" ~ "NEIL1-KO",
      genotype_code == "neil2" ~ "NEIL2-KO",
      TRUE ~ NA_character_
    ),
    
    region = case_when(
      region_code == "con" ~ "contralateral",
      region_code == "ips" ~ "ipsilateral",
      TRUE ~ NA_character_
    ),
    
    tissue = "hippocampus"
  ) %>%
  select(
    count_column,
    genotype,
    region,
    biological_replicate,
    time,
    tissue
  )

#-----------------------------------------------------------
# 6. Validate metadata parsing
#-----------------------------------------------------------

cat("\nPrepared metadata:\n")
print(metadata)

if (anyNA(metadata$genotype)) {
  stop("One or more genotypes could not be parsed.")
}

if (anyNA(metadata$region)) {
  stop("One or more regions could not be parsed.")
}

if (anyNA(metadata$biological_replicate)) {
  stop("One or more replicate numbers could not be parsed.")
}

if (anyNA(metadata$time)) {
  stop("One or more time points could not be parsed.")
}

if (anyDuplicated(metadata$count_column)) {
  stop("Duplicate sample-column names were detected.")
}

cat("\nSamples by genotype:\n")
print(table(metadata$genotype))

cat("\nSamples by region:\n")
print(table(metadata$region))

cat("\nSamples by time:\n")
print(table(metadata$time))

cat("\nGenotype x region x time:\n")
print(
  with(
    metadata,
    table(
      genotype,
      region,
      time
    )
  )
)

#-----------------------------------------------------------
# 7. Construct expression matrix
#-----------------------------------------------------------

expression_table <- counts_raw %>%
  select(all_of(sample_columns)) %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  )

expression_matrix <- as.matrix(expression_table)

rownames(expression_matrix) <- counts_raw$gene_id
colnames(expression_matrix) <- sample_columns

#-----------------------------------------------------------
# 8. Validate expression matrix
#-----------------------------------------------------------

cat("\nExpression matrix dimensions:\n")
print(dim(expression_matrix))

cat("\nExpression-value summary:\n")
print(summary(as.vector(expression_matrix)))

cat("\nMissing expression values:\n")
print(sum(is.na(expression_matrix)))

cat("\nAre all values integers?\n")
print(
  all(
    expression_matrix == round(expression_matrix),
    na.rm = TRUE
  )
)

if (!identical(
  colnames(expression_matrix),
  metadata$count_column
)) {
  stop(
    "Expression matrix columns and metadata rows ",
    "are not aligned."
  )
}

if (anyDuplicated(rownames(expression_matrix))) {
  warning(
    "Duplicate gene IDs are present. ",
    "These may need to be handled before differential expression."
  )
}

#-----------------------------------------------------------
# 9. Save prepared metadata
#-----------------------------------------------------------

dir.create(
  "metadata",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  metadata,
  output_metadata,
  row.names = FALSE
)

cat(
  "\nPrepared metadata saved to:\n",
  output_metadata,
  "\n"
)

#-----------------------------------------------------------
# 10. Save combined prepared object
#-----------------------------------------------------------

prepared_data <- list(
  expression = expression_matrix,
  metadata = metadata,
  gene_id = counts_raw$gene_id
)

saveRDS(
  prepared_data,
  output_rds
)

cat(
  "\nPrepared expression object saved to:\n",
  output_rds,
  "\n"
)

cat("\nScript 03 complete.\n")