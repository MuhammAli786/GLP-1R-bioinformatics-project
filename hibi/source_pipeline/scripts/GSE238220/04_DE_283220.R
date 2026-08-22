# GSE238220
# Differential expression using DESeq2
#
# Comparisons:
#   Microglia 1d HI vs Control
#   Microglia 3d HI vs Control
#
# Positive log2FoldChange = higher in HI
# Negative log2FoldChange = higher in Control

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Paths and thresholds

input_file <- "data/GSE238220_prepared_data.rds"

output_directory <- "results/GSE238220/DE"

annotation_file <- "metadata/GSE238220_gene_annotation.csv"

fdr_threshold <- 0.05
logfc_threshold <- 0.2

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Prepared data file not found:\n",
    input_file,
    "\nRun Script 03 first."
  )
}

# Load prepared data and rebuild analysis metadata

prepared_data <- readRDS(
  "data/GSE238220_prepared_data.rds"
)

count_matrix <- prepared_data$counts

metadata <- as.data.frame(
  prepared_data$metadata,
  stringsAsFactors = FALSE
)

# Ensure optional raw columns exist
for (column_name in c(
  "title",
  "cell_type_raw",
  "treatment_raw",
  "time_raw",
  "sex_raw"
)) {
  
  if (!column_name %in% colnames(metadata)) {
    metadata[[column_name]] <- NA_character_
  }
}

# Combine all available text for robust parsing
metadata <- metadata %>%
  mutate(
    sample_id = trimws(as.character(sample_id)),
    title = trimws(as.character(title)),
    
    cell_text = paste(
      title,
      cell_type_raw,
      cell_type
    ),
    
    treatment_text = paste(
      title,
      treatment_raw,
      treatment
    ),
    
    time_text = paste(
      title,
      time_raw,
      time
    ),
    
    sex_text = paste(
      title,
      sex_raw,
      sex
    ),
    
    # Microglial titles contain "Mµg" or metadata says Microglia.
    # Macrophage titles contain "Mφ" or metadata says Macrophages.
    cell_type = case_when(
      grepl(
        "microglia|Mµg",
        cell_text,
        ignore.case = TRUE
      ) ~ "Microglia",
      
      grepl(
        "macrophage|Mφ",
        cell_text,
        ignore.case = TRUE
      ) ~ "Macrophage",
      
      TRUE ~ NA_character_
    ),
    
    # Nµg means untreated/normal microglia.
    # HIµg means HI microglia.
    treatment = case_when(
      grepl(
        "hypoxia|treatment:\\s*hypoxia|HIµg|HI",
        treatment_text,
        ignore.case = TRUE
      ) ~ "HI",
      
      grepl(
        "treatment:\\s*none|^none$|Nµg|control",
        treatment_text,
        ignore.case = TRUE
      ) ~ "Control",
      
      TRUE ~ NA_character_
    ),
    
    time = case_when(
      grepl(
        "1 day|1d|f1|m1",
        time_text,
        ignore.case = TRUE
      ) ~ "1d",
      
      grepl(
        "3 days|3 day|3d|f3|m3",
        time_text,
        ignore.case = TRUE
      ) ~ "3d",
      
      TRUE ~ NA_character_
    ),
    
    sex = case_when(
      grepl(
        "sex:\\s*female|_f[13]",
        sex_text,
        ignore.case = TRUE
      ) ~ "Female",
      
      grepl(
        "sex:\\s*male|_m[13]",
        sex_text,
        ignore.case = TRUE
      ) ~ "Male",
      
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    sample_id,
    title,
    cell_type,
    treatment,
    time,
    sex,
    everything(),
    -cell_text,
    -treatment_text,
    -time_text,
    -sex_text
  )

cat("\nReconstructed metadata groups:\n")

group_table <- with(
  metadata,
  table(
    cell_type,
    treatment,
    time,
    useNA = "ifany"
  )
)

print(group_table)

cat("\nReconstructed microglial samples:\n")

print(
  metadata %>%
    filter(cell_type == "Microglia") %>%
    select(
      sample_id,
      title,
      sex,
      treatment,
      time
    ),
  row.names = FALSE
)

# Hard validation before running DESeq2
expected_groups <- metadata %>%
  dplyr::filter(
    cell_type == "Microglia",
    treatment %in% c("Control", "HI"),
    time %in% c("1d", "3d")
  ) %>%
  dplyr::group_by(
    time,
    treatment
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

cat("\nMicroglial sample counts used for analysis:\n")
print(expected_groups)

required_combinations <- expand.grid(
  time = c("1d", "3d"),
  treatment = c("Control", "HI"),
  stringsAsFactors = FALSE
)

missing_combinations <- dplyr::anti_join(
  required_combinations,
  expected_groups,
  by = c("time", "treatment")
)

if (nrow(missing_combinations) > 0) {
  stop(
    "Metadata reconstruction failed for one or more groups.\n",
    paste(
      apply(
        missing_combinations,
        1,
        paste,
        collapse = " / "
      ),
      collapse = "\n"
    )
  )
}

# Clean gene identifiers

gene_ids <- rownames(count_matrix)

if (is.null(gene_ids)) {
  stop("The count matrix has no row names.")
}

gene_ids <- trimws(
  as.character(gene_ids)
)

valid_gene <- !is.na(gene_ids) &
  gene_ids != "" &
  gene_ids != "NA"

count_matrix <- count_matrix[
  valid_gene,
  ,
  drop = FALSE
]

gene_ids <- gene_ids[valid_gene]

# Detect identifier type
ensembl_fraction <- mean(
  grepl("^ENSMUSG", gene_ids)
)

entrez_fraction <- mean(
  grepl("^[0-9]+$", gene_ids)
)

if (ensembl_fraction > 0.5) {
  
  gene_keytype <- "ENSEMBL"
  
  gene_ids_clean <- sub(
    "\\.[0-9]+$",
    "",
    gene_ids
  )
  
} else if (entrez_fraction > 0.5) {
  
  gene_keytype <- "ENTREZID"
  gene_ids_clean <- gene_ids
  
} else {
  
  gene_keytype <- "SYMBOL"
  gene_ids_clean <- gene_ids
}

cat(
  "\nDetected gene identifier type:",
  gene_keytype,
  "\n"
)

rownames(count_matrix) <- gene_ids_clean

# Collapse duplicate identifiers by summing counts
if (anyDuplicated(rownames(count_matrix))) {
  
  cat(
    "\nDuplicate gene identifiers detected.",
    "\nCollapsing duplicate rows by sum.\n"
  )
  
  count_matrix <- rowsum(
    count_matrix,
    group = rownames(count_matrix),
    reorder = FALSE,
    na.rm = TRUE
  )
  
  storage.mode(count_matrix) <- "integer"
}

# Create gene annotation table

annotation <- tibble(
  input_id = rownames(count_matrix)
)

valid_keys <- intersect(
  annotation$input_id,
  AnnotationDbi::keys(
    org.Mm.eg.db,
    keytype = gene_keytype
  )
)

annotation$SYMBOL <- NA_character_
annotation$ENTREZID <- NA_character_
annotation$ENSEMBL <- NA_character_
annotation$GENENAME <- NA_character_

if (length(valid_keys) > 0) {
  
  symbol_map <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = valid_keys,
    keytype = gene_keytype,
    column = "SYMBOL",
    multiVals = "first"
  )
  
  entrez_map <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = valid_keys,
    keytype = gene_keytype,
    column = "ENTREZID",
    multiVals = "first"
  )
  
  ensembl_map <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = valid_keys,
    keytype = gene_keytype,
    column = "ENSEMBL",
    multiVals = "first"
  )
  
  genename_map <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = valid_keys,
    keytype = gene_keytype,
    column = "GENENAME",
    multiVals = "first"
  )
  
  match_index <- match(
    valid_keys,
    annotation$input_id
  )
  
  annotation$SYMBOL[match_index] <-
    unname(symbol_map[valid_keys])
  
  annotation$ENTREZID[match_index] <-
    unname(entrez_map[valid_keys])
  
  annotation$ENSEMBL[match_index] <-
    unname(ensembl_map[valid_keys])
  
  annotation$GENENAME[match_index] <-
    unname(genename_map[valid_keys])
}

if (gene_keytype == "SYMBOL") {
  annotation$SYMBOL[
    is.na(annotation$SYMBOL)
  ] <- annotation$input_id[
    is.na(annotation$SYMBOL)
  ]
}

if (gene_keytype == "ENTREZID") {
  annotation$ENTREZID[
    is.na(annotation$ENTREZID)
  ] <- annotation$input_id[
    is.na(annotation$ENTREZID)
  ]
}

if (gene_keytype == "ENSEMBL") {
  annotation$ENSEMBL[
    is.na(annotation$ENSEMBL)
  ] <- annotation$input_id[
    is.na(annotation$ENSEMBL)
  ]
}

write_csv(
  annotation,
  annotation_file
)

cat(
  "\nMapped Entrez IDs:",
  sum(!is.na(annotation$ENTREZID)),
  "\n"
)

# Function to run one comparison

run_deseq_comparison <- function(selected_time) {
  
  comparison_name <- paste0(
    "Microglia_",
    selected_time,
    "_HI_vs_Control"
  )
  
  cat("\n\n=================================================\n")
  cat("Running:", comparison_name, "\n")
  cat("=================================================\n")
  
  # Select microglia at the requested time
  
  sample_metadata <- metadata %>%
    filter(
      grepl(
        "^Microglia$",
        cell_type,
        ignore.case = TRUE
      ),
      grepl(
        paste0("^", selected_time, "$"),
        time,
        ignore.case = TRUE
      ),
      treatment %in% c(
        "Control",
        "HI"
      )
    )
  
  cat("\nSelected samples:\n")
  
  print(
    sample_metadata %>%
      select(
        sample_id,
        title,
        sex,
        treatment,
        time
      )
  )
  
  treatment_counts <- table(
    sample_metadata$treatment
  )
  
  cat("\nTreatment counts:\n")
  print(treatment_counts)
  
  cat("\nSex by treatment:\n")
  
  print(
    table(
      sample_metadata$sex,
      sample_metadata$treatment,
      useNA = "ifany"
    )
  )
  
  if (
    !"Control" %in% names(treatment_counts) ||
    !"HI" %in% names(treatment_counts)
  ) {
    stop(
      comparison_name,
      ": both treatment groups are required."
    )
  }
  
  if (
    treatment_counts[["Control"]] < 2 ||
    treatment_counts[["HI"]] < 2
  ) {
    stop(
      comparison_name,
      ": at least two samples per group are required."
    )
  }
  
  # Match count matrix to metadata
  
  selected_counts <- count_matrix[
    ,
    sample_metadata$sample_id,
    drop = FALSE
  ]
  
  sample_metadata <- sample_metadata[
    match(
      colnames(selected_counts),
      sample_metadata$sample_id
    ),
    ,
    drop = FALSE
  ]
  
  rownames(sample_metadata) <- sample_metadata$sample_id
  
  if (!identical(
    rownames(sample_metadata),
    colnames(selected_counts)
  )) {
    stop(
      comparison_name,
      ": sample metadata and count matrix are not aligned."
    )
  }
  
  # Filter low-count genes
  
  keep_gene <- rowSums(
    selected_counts >= 10
  ) >= 2
  
  selected_counts <- selected_counts[
    keep_gene,
    ,
    drop = FALSE
  ]
  
  cat(
    "\nGenes before filtering:",
    nrow(count_matrix),
    "\n"
  )
  
  cat(
    "Genes after filtering:",
    nrow(selected_counts),
    "\n"
  )
  
  if (nrow(selected_counts) == 0) {
    stop(
      comparison_name,
      ": no genes remained after filtering."
    )
  }
  
  # Prepare DESeq2 design
  
  sample_metadata$treatment <- factor(
    sample_metadata$treatment,
    levels = c(
      "Control",
      "HI"
    )
  )
  
  sample_metadata$sex <- factor(
    sample_metadata$sex
  )
  
  use_sex <- all(
    !is.na(sample_metadata$sex)
  ) &&
    length(unique(sample_metadata$sex)) > 1
  
  if (use_sex) {
    
    design_formula <- ~ sex + treatment
    
    design_matrix <- model.matrix(
      design_formula,
      data = sample_metadata
    )
    
    if (
      qr(design_matrix)$rank <
      ncol(design_matrix)
    ) {
      
      warning(
        comparison_name,
        ": sex and treatment are confounded. ",
        "Using design = ~ treatment."
      )
      
      design_formula <- ~ treatment
    }
    
  } else {
    
    design_formula <- ~ treatment
  }
  
  cat(
    "\nDESeq2 design:",
    paste(
      deparse(design_formula),
      collapse = ""
    ),
    "\n"
  )
  
  # Run DESeq2
  
  dds <- DESeqDataSetFromMatrix(
    countData = selected_counts,
    colData = sample_metadata,
    design = design_formula
  )
  
  dds <- DESeq(
    dds,
    quiet = FALSE
  )
  
  result_object <- results(
    dds,
    contrast = c(
      "treatment",
      "HI",
      "Control"
    ),
    alpha = fdr_threshold
  )
  
  # Format results
  
  result <- as.data.frame(
    result_object
  ) %>%
    rownames_to_column("input_id") %>%
    left_join(
      annotation,
      by = "input_id"
    ) %>%
    mutate(
      comparison = comparison_name,
      
      direction = case_when(
        !is.na(padj) &
          padj < fdr_threshold &
          log2FoldChange >= logfc_threshold ~ "Up",
        
        !is.na(padj) &
          padj < fdr_threshold &
          log2FoldChange <= -logfc_threshold ~ "Down",
        
        TRUE ~ "Not_significant"
      ),
      
      significant = direction != "Not_significant"
    ) %>%
    select(
      comparison,
      input_id,
      SYMBOL,
      ENTREZID,
      ENSEMBL,
      GENENAME,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      direction,
      significant
    ) %>%
    arrange(
      is.na(padj),
      padj,
      pvalue
    )
  
  deg <- result %>%
    filter(significant)
  
  # Save result tables
  
  write_csv(
    result,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_all_genes.csv"
      )
    )
  )
  
  write_csv(
    deg,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_DEGs.csv"
      )
    )
  )
  
  # Save gene-symbol lists
  
  write_lines(
    deg %>%
      filter(
        !is.na(SYMBOL),
        SYMBOL != ""
      ) %>%
      distinct(SYMBOL) %>%
      pull(SYMBOL),
    
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_gene_symbols.txt"
      )
    )
  )
  
  write_lines(
    deg %>%
      filter(
        direction == "Up",
        !is.na(SYMBOL),
        SYMBOL != ""
      ) %>%
      distinct(SYMBOL) %>%
      pull(SYMBOL),
    
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_up_gene_symbols.txt"
      )
    )
  )
  
  write_lines(
    deg %>%
      filter(
        direction == "Down",
        !is.na(SYMBOL),
        SYMBOL != ""
      ) %>%
      distinct(SYMBOL) %>%
      pull(SYMBOL),
    
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_down_gene_symbols.txt"
      )
    )
  )
  
  # Save Entrez ID list
  
  write_lines(
    deg %>%
      filter(
        !is.na(ENTREZID),
        ENTREZID != ""
      ) %>%
      distinct(ENTREZID) %>%
      pull(ENTREZID) %>%
      as.character(),
    
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_entrez_ids.txt"
      )
    )
  )
  
  # Save DESeq2 object and selected samples
  
  saveRDS(
    dds,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_dds.rds"
      )
    )
  )
  
  write_csv(
    sample_metadata %>%
      rownames_to_column(
        "count_matrix_sample"
      ),
    
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_selected_samples.csv"
      )
    )
  )
  
  # Summary
  
  summary_row <- tibble(
    comparison = comparison_name,
    
    control_samples = sum(
      sample_metadata$treatment == "Control"
    ),
    
    HI_samples = sum(
      sample_metadata$treatment == "HI"
    ),
    
    tested_genes = nrow(result),
    
    significant_DEGs = nrow(deg),
    
    upregulated = sum(
      deg$direction == "Up"
    ),
    
    downregulated = sum(
      deg$direction == "Down"
    ),
    
    mapped_entrez = sum(
      !is.na(deg$ENTREZID) &
        deg$ENTREZID != ""
    ),
    
    design = paste(
      deparse(design_formula),
      collapse = ""
    )
  )
  
  cat("\nComparison summary:\n")
  print(summary_row)
  
  summary_row
}

# Run 1-day and 3-day comparisons

summary_results <- bind_rows(
  run_deseq_comparison("1d"),
  run_deseq_comparison("3d")
)

# Save overall summary

write_csv(
  summary_results,
  file.path(
    output_directory,
    "GSE238220_DE_summary.csv"
  )
)

cat("\nFinal differential-expression summary:\n")
print(summary_results)

cat("\nScript 04 complete.\n")
