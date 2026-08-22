# GSE236133 differential expression using normalized counts; input: data/GSE236133_prepared_data.rds -> output: results/GSE236133/DE/
# Input gene identifiers are mouse SYMBOLS. Contrast: positive logFC = higher in ipsilateral, negative = higher in contralateral

library(limma)
library(dplyr)
library(readr)
library(tibble)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Paths and thresholds: logFC >= 0.2, FDR < 0.05
input_file <- "data/GSE236133_prepared_data.rds"

results_directory <- "results/GSE236133/DE"
metadata_directory <- "metadata"

logfc_threshold <- 0.2
fdr_threshold <- 0.05

dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  metadata_directory,
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

# Load prepared data
prepared_data <- readRDS(input_file)

expression_matrix <- prepared_data$expression
metadata <- as.data.frame(prepared_data$metadata)

if (!identical(
  colnames(expression_matrix),
  as.character(metadata$count_column)
)) {
  stop(
    "Expression columns and metadata rows are not aligned."
  )
}

cat("\nOriginal expression dimensions:\n")
print(dim(expression_matrix))

# Clean gene symbols
gene_symbols <- rownames(expression_matrix)

if (is.null(gene_symbols)) {
  stop("Expression matrix has no gene identifiers.")
}

gene_symbols <- trimws(as.character(gene_symbols))

valid_gene <- !is.na(gene_symbols) &
  gene_symbols != "" &
  gene_symbols != "NA"

expression_matrix <- expression_matrix[
  valid_gene,
  ,
  drop = FALSE
]

gene_symbols <- gene_symbols[valid_gene]

rownames(expression_matrix) <- gene_symbols

cat(
  "\nRemoved",
  sum(!valid_gene),
  "rows with missing or blank gene symbols.\n"
)

# Collapse duplicate gene symbols by mean expression
if (anyDuplicated(rownames(expression_matrix))) {
  
  cat(
    "\nDuplicate gene symbols detected.",
    "\nCollapsing duplicate rows by mean expression.\n"
  )
  
  expression_matrix <- limma::avereps(
    expression_matrix,
    ID = rownames(expression_matrix)
  )
}

cat("\nExpression dimensions after cleaning:\n")
print(dim(expression_matrix))

# Validate and transform expression values
storage.mode(expression_matrix) <- "numeric"

if (anyNA(expression_matrix)) {
  warning(
    "Missing expression values detected: ",
    sum(is.na(expression_matrix))
  )
}

if (any(expression_matrix < 0, na.rm = TRUE)) {
  stop(
    "Negative values were found. ",
    "Do not apply log2(x + 1) until the data are inspected."
  )
}

cat("\nNormalized-count summary before transformation:\n")
print(summary(as.vector(expression_matrix)))

# These are normalized counts, so use log2(x + 1)
expression_log2 <- log2(expression_matrix + 1)

cat("\nLog2-expression summary:\n")
print(summary(as.vector(expression_log2)))

saveRDS(
  expression_log2,
  "data/GSE236133_log2_expression.rds"
)

# Filter very low-expression genes: keep genes with normalized expression > 1 in at least two samples
keep_gene <- rowSums(
  expression_matrix > 1,
  na.rm = TRUE
) >= 2

expression_filtered <- expression_log2[
  keep_gene,
  ,
  drop = FALSE
]

cat(
  "\nGenes before filtering:",
  nrow(expression_log2),
  "\n"
)

cat(
  "Genes after filtering:",
  nrow(expression_filtered),
  "\n"
)

if (nrow(expression_filtered) == 0) {
  stop("No genes remained after expression filtering.")
}

# Map gene SYMBOLS to Entrez and Ensembl IDs
symbols <- rownames(expression_filtered)

gene_annotation <- tibble(
  SYMBOL = symbols
) %>%
  mutate(
    ENTREZID = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = SYMBOL,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    ),
    
    ENSEMBL = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = SYMBOL,
      keytype = "SYMBOL",
      column = "ENSEMBL",
      multiVals = "first"
    ),
    
    GENENAME = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = SYMBOL,
      keytype = "SYMBOL",
      column = "GENENAME",
      multiVals = "first"
    )
  )

write_csv(
  gene_annotation,
  file.path(
    metadata_directory,
    "GSE236133_gene_annotation.csv"
  )
)

cat(
  "\nTotal unique symbols:",
  nrow(gene_annotation),
  "\n"
)

cat(
  "Mapped Entrez IDs:",
  sum(!is.na(gene_annotation$ENTREZID)),
  "\n"
)

cat(
  "Mapped Ensembl IDs:",
  sum(!is.na(gene_annotation$ENSEMBL)),
  "\n"
)

# Define comparisons: WT, NEIL1-KO, NEIL2-KO at 3h and 6h
comparisons <- tribble(
  ~genotype,  ~time, ~comparison_name,
  
  "WT",
  "3h",
  "WT_3h_ipsilateral_vs_contralateral",
  
  "WT",
  "6h",
  "WT_6h_ipsilateral_vs_contralateral",
  
  "NEIL1-KO",
  "3h",
  "NEIL1KO_3h_ipsilateral_vs_contralateral",
  
  "NEIL1-KO",
  "6h",
  "NEIL1KO_6h_ipsilateral_vs_contralateral",
  
  "NEIL2-KO",
  "3h",
  "NEIL2KO_3h_ipsilateral_vs_contralateral",
  
  "NEIL2-KO",
  "6h",
  "NEIL2KO_6h_ipsilateral_vs_contralateral"
)

# Function to run one limma comparison
run_limma_comparison <- function(
    selected_genotype,
    selected_time,
    comparison_name
) {
  
  cat("\n=================================================\n")
  cat("Running:", comparison_name, "\n")
  cat("=================================================\n")
  
  sample_metadata <- metadata %>%
    filter(
      genotype == selected_genotype,
      time == selected_time
    )
  
  cat("\nSelected samples:\n")
  
  print(
    sample_metadata[
      ,
      c(
        "count_column",
        "genotype",
        "region",
        "biological_replicate",
        "time"
      ),
      drop = FALSE
    ]
  )
  
  if (nrow(sample_metadata) != 4) {
    stop(
      comparison_name,
      ": expected 4 samples, but found ",
      nrow(sample_metadata),
      "."
    )
  }
  
  region_counts <- table(sample_metadata$region)
  
  cat("\nRegion counts:\n")
  print(region_counts)
  
  if (
    !all(
      c("contralateral", "ipsilateral") %in%
      names(region_counts)
    )
  ) {
    stop(
      comparison_name,
      ": both regions are not represented."
    )
  }
  
  if (
    region_counts["contralateral"] != 2 ||
    region_counts["ipsilateral"] != 2
  ) {
    stop(
      comparison_name,
      ": expected two samples in each region."
    )
  }
  
  selected_expression <- expression_filtered[
    ,
    sample_metadata$count_column,
    drop = FALSE
  ]
  
  # Contralateral is the reference level
  group <- factor(
    sample_metadata$region,
    levels = c(
      "contralateral",
      "ipsilateral"
    )
  )
  
  design <- model.matrix(
    ~ 0 + group
  )
  
  colnames(design) <- levels(group)
  
  contrast_matrix <- makeContrasts(
    ipsilateral - contralateral,
    levels = design
  )
  
  fit <- lmFit(
    selected_expression,
    design
  )
  
  fit <- contrasts.fit(
    fit,
    contrast_matrix
  )
  
  fit <- eBayes(
    fit,
    trend = TRUE,
    robust = TRUE
  )
  
  result <- topTable(
    fit,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  ) %>%
    rownames_to_column("SYMBOL") %>%
    left_join(
      gene_annotation,
      by = "SYMBOL"
    ) %>%
    mutate(
      comparison = comparison_name,
      
      direction = case_when(
        adj.P.Val < fdr_threshold &
          logFC >= logfc_threshold ~ "Up",
        
        adj.P.Val < fdr_threshold &
          logFC <= -logfc_threshold ~ "Down",
        
        TRUE ~ "Not_significant"
      ),
      
      significant = direction != "Not_significant"
    ) %>%
    select(
      comparison,
      SYMBOL,
      ENTREZID,
      ENSEMBL,
      GENENAME,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      direction,
      significant
    )
  
  deg <- result %>%
    filter(significant)
  
  # Save all genes and significant DEGs
  write_csv(
    result,
    file.path(
      results_directory,
      paste0(
        comparison_name,
        "_all_genes.csv"
      )
    )
  )
  
  write_csv(
    deg,
    file.path(
      results_directory,
      paste0(
        comparison_name,
        "_DEGs.csv"
      )
    )
  )
  
  # Save gene lists
  symbol_list <- deg %>%
    filter(
      !is.na(SYMBOL),
      SYMBOL != ""
    ) %>%
    distinct(SYMBOL) %>%
    pull(SYMBOL)
  
  entrez_list <- deg %>%
    filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID) %>%
    as.character()
  
  write_lines(
    symbol_list,
    file.path(
      results_directory,
      paste0(
        comparison_name,
        "_gene_symbols.txt"
      )
    )
  )
  
  write_lines(
    entrez_list,
    file.path(
      results_directory,
      paste0(
        comparison_name,
        "_entrez_ids.txt"
      )
    )
  )
  
  # Summary
  summary_row <- tibble(
    comparison = comparison_name,
    genotype = selected_genotype,
    time = selected_time,
    total_tested = nrow(result),
    significant_DEGs = nrow(deg),
    upregulated = sum(deg$direction == "Up"),
    downregulated = sum(deg$direction == "Down"),
    mapped_entrez_DEGs = sum(
      !is.na(deg$ENTREZID) &
        deg$ENTREZID != ""
    )
  )
  
  cat("\nDEG summary:\n")
  print(summary_row)
  
  return(summary_row)
}

# Run all comparisons
summary_results <- bind_rows(
  lapply(
    seq_len(nrow(comparisons)),
    function(i) {
      
      run_limma_comparison(
        selected_genotype = comparisons$genotype[i],
        selected_time = comparisons$time[i],
        comparison_name = comparisons$comparison_name[i]
      )
    }
  )
)

# Save summary
write_csv(
  summary_results,
  file.path(
    results_directory,
    "GSE236133_DE_summary.csv"
  )
)

cat("\nFinal differential-expression summary:\n")
print(summary_results)

cat("\nScript 04 complete.\n")