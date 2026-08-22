# Create GO cnet plot for GSE144456 P5 HI vs control at 3h
# Input: results/GSE144456_P5_3h_GO_BP.rds, results/GSE144456_P5_3h_DEGs_annotated.csv -> Output: figures/GSE144456_P5_3h_cnetplot.pdf/.png

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

ego_file <- "results/GSE144456_P5_3h_GO_BP.rds"

if (!file.exists(ego_file)) {
  stop(
    "GO enrichment object not found: ",
    ego_file,
    "\nRun script 06 first."
  )
}

ego <- readRDS(ego_file)
go_results <- as.data.frame(ego)

cat("Significant GO BP terms:", nrow(go_results), "\n")

if (nrow(go_results) == 0) {
  stop(
    "No significant GO BP terms were found, ",
    "so a cnetplot cannot be created."
  )
}

deg_file <- "results/GSE144456_P5_3h_DEGs_annotated.csv"

if (!file.exists(deg_file)) {
  stop(
    "Annotated DEG file not found: ",
    deg_file,
    "\nRun script 04 first."
  )
}

deg <- read.csv(
  deg_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

possible_symbol_columns <- c(
  "GENE_SYMBOL",
  "Gene.Symbol",
  "Gene symbol",
  "Gene.symbol",
  "SYMBOL",
  "Symbol",
  "gene_symbol",
  "GeneSymbol"
)

symbol_column <- possible_symbol_columns[
  possible_symbol_columns %in% colnames(deg)
][1]

if (is.na(symbol_column)) {
  
  symbol_matches <- grep(
    "gene.*symbol|symbol",
    colnames(deg),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(symbol_matches) > 0) {
    symbol_column <- symbol_matches[1]
  }
}

if (is.na(symbol_column)) {
  
  cat("\nAvailable columns:\n")
  print(colnames(deg))
  
  stop("No gene-symbol column was detected.")
}

cat("Using symbol column:", symbol_column, "\n")

clean_symbols <- function(x) {
  
  x <- as.character(x)
  
  x <- sub("\\s*///.*$", "", x)
  x <- sub("\\s*//.*$", "", x)
  x <- sub("\\s*;.*$", "", x)
  x <- sub("\\s*\\|.*$", "", x)
  
  x <- trimws(x)
  
  x[
    is.na(x) |
      x == "" |
      x == "---" |
      x == "NA"
  ] <- NA_character_
  
  x
}

deg$clean_symbol <- clean_symbols(
  deg[[symbol_column]]
)

fold_change_data <- deg[
  !is.na(deg$clean_symbol) &
    !is.na(deg$logFC),
  c("clean_symbol", "logFC"),
  drop = FALSE
]

fold_change_data <- fold_change_data[
  order(
    abs(fold_change_data$logFC),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

fold_change_data <- fold_change_data[
  !duplicated(fold_change_data$clean_symbol),
  ,
  drop = FALSE
]

geneList <- fold_change_data$logFC
names(geneList) <- fold_change_data$clean_symbol

geneList <- sort(
  geneList,
  decreasing = TRUE
)

cat("Named fold-change values:", length(geneList), "\n")

go_genes <- unique(
  unlist(
    strsplit(
      go_results$geneID,
      "/",
      fixed = TRUE
    )
  )
)

overlap_count <- length(
  intersect(
    names(geneList),
    go_genes
  )
)

cat(
  "Gene-symbol overlap with GO object:",
  overlap_count,
  "\n"
)

if (overlap_count == 0) {
  stop(
    "No overlap exists between the fold-change names ",
    "and genes stored in the GO enrichment object."
  )
}

simplify_terms <- FALSE

if (simplify_terms) {
  
  ego_plot <- simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
} else {
  
  ego_plot <- ego
}

show_n <- min(
  10,
  nrow(as.data.frame(ego_plot))
)

cnet_plot <- cnetplot(
  ego_plot,
  showCategory = show_n,
  foldChange = geneList,
  node_label = "all"
)

pdf(
  "figures/GSE144456_P5_3h_cnetplot.pdf",
  width = 14,
  height = 12
)

print(cnet_plot)

dev.off()

ggsave(
  filename = "figures/GSE144456_P5_3h_cnetplot.png",
  plot = cnet_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300
)

cat(
  "\nCnet plot complete using ",
  show_n,
  " GO terms.\n",
  sep = ""
)