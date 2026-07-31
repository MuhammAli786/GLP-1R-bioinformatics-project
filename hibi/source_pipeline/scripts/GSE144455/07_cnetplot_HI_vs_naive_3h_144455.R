############################################################
# GSE144455 - GO Cnet Plot
# Comparison: HI + PBS versus Naive + PBS at 3 hours
############################################################

library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)

#-----------------------------------------------------------
# 1. Load GO enrichment object
#-----------------------------------------------------------

ego_file <-
  "results/GSE144455_HI_vs_naive_3h_GO_BP.rds"

if (!file.exists(ego_file)) {
  stop(
    "GO enrichment object not found: ",
    ego_file,
    "\nRun script 06 first."
  )
}

ego <- readRDS(ego_file)

if (nrow(as.data.frame(ego)) == 0) {
  stop(
    "The GO enrichment object contains no significant terms, ",
    "so a cnet plot cannot be created."
  )
}

#-----------------------------------------------------------
# 2. Load annotated DEG results
#-----------------------------------------------------------

deg_file <-
  "results/GSE144455_HI_vs_naive_3h_DEGs_annotated.csv"

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

#-----------------------------------------------------------
# 3. Detect the gene-symbol column
#-----------------------------------------------------------

possible_symbol_columns <- c(
  "Gene.Symbol",
  "Gene symbol",
  "Gene.symbol",
  "GENE_SYMBOL",
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
  
  stop(
    "No gene-symbol column was detected."
  )
}

cat(
  "Using symbol column:",
  symbol_column,
  "\n"
)

#-----------------------------------------------------------
# 4. Clean symbols and create fold-change vector
#-----------------------------------------------------------

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

# When multiple probes map to the same symbol,
# retain the probe with the largest absolute logFC.
fold_change_data <- fold_change_data[
  order(
    abs(fold_change_data$logFC),
    decreasing = TRUE
  ),
]

fold_change_data <- fold_change_data[
  !duplicated(fold_change_data$clean_symbol),
]

geneList <- fold_change_data$logFC
names(geneList) <- fold_change_data$clean_symbol

geneList <- sort(
  geneList,
  decreasing = TRUE
)

cat(
  "Named fold-change values:",
  length(geneList),
  "\n"
)

#-----------------------------------------------------------
# 5. Check identifier overlap
#-----------------------------------------------------------

go_genes <- unique(
  unlist(
    strsplit(
      as.data.frame(ego)$geneID,
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
    "No overlap exists between geneList names and genes stored in ego. ",
    "Check whether ego is readable and whether the symbol column is correct."
  )
}

#-----------------------------------------------------------
# 6. Build cnet plot
#-----------------------------------------------------------

show_n <- min(
  15,
  nrow(as.data.frame(ego))
)

cnet_plot <- cnetplot(
  ego,
  showCategory = show_n,
  foldChange = geneList,
  node_label = "all"
)

#-----------------------------------------------------------
# 7. Save plot
#-----------------------------------------------------------

pdf(
  "figures/GSE144455_HI_vs_naive_3h_cnetplot.pdf",
  width = 16,
  height = 14
)

print(cnet_plot)

dev.off()

ggsave(
  filename =
    "figures/GSE144455_HI_vs_naive_3h_cnetplot.png",
  plot = cnet_plot,
  width = 13,
  height = 11,
  units = "in",
  dpi = 300
)

cat("\nCnet plot complete.\n")