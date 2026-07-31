############################################################
# GSE144455 - Volcano Plot
# Comparison: HI + PBS versus Naive + PBS at 3 hours
############################################################

library(EnhancedVolcano)
library(ggplot2)

#-----------------------------------------------------------
# 1. Load all annotated limma results
#-----------------------------------------------------------

results_file <-
  "results/GSE144455_HI_vs_naive_3h_all_annotated.csv"

if (!file.exists(results_file)) {
  stop(
    "Results file not found: ", results_file,
    "\nRun script 04 first."
  )
}

de_results <- read.csv(
  results_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#-----------------------------------------------------------
# 2. Check required columns
#-----------------------------------------------------------

required_columns <- c(
  "probe_id",
  "logFC",
  "adj.P.Val"
)

missing_columns <- setdiff(
  required_columns,
  colnames(de_results)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Replace missing adjusted p-values with 1 for plotting
de_results$adj.P.Val[
  is.na(de_results$adj.P.Val)
] <- 1

# Prevent zero adjusted p-values from producing Inf
smallest_nonzero <- min(
  de_results$adj.P.Val[
    de_results$adj.P.Val > 0
  ],
  na.rm = TRUE
)

de_results$adj.P.Val[
  de_results$adj.P.Val == 0
] <- smallest_nonzero

#-----------------------------------------------------------
# 3. Find an annotation column for plot labels
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
  possible_symbol_columns %in% colnames(de_results)
][1]

if (!is.na(symbol_column)) {
  
  plot_labels <- as.character(
    de_results[[symbol_column]]
  )
  
  # For entries containing multiple annotations,
  # retain only the first symbol.
  plot_labels <- sub(
    "\\s*///.*$",
    "",
    plot_labels
  )
  
  plot_labels <- sub(
    "\\s*//.*$",
    "",
    plot_labels
  )
  
  plot_labels <- sub(
    "\\s*;.*$",
    "",
    plot_labels
  )
  
  missing_labels <-
    is.na(plot_labels) |
    plot_labels == "" |
    plot_labels == "---"
  
  plot_labels[missing_labels] <-
    as.character(
      de_results$probe_id[missing_labels]
    )
  
} else {
  
  message(
    "No gene-symbol column was detected. ",
    "Probe IDs will be used as labels."
  )
  
  plot_labels <- as.character(
    de_results$probe_id
  )
}

de_results$plot_label <- plot_labels

#-----------------------------------------------------------
# 4. Identify genes to label
#-----------------------------------------------------------

# Label the 10 most statistically significant DEGs
significant_results <- de_results[
  abs(de_results$logFC) >= 0.2 &
    de_results$adj.P.Val < 0.05,
  ,
  drop = FALSE
]

significant_results <- significant_results[
  order(significant_results$adj.P.Val),
  ,
  drop = FALSE
]

genes_to_label <- head(
  significant_results$plot_label,
  10
)

#-----------------------------------------------------------
# 5. Create volcano plot
#-----------------------------------------------------------

volcano_plot <- EnhancedVolcano(
  de_results,
  
  lab = de_results$plot_label,
  x = "logFC",
  y = "adj.P.Val",
  
  selectLab = genes_to_label,
  
  title = "GSE144455: HI versus Naive at 3 hours",
  subtitle = "Mouse forebrain, PBS-treated groups",
  
  caption = paste0(
    "|logFC| ≥ 0.2; adjusted P < 0.05; ",
    nrow(significant_results),
    " DEGs"
  ),
  
  pCutoff = 0.05,
  FCcutoff = 0.2,
  
  pointSize = 2.0,
  labSize = 3.5,
  
  xlab = expression(Log[2] ~ "fold change"),
  ylab = expression(-Log[10] ~ "adjusted P-value"),
  
  legendPosition = "right",
  
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  
  max.overlaps = 20
)

#-----------------------------------------------------------
# 6. Save the volcano plot
#-----------------------------------------------------------

pdf(
  "figures/GSE144455_HI_vs_naive_3h_volcano.pdf",
  width = 10,
  height = 8
)

print(volcano_plot)

dev.off()

ggsave(
  filename =
    "figures/GSE144455_HI_vs_naive_3h_volcano.png",
  plot = volcano_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

cat("\nVolcano plot complete.\n")
cat(
  "Significant DEGs:",
  nrow(significant_results),
  "\n"
)