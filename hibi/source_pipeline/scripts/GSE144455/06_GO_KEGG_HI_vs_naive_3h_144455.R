############################################################
# GSE144455 - GO and KEGG Enrichment
# Comparison: HI + PBS versus Naive + PBS at 3 hours
############################################################

library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)

#-----------------------------------------------------------
# 1. Load annotated limma results
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
# 2. Locate the gene-symbol column
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

# Try a broader search if no exact match was found
if (is.na(symbol_column)) {
  
  symbol_matches <- grep(
    "gene.*symbol|symbol",
    colnames(de_results),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(symbol_matches) > 0) {
    symbol_column <- symbol_matches[1]
  }
}

if (is.na(symbol_column)) {
  
  cat("\nAvailable annotation columns:\n")
  print(colnames(de_results))
  
  stop(
    "No gene-symbol column was detected. ",
    "Inspect the printed column names."
  )
}

cat(
  "\nUsing gene-symbol column:",
  symbol_column,
  "\n"
)

#-----------------------------------------------------------
# 3. Clean gene symbols
#-----------------------------------------------------------

clean_symbols <- function(x) {
  
  x <- as.character(x)
  
  # Keep the first annotation when several are present
  x <- sub("\\s*///.*$", "", x)
  x <- sub("\\s*//.*$", "", x)
  x <- sub("\\s*;.*$", "", x)
  x <- sub("\\s*\\|.*$", "", x)
  
  x <- trimws(x)
  
  x[
    x == "" |
      x == "---" |
      x == "NA"
  ] <- NA_character_
  
  x
}

de_results$clean_symbol <- clean_symbols(
  de_results[[symbol_column]]
)

#-----------------------------------------------------------
# 4. Define significant DEGs
#-----------------------------------------------------------

deg <- de_results[
  !is.na(de_results$clean_symbol) &
    !is.na(de_results$logFC) &
    !is.na(de_results$adj.P.Val) &
    abs(de_results$logFC) >= 0.2 &
    de_results$adj.P.Val < 0.05,
  ,
  drop = FALSE
]

cat(
  "\nSignificant annotated probes:",
  nrow(deg),
  "\n"
)

cat(
  "Unique significant gene symbols:",
  length(unique(deg$clean_symbol)),
  "\n"
)

#-----------------------------------------------------------
# 5. Convert gene symbols to Entrez IDs
#-----------------------------------------------------------

deg_conversion <- bitr(
  unique(deg$clean_symbol),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# Create background from all tested and annotated genes
background_symbols <- unique(
  na.omit(de_results$clean_symbol)
)

background_conversion <- bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

entrez_deg <- unique(
  deg_conversion$ENTREZID
)

entrez_background <- unique(
  background_conversion$ENTREZID
)

cat(
  "Mapped DEG Entrez IDs:",
  length(entrez_deg),
  "\n"
)

cat(
  "Mapped background Entrez IDs:",
  length(entrez_background),
  "\n"
)

if (length(entrez_deg) < 10) {
  stop(
    "Fewer than 10 DEG Entrez IDs were mapped. ",
    "Check the annotation column and gene symbols."
  )
}

# Save the symbol-to-Entrez conversion
write.csv(
  deg_conversion,
  paste0(
    "results/",
    "GSE144455_HI_vs_naive_3h_",
    "DEG_symbol_to_Entrez.csv"
  ),
  row.names = FALSE
)

#-----------------------------------------------------------
# 6. GO Biological Process enrichment
#-----------------------------------------------------------

ego <- enrichGO(
  gene = entrez_deg,
  universe = entrez_background,
  
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  
  ont = "BP",
  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  
  readable = TRUE
)

go_results <- as.data.frame(ego)

cat(
  "\nSignificant GO BP terms:",
  nrow(go_results),
  "\n"
)

write.csv(
  go_results,
  "results/GSE144455_HI_vs_naive_3h_GO_BP.csv",
  row.names = FALSE
)

saveRDS(
  ego,
  "results/GSE144455_HI_vs_naive_3h_GO_BP.rds"
)

#-----------------------------------------------------------
# 7. Save GO plots only when results exist
#-----------------------------------------------------------

if (nrow(go_results) > 0) {
  
  number_to_show <- min(
    15,
    nrow(go_results)
  )
  
  go_dotplot <- dotplot(
    ego,
    showCategory = number_to_show,
    title = paste(
      "GSE144455 GO Biological Process",
      "HI versus Naive, 3 hours"
    )
  )
  
  go_barplot <- barplot(
    ego,
    showCategory = number_to_show,
    title = paste(
      "GSE144455 GO Biological Process",
      "HI versus Naive, 3 hours"
    ),
    font.size = 10
  )
  
  ggsave(
    "figures/GSE144455_HI_vs_naive_3h_GO_BP_dotplot.pdf",
    plot = go_dotplot,
    width = 11,
    height = 8
  )
  
  ggsave(
    "figures/GSE144455_HI_vs_naive_3h_GO_BP_dotplot.png",
    plot = go_dotplot,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    "figures/GSE144455_HI_vs_naive_3h_GO_BP_barplot.pdf",
    plot = go_barplot,
    width = 11,
    height = 8
  )
  
} else {
  
  message(
    "No significant GO BP terms were found. ",
    "No GO plot was created."
  )
}

#-----------------------------------------------------------
# 8. KEGG pathway enrichment
#-----------------------------------------------------------

ekegg <- enrichKEGG(
  gene = entrez_deg,
  universe = entrez_background,
  
  organism = "mmu",
  keyType = "ncbi-geneid",
  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Convert Entrez IDs to readable mouse symbols
if (nrow(as.data.frame(ekegg)) > 0) {
  
  ekegg <- setReadable(
    ekegg,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID"
  )
}

kegg_results <- as.data.frame(ekegg)

cat(
  "\nSignificant KEGG pathways:",
  nrow(kegg_results),
  "\n"
)

write.csv(
  kegg_results,
  "results/GSE144455_HI_vs_naive_3h_KEGG.csv",
  row.names = FALSE
)

saveRDS(
  ekegg,
  "results/GSE144455_HI_vs_naive_3h_KEGG.rds"
)

#-----------------------------------------------------------
# 9. Save KEGG plots only when results exist
#-----------------------------------------------------------

if (nrow(kegg_results) > 0) {
  
  number_to_show <- min(
    15,
    nrow(kegg_results)
  )
  
  kegg_dotplot <- dotplot(
    ekegg,
    showCategory = number_to_show,
    title = paste(
      "GSE144455 KEGG pathways",
      "HI versus Naive, 3 hours"
    )
  )
  
  kegg_barplot <- barplot(
    ekegg,
    showCategory = number_to_show,
    title = paste(
      "GSE144455 KEGG pathways",
      "HI versus Naive, 3 hours"
    ),
    font.size = 10
  )
  
  ggsave(
    "figures/GSE144455_HI_vs_naive_3h_KEGG_dotplot.pdf",
    plot = kegg_dotplot,
    width = 11,
    height = 8
  )
  
  ggsave(
    "figures/GSE144455_HI_vs_naive_3h_KEGG_dotplot.png",
    plot = kegg_dotplot,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    "figures/GSE144455_HI_vs_naive_3h_KEGG_barplot.pdf",
    plot = kegg_barplot,
    width = 11,
    height = 8
  )
  
} else {
  
  message(
    "No significant KEGG pathways were found. ",
    "No KEGG plot was created."
  )
}

#-----------------------------------------------------------
# 10. Final summary
#-----------------------------------------------------------

cat("\n============================================\n")
cat("GSE144455 enrichment analysis complete\n")
cat("============================================\n")

cat(
  "Significant DEGs:",
  nrow(deg),
  "\n"
)

cat(
  "Mapped DEG Entrez IDs:",
  length(entrez_deg),
  "\n"
)

cat(
  "Significant GO BP terms:",
  nrow(go_results),
  "\n"
)

cat(
  "Significant KEGG pathways:",
  nrow(kegg_results),
  "\n"
)