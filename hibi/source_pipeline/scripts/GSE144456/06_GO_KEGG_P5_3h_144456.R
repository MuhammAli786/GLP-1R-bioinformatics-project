# GO and KEGG enrichment for GSE144456 P5 HI vs control at 3h
# Input: results/GSE144456_P5_3h_all_annotated.csv -> Output: results/GSE144456_P5_3h_GO_BP.csv/.rds, results/GSE144456_P5_3h_KEGG.csv/.rds, figures/GO_BP_*.pdf/.png, figures/KEGG_*.pdf/.png

library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)

results_file <- "results/GSE144456_P5_3h_all_annotated.csv"

if (!file.exists(results_file)) {
  stop(
    "Results file not found: ",
    results_file,
    "\nRun script 04 first."
  )
}

de_results <- read.csv(
  results_file,
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
  possible_symbol_columns %in% colnames(de_results)
][1]

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
  
  cat("\nAvailable columns:\n")
  print(colnames(de_results))
  
  stop(
    "No gene-symbol column was detected."
  )
}

cat(
  "\nUsing gene-symbol column:",
  symbol_column,
  "\n"
)

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

de_results$clean_symbol <- clean_symbols(
  de_results[[symbol_column]]
)

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
  "Unique significant symbols:",
  length(unique(deg$clean_symbol)),
  "\n"
)

if (nrow(deg) == 0) {
  stop(
    "No significant annotated DEGs are available for enrichment."
  )
}

deg_conversion <- bitr(
  unique(deg$clean_symbol),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

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

write.csv(
  deg_conversion,
  "results/GSE144456_P5_3h_DEG_symbol_to_Entrez.csv",
  row.names = FALSE
)

if (length(entrez_deg) < 10) {
  warning(
    "Fewer than 10 DEG Entrez IDs were mapped. ",
    "Enrichment may return no significant terms."
  )
}

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
  "results/GSE144456_P5_3h_GO_BP.csv",
  row.names = FALSE
)

saveRDS(
  ego,
  "results/GSE144456_P5_3h_GO_BP.rds"
)

if (nrow(go_results) > 0) {
  
  show_n <- min(
    15,
    nrow(go_results)
  )
  
  go_dotplot <- dotplot(
    ego,
    showCategory = show_n,
    title = "GSE144456 GO BP: P5 HI versus Control, 3 hours"
  )
  
  go_barplot <- barplot(
    ego,
    showCategory = show_n,
    title = "GSE144456 GO BP: P5 HI versus Control, 3 hours",
    font.size = 10
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_GO_BP_dotplot.pdf",
    plot = go_dotplot,
    width = 11,
    height = 8
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_GO_BP_dotplot.png",
    plot = go_dotplot,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_GO_BP_barplot.pdf",
    plot = go_barplot,
    width = 11,
    height = 8
  )
  
} else {
  
  message(
    "No significant GO BP terms were found. ",
    "No GO plots were created."
  )
}

ekegg <- enrichKEGG(
  gene = entrez_deg,
  universe = entrez_background,
  
  organism = "mmu",
  keyType = "ncbi-geneid",
  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

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
  "results/GSE144456_P5_3h_KEGG.csv",
  row.names = FALSE
)

saveRDS(
  ekegg,
  "results/GSE144456_P5_3h_KEGG.rds"
)

if (nrow(kegg_results) > 0) {
  
  show_n <- min(
    15,
    nrow(kegg_results)
  )
  
  kegg_dotplot <- dotplot(
    ekegg,
    showCategory = show_n,
    title = "GSE144456 KEGG: P5 HI versus Control, 3 hours"
  )
  
  kegg_barplot <- barplot(
    ekegg,
    showCategory = show_n,
    title = "GSE144456 KEGG: P5 HI versus Control, 3 hours",
    font.size = 10
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_KEGG_dotplot.pdf",
    plot = kegg_dotplot,
    width = 11,
    height = 8
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_KEGG_dotplot.png",
    plot = kegg_dotplot,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_KEGG_barplot.pdf",
    plot = kegg_barplot,
    width = 11,
    height = 8
  )
  
} else {
  
  message(
    "No significant KEGG pathways were found. ",
    "No KEGG plots were created."
  )
}

cat("\n============================================\n")
cat("GSE144456 enrichment analysis complete\n")
cat("============================================\n")

cat(
  "Significant annotated probes:",
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