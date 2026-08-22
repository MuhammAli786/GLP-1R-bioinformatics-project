# Reactome enrichment for GSE144456 P5 HI vs control at 3h
# Input: results/GSE144456_P5_3h_DEG_symbol_to_Entrez.csv -> Output: results/GSE144456_P5_3h_Reactome.csv/.rds, figures/Reactome_*.pdf/.png

library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)

conversion_file <-
  "results/GSE144456_P5_3h_DEG_symbol_to_Entrez.csv"

if (!file.exists(conversion_file)) {
  stop(
    "Entrez conversion file not found: ",
    conversion_file,
    "\nRun script 06 first."
  )
}

deg_conversion <- read.csv(
  conversion_file,
  stringsAsFactors = FALSE
)

if (!"ENTREZID" %in% colnames(deg_conversion)) {
  stop("The conversion file does not contain an ENTREZID column.")
}

entrez_deg <- unique(
  na.omit(as.character(deg_conversion$ENTREZID))
)

cat("Mapped DEG Entrez IDs:", length(entrez_deg), "\n")

if (length(entrez_deg) == 0) {
  stop("No Entrez IDs are available for Reactome enrichment.")
}

if (length(entrez_deg) < 10) {
  warning(
    "Fewer than 10 Entrez IDs are available. ",
    "Reactome may return no significant pathways."
  )
}

reactome <- enrichPathway(
  gene = entrez_deg,
  organism = "mouse",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

reactome_results <- as.data.frame(reactome)

cat(
  "Significant Reactome pathways:",
  nrow(reactome_results),
  "\n"
)

write.csv(
  reactome_results,
  "results/GSE144456_P5_3h_Reactome.csv",
  row.names = FALSE
)

saveRDS(
  reactome,
  "results/GSE144456_P5_3h_Reactome.rds"
)

if (nrow(reactome_results) > 0) {
  
  show_n <- min(
    15,
    nrow(reactome_results)
  )
  
  reactome_dotplot <- dotplot(
    reactome,
    showCategory = show_n,
    title = "GSE144456 Reactome: P5 HI versus Control, 3 hours"
  )
  
  reactome_barplot <- barplot(
    reactome,
    showCategory = show_n,
    title = "GSE144456 Reactome: P5 HI versus Control, 3 hours",
    font.size = 10
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_Reactome_dotplot.pdf",
    plot = reactome_dotplot,
    width = 11,
    height = 8
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_Reactome_dotplot.png",
    plot = reactome_dotplot,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    "figures/GSE144456_P5_3h_Reactome_barplot.pdf",
    plot = reactome_barplot,
    width = 11,
    height = 8
  )
  
} else {
  
  message(
    "No significant Reactome pathways were found. ",
    "No Reactome plots were created."
  )
}

cat("\nReactome enrichment complete.\n")