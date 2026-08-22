# Reactome enrichment for GSE23317 3h HI vs Sham
# Input: (requires `entrez` variable from script 06) -> Output: results/Reactome_GSE23317_3h.csv, figures/Reactome_dotplot.pdf, figures/Reactome_barplot_GSE23317.pdf

library(ReactomePA)
library(enrichplot)

reactome <- enrichPathway(
  gene = entrez,
  organism = "mouse",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

head(as.data.frame(reactome))

reactome_results <- as.data.frame(reactome)

write.csv(
  reactome_results,
  "results/Reactome_GSE23317_3h.csv",
  row.names = FALSE
)

pdf(
  "figures/Reactome_dotplot.pdf",
  width = 10,
  height = 8
)

dotplot(reactome, showCategory = 20)

dev.off()

pdf(
  "figures/Reactome_barplot_GSE23317.pdf",
  width = 10,
  height = 8
)

barplot(reactome, showCategory = 20)

dev.off()
