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
  "results/Reactome_GSE23319_3h.csv",
  row.names = FALSE
)

#Dotplot
pdf(
  "figures/Reactome_dotplot.pdf",
  width = 10,
  height = 8
)

dotplot(reactome, showCategory = 20)

dev.off()

#Barplot
pdf(
  "figures/Reactome_barplot_GSE23319.pdf",
  width = 10,
  height = 8
)

barplot(reactome, showCategory = 20)

dev.off()
