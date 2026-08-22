# Create volcano plot for GSE23317 cortex 3h HI vs Sham
# Input: results/GSE23317_3h_all_genes.csv -> Output: figures/GSE23317_3h_volcano.pdf

library(EnhancedVolcano)

results <- read.csv(
  "results/GSE23317_3h_all_genes.csv",
  row.names = 1
)

pdf("figures/GSE23317_3h_volcano.pdf",
    width = 8,
    height = 6)

EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 0.2,
  title = "GSE23317 Cortex 3h",
  subtitle = "HI vs Sham"
)

dev.off()