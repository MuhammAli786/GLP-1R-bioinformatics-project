library(EnhancedVolcano)

pdf(
  "figures/GSE23319_volcano.pdf",
  width = 8,
  height = 6
)

EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 0.2,
  title = "GSE23319\n3 h HI vs Sham"
)

dev.off()