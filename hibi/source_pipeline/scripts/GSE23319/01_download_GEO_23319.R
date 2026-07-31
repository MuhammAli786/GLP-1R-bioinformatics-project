############################################################
## Download GSE23319
############################################################

library(GEOquery)

dir.create("data", showWarnings = FALSE)

gse <- getGEO(
  "GSE23319",
  GSEMatrix = TRUE
)

saveRDS(
  gse,
  file = "data/GSE23319.rds"
)

cat("GSE23319 downloaded successfully!\n")