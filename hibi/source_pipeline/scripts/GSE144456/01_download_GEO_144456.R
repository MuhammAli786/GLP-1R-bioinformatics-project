############################################################
# GSE144456 - Download GEO data
############################################################

library(GEOquery)

gse <- getGEO(
  "GSE144456",
  GSEMatrix = TRUE
)

cat("Number of ExpressionSet objects:", length(gse), "\n")

for (i in seq_along(gse)) {
  cat(
    "ExpressionSet", i,
    "- platform:", annotation(gse[[i]]),
    "- dimensions:", paste(dim(exprs(gse[[i]])), collapse = " x "),
    "\n"
  )
}

saveRDS(
  gse,
  "data/GSE144456.rds"
)

cat("GSE144456 download complete.\n")