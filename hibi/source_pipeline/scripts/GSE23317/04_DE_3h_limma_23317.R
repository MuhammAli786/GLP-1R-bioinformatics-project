# Differential expression for GSE23317 cortex 3h HI vs Sham using limma
# Input: data/GSE23317.rds -> Output: results/GSE23317_3h_all_genes.csv, results/GSE23317_3h_DEGs.csv

library(GEOquery)
library(limma)

gse <- readRDS("data/GSE23317.rds")
eset <- gse[[1]]

expr <- exprs(eset)

pheno <- pData(eset)

keep <- grepl("3h", pheno$title) &
  (grepl("Sham", pheno$title) | grepl("HI-", pheno$title))

expr_sub <- expr[, keep]
pheno_sub <- pheno[keep, ]

pheno_sub$title

group <- ifelse(grepl("Sham", pheno_sub$title), "Sham", "HI")
group <- factor(group, levels = c("Sham", "HI"))

design <- model.matrix(~ group)

fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)

results <- topTable(
  fit,
  coef = "groupHI",
  number = Inf,
  adjust.method = "BH"
)

write.csv(
  results,
  "results/GSE23317_3h_all_genes.csv"
)

deg <- subset(
  results,
  abs(logFC) >= 0.2 &
    adj.P.Val < 0.05
)

write.csv(
  deg,
  "results/GSE23317_3h_DEGs.csv"
)

cat("Number of DEGs:", nrow(deg), "\n")