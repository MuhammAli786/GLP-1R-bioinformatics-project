library(GEOquery)
library(limma)

# Load dataset
gse <- readRDS("data/GSE23317.rds")
eset <- gse[[1]]

# Expression matrix
expr <- exprs(eset)

# Sample metadata
pheno <- pData(eset)

# Keep only 3h Sham and HI samples
keep <- grepl("3h", pheno$title) &
  (grepl("Sham", pheno$title) | grepl("HI-", pheno$title))

expr_sub <- expr[, keep]
pheno_sub <- pheno[keep, ]

# Check selected samples
pheno_sub$title

# Create group variable
group <- ifelse(grepl("Sham", pheno_sub$title), "Sham", "HI")
group <- factor(group, levels = c("Sham", "HI"))

# Design matrix
design <- model.matrix(~ group)

# Fit linear model
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)

# Differential expression
results <- topTable(
  fit,
  coef = "groupHI",
  number = Inf,
  adjust.method = "BH"
)

# Save all genes
write.csv(
  results,
  "results/GSE23317_3h_all_genes.csv"
)

# Filter DEGs
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