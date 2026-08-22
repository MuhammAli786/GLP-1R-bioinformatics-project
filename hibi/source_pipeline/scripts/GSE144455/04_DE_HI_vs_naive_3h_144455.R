# Differential expression for GSE144455: HI+PBS vs Naive+PBS at 3h using intercept-only limma
# Input: data/GSE144455.rds, metadata/GSE144455_HI_vs_naive_PBS_3h_samples.rds -> Output: results/GSE144455_HI_vs_naive_3h_*.csv/.rds

library(GEOquery)
library(limma)

if (file.exists("data/GSE144455.rds")) {
  
  gse <- readRDS("data/GSE144455.rds")
  
} else if (!exists("gse")) {
  
  gse <- getGEO("GSE144455", GSEMatrix = TRUE)
  
}

eset <- gse[[1]]
expr <- exprs(eset)
pheno <- pData(eset)

selected_metadata <- readRDS(
  "metadata/GSE144455_HI_vs_naive_PBS_3h_samples.rds"
)

selected_samples <- selected_metadata$sample_id

if (!all(selected_samples %in% colnames(expr))) {
  stop("One or more selected samples are missing from the expression matrix.")
}

expr_selected <- expr[, selected_samples, drop = FALSE]

cat("\nSelected samples:\n")
print(colnames(expr_selected))

cat("\nExpression matrix dimensions:\n")
print(dim(expr_selected))

cat("\nExpression-value summary:\n")
print(summary(as.vector(expr_selected)))

keep <- rowSums(expr_selected != 0, na.rm = TRUE) > 0

expr_filtered <- expr_selected[keep, , drop = FALSE]

cat("\nProbes before filtering:", nrow(expr_selected), "\n")
cat("Probes after filtering:", nrow(expr_filtered), "\n")

design <- matrix(
  1,
  nrow = ncol(expr_filtered),
  ncol = 1
)

rownames(design) <- colnames(expr_filtered)
colnames(design) <- "HI_vs_Naive_3h"

print(design)

fit <- lmFit(expr_filtered, design)
fit <- eBayes(fit)

de_all <- topTable(
  fit,
  coef = "HI_vs_Naive_3h",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

de_all$probe_id <- rownames(de_all)

de_all <- de_all[, c(
  "probe_id",
  setdiff(colnames(de_all), "probe_id")
)]

feature_annotation <- fData(eset)
feature_annotation$probe_id <- rownames(feature_annotation)

de_annotated <- merge(
  de_all,
  feature_annotation,
  by = "probe_id",
  all.x = TRUE,
  sort = FALSE
)

de_annotated <- de_annotated[
  match(de_all$probe_id, de_annotated$probe_id),
]

deg <- de_annotated[
  !is.na(de_annotated$logFC) &
    !is.na(de_annotated$adj.P.Val) &
    abs(de_annotated$logFC) >= 0.2 &
    de_annotated$adj.P.Val < 0.05,
]

deg$direction <- ifelse(
  deg$logFC > 0,
  "Up",
  "Down"
)

cat("\nTotal significant DEGs:", nrow(deg), "\n")
cat("\nDirection counts:\n")
print(table(deg$direction))

write.csv(
  de_all,
  "results/GSE144455_HI_vs_naive_3h_all_limma_results.csv",
  row.names = FALSE
)

write.csv(
  de_annotated,
  "results/GSE144455_HI_vs_naive_3h_all_annotated.csv",
  row.names = FALSE
)

write.csv(
  deg,
  "results/GSE144455_HI_vs_naive_3h_DEGs_annotated.csv",
  row.names = FALSE
)

saveRDS(
  fit,
  "results/GSE144455_HI_vs_naive_3h_limma_fit.rds"
)

saveRDS(
  deg,
  "results/GSE144455_HI_vs_naive_3h_DEGs_annotated.rds"
)

cat("\nDifferential-expression analysis complete.\n")