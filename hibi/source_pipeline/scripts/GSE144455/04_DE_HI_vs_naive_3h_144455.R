############################################################
# GSE144455 - Differential Expression
# Comparison: HI + PBS versus Naive + PBS at 3 hours
#
# Two-colour Agilent array:
# Each selected array already contains a processed log-ratio.
############################################################

library(GEOquery)
library(limma)

#-----------------------------------------------------------
# 1. Load GEO data
#-----------------------------------------------------------

if (file.exists("data/GSE144455.rds")) {
  
  gse <- readRDS("data/GSE144455.rds")
  
} else if (!exists("gse")) {
  
  gse <- getGEO("GSE144455", GSEMatrix = TRUE)
  
}

eset <- gse[[1]]
expr <- exprs(eset)
pheno <- pData(eset)

#-----------------------------------------------------------
# 2. Load the samples selected in script 03
#-----------------------------------------------------------

selected_metadata <- readRDS(
  "metadata/GSE144455_HI_vs_naive_PBS_3h_samples.rds"
)

selected_samples <- selected_metadata$sample_id

# Confirm that all selected samples exist in the matrix
if (!all(selected_samples %in% colnames(expr))) {
  stop("One or more selected samples are missing from the expression matrix.")
}

# Extract the three HI-versus-naive arrays
expr_selected <- expr[, selected_samples, drop = FALSE]

cat("\nSelected samples:\n")
print(colnames(expr_selected))

cat("\nExpression matrix dimensions:\n")
print(dim(expr_selected))

cat("\nExpression-value summary:\n")
print(summary(as.vector(expr_selected)))

#-----------------------------------------------------------
# 3. Remove probes with no variation/information
#-----------------------------------------------------------

# Remove probes that are zero in all selected arrays
keep <- rowSums(expr_selected != 0, na.rm = TRUE) > 0

expr_filtered <- expr_selected[keep, , drop = FALSE]

cat("\nProbes before filtering:", nrow(expr_selected), "\n")
cat("Probes after filtering:", nrow(expr_filtered), "\n")

#-----------------------------------------------------------
# 4. Intercept-only limma model
#-----------------------------------------------------------

# Each column is already an HI-versus-naive log-ratio.
# The intercept estimates the average log-ratio across arrays.
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

#-----------------------------------------------------------
# 5. Extract all differential-expression results
#-----------------------------------------------------------

de_all <- topTable(
  fit,
  coef = "HI_vs_Naive_3h",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

de_all$probe_id <- rownames(de_all)

# Put probe ID first
de_all <- de_all[, c(
  "probe_id",
  setdiff(colnames(de_all), "probe_id")
)]

#-----------------------------------------------------------
# 6. Add GPL feature annotation
#-----------------------------------------------------------

feature_annotation <- fData(eset)
feature_annotation$probe_id <- rownames(feature_annotation)

de_annotated <- merge(
  de_all,
  feature_annotation,
  by = "probe_id",
  all.x = TRUE,
  sort = FALSE
)

# Restore the limma result ordering after merge
de_annotated <- de_annotated[
  match(de_all$probe_id, de_annotated$probe_id),
]

#-----------------------------------------------------------
# 7. Apply DEG thresholds
#-----------------------------------------------------------

# Same threshold used for GSE23317 and GSE23319
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

#-----------------------------------------------------------
# 8. Save results
#-----------------------------------------------------------

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