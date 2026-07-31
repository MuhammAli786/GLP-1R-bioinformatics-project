############################################################
# GSE144456 - Differential Expression
# Comparison: P5 hypoxia-ischemia versus control at 3 hours
#
# Two-colour microarray:
# Each selected array already contains an ischemic-versus-
# control log-ratio, so an intercept-only limma model is used.
############################################################

library(GEOquery)
library(limma)

#-----------------------------------------------------------
# 1. Load GEO dataset
#-----------------------------------------------------------

gse_file <- "data/GSE144456.rds"

if (!file.exists(gse_file)) {
  stop(
    "GSE144456 RDS file not found: ",
    gse_file,
    "\nRun script 01 first."
  )
}

gse <- readRDS(gse_file)

eset <- gse[[1]]
expr <- exprs(eset)

cat("\nPlatform:\n")
print(annotation(eset))

cat("\nFull expression-matrix dimensions:\n")
print(dim(expr))

#-----------------------------------------------------------
# 2. Load samples selected in script 03
#-----------------------------------------------------------

metadata_file <-
  "metadata/GSE144456_P5_3h_samples.rds"

if (!file.exists(metadata_file)) {
  stop(
    "Prepared metadata file not found: ",
    metadata_file,
    "\nRun script 03 first."
  )
}

selected_metadata <- readRDS(metadata_file)

selected_samples <- as.character(
  selected_metadata$sample_id
)

if (!all(selected_samples %in% colnames(expr))) {
  missing_samples <- setdiff(
    selected_samples,
    colnames(expr)
  )
  
  stop(
    "These selected samples are missing from the ",
    "expression matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

expr_selected <- expr[
  ,
  selected_samples,
  drop = FALSE
]

cat("\nSelected samples:\n")
print(colnames(expr_selected))

cat("\nSelected expression-matrix dimensions:\n")
print(dim(expr_selected))

cat("\nExpression-value summary:\n")
print(summary(as.vector(expr_selected)))

#-----------------------------------------------------------
# 3. Check that values look like processed log-ratios
#-----------------------------------------------------------

if (
  min(expr_selected, na.rm = TRUE) >= 0
) {
  warning(
    "All selected expression values are non-negative. ",
    "Confirm that the GEO matrix contains log-ratios ",
    "before using an intercept-only model."
  )
}

#-----------------------------------------------------------
# 4. Filter uninformative probes
#-----------------------------------------------------------

# Remove probes that are all missing
keep_not_all_na <- rowSums(
  !is.na(expr_selected)
) > 0

expr_filtered <- expr_selected[
  keep_not_all_na,
  ,
  drop = FALSE
]

# Remove probes equal to zero across all three arrays
keep_nonzero <- rowSums(
  expr_filtered != 0,
  na.rm = TRUE
) > 0

expr_filtered <- expr_filtered[
  keep_nonzero,
  ,
  drop = FALSE
]

cat(
  "\nProbes before filtering:",
  nrow(expr_selected),
  "\n"
)

cat(
  "Probes after filtering:",
  nrow(expr_filtered),
  "\n"
)

#-----------------------------------------------------------
# 5. Fit intercept-only limma model
#-----------------------------------------------------------

# Channel 1 = control
# Channel 2 = ischemic
#
# Assuming GEO stores the processed ratio as:
# log2(channel 2 / channel 1)
#
# Positive logFC = higher expression in ischemic tissue
# Negative logFC = lower expression in ischemic tissue

design <- matrix(
  1,
  nrow = ncol(expr_filtered),
  ncol = 1
)

rownames(design) <- colnames(expr_filtered)
colnames(design) <- "HI_vs_Control_P5_3h"

cat("\nDesign matrix:\n")
print(design)

fit <- lmFit(
  expr_filtered,
  design
)

fit <- eBayes(fit)

#-----------------------------------------------------------
# 6. Extract all limma results
#-----------------------------------------------------------

de_all <- topTable(
  fit,
  coef = "HI_vs_Control_P5_3h",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

de_all$probe_id <- rownames(de_all)

de_all <- de_all[
  ,
  c(
    "probe_id",
    setdiff(
      colnames(de_all),
      "probe_id"
    )
  )
]

#-----------------------------------------------------------
# 7. Add platform annotation
#-----------------------------------------------------------

feature_annotation <- fData(eset)

feature_annotation$probe_id <-
  rownames(feature_annotation)

de_annotated <- merge(
  de_all,
  feature_annotation,
  by = "probe_id",
  all.x = TRUE,
  sort = FALSE
)

# Restore limma ranking after merge
de_annotated <- de_annotated[
  match(
    de_all$probe_id,
    de_annotated$probe_id
  ),
  ,
  drop = FALSE
]

# Verify ordering was preserved
stopifnot(
  identical(
    de_annotated$probe_id,
    de_all$probe_id
  )
)

#-----------------------------------------------------------
# 8. Apply DEG thresholds
#-----------------------------------------------------------

# Same thresholds used in the previous datasets:
# |logFC| >= 0.2
# adjusted P-value < 0.05

deg <- de_annotated[
  !is.na(de_annotated$logFC) &
    !is.na(de_annotated$adj.P.Val) &
    abs(de_annotated$logFC) >= 0.2 &
    de_annotated$adj.P.Val < 0.05,
  ,
  drop = FALSE
]

deg$direction <- ifelse(
  deg$logFC > 0,
  "Up",
  "Down"
)

cat("\n============================================\n")
cat("Differential-expression summary\n")
cat("============================================\n")

cat(
  "Total significant DEGs:",
  nrow(deg),
  "\n"
)

cat("\nDirection counts:\n")
print(
  table(
    deg$direction
  )
)

#-----------------------------------------------------------
# 9. Save results
#-----------------------------------------------------------

write.csv(
  de_all,
  "results/GSE144456_P5_3h_all_limma_results.csv",
  row.names = FALSE
)

write.csv(
  de_annotated,
  "results/GSE144456_P5_3h_all_annotated.csv",
  row.names = FALSE
)

write.csv(
  deg,
  "results/GSE144456_P5_3h_DEGs_annotated.csv",
  row.names = FALSE
)

saveRDS(
  fit,
  "results/GSE144456_P5_3h_limma_fit.rds"
)

saveRDS(
  de_annotated,
  "results/GSE144456_P5_3h_all_annotated.rds"
)

saveRDS(
  deg,
  "results/GSE144456_P5_3h_DEGs_annotated.rds"
)

cat("\nDifferential-expression analysis complete.\n")