#!/usr/bin/env Rscript
# GEO2R-style limma differential expression for Illumina GPL6885 arrays
# (GSE71850/Acc10, GSE41345/Acc11, GSE39586/Acc12).
#
# Workflow: log2 auto-transform -> lmFit(~0+group) -> makeContrasts ->
# contrasts.fit -> eBayes(0.01) -> topTable(adjust="fdr", number=Inf)
#
# Inputs  (per dataset, prepared by the Python driver):
#   <work>/<tag>_exprs.tsv      probe x sample matrix (ID_REF + GSMs)
#   <work>/<tag>_groups.tsv     sample <tab> group
#   <work>/<tag>_contrasts.tsv  one contrast per line, e.g. "Exendin4-Control"
# Output:
#   <work>/<tag>__<contrast>__topTable.tsv  (ID, logFC, adj.P.Val, P.Value, t, B)

suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
work <- args[1]
tag  <- args[2]

ex <- as.matrix(read.delim(file.path(work, paste0(tag, "_exprs.tsv")),
                           row.names = 1, check.names = FALSE))
grp <- read.delim(file.path(work, paste0(tag, "_groups.tsv")),
                  header = TRUE, check.names = FALSE)
cts <- readLines(file.path(work, paste0(tag, "_contrasts.tsv")))
cts <- cts[nzchar(cts)]

# order expression columns to match group table
ex <- ex[, as.character(grp$sample), drop = FALSE]

# GEO2R log2 auto-transform
qx <- as.numeric(quantile(ex, c(0, .25, .5, .75, .99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[ex <= 0] <- NaN; ex <- log2(ex) }

# drop probes with missing values (as GEO2R: complete.cases)
ex <- ex[complete.cases(ex), , drop = FALSE]

gs <- factor(grp$group)
design <- model.matrix(~ 0 + gs)
colnames(design) <- levels(gs)

fit <- lmFit(ex, design)
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

for (i in seq_along(cts)) {
  tT <- topTable(fit2, coef = i, adjust = "fdr", sort.by = "B", number = Inf)
  tT$ID <- rownames(tT)
  out <- tT[, c("ID", "logFC", "P.Value", "adj.P.Val", "t", "B")]
  safe <- gsub("[^A-Za-z0-9]+", "_", cts[i])
  write.table(out, file = file.path(work, paste0(tag, "__", safe, "__topTable.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("wrote", tag, cts[i], nrow(out), "probes\n")
}
