# Create GO cnet plot for GSE23317 3h HI vs Sham
# Input: requires `ego` and `deg_annot` from script 06 -> Output: figures/GSE23317_cnetplot.pdf

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

if (!exists("ego")) {
  stop("The GO enrichment object 'ego' was not found. Run the GO enrichment script first.")
}

if (!exists("deg_annot")) {
  stop("The annotated DEG table 'deg_annot' was not found. Run the GO enrichment script first.")
}

geneList <- deg_annot$logFC
names(geneList) <- deg_annot$Gene.Symbol

geneList <- sort(geneList, decreasing = TRUE)

pdf("figures/GSE23317_cnetplot.pdf",
    width = 12,
    height = 10)

cnetplot(
  ego,
  showCategory = 5,
  foldChange = geneList,
  node_label = "all"
)

dev.off()