############################################################
# GSE23317 - GO Cnet Plot
############################################################

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Check that the GO enrichment object exists
if (!exists("ego")) {
  stop("The GO enrichment object 'ego' was not found. Run the GO enrichment script first.")
}

# Check that the annotated DEG table exists
if (!exists("deg_annot")) {
  stop("The annotated DEG table 'deg_annot' was not found. Run the GO enrichment script first.")
}

# Create named fold-change vector
geneList <- deg_annot$logFC
names(geneList) <- deg_annot$Gene.Symbol

geneList <- sort(geneList, decreasing = TRUE)

# Save figure
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