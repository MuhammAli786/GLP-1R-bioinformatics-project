# GSE23319 GO BP and KEGG enrichment; input: results/GSE23319_3h_DEGs.csv -> output: results/GO_BP_GSE23319_3h.csv, results/KEGG_GSE23319_3h.csv

library(clusterProfiler)
library(org.Mm.eg.db)

# Load dataset
gse <- readRDS("data/GSE23319.rds")
eset <- gse[[1]]

# Platform annotation
annot <- fData(eset)

# DEG table
deg <- read.csv(
  "results/GSE23319_3h_DEGs.csv",
  row.names = 1
)

# Add probe IDs as a column
deg$ID <- rownames(deg)

# Merge DEGs with annotation
deg_annot <- merge(
  deg,
  annot[, c("ID", "Symbol", "Entrez_Gene_ID")],
  by = "ID"
)

# Remove missing annotations
deg_annot <- deg_annot[
  deg_annot$Symbol != "" &
    !is.na(deg_annot$Symbol),
]

# Remove duplicate genes
deg_annot <- deg_annot[
  !duplicated(deg_annot$Symbol),
]

# Convert Entrez IDs to numeric
entrez <- unique(as.character(deg_annot$Entrez_Gene_ID))
entrez <- entrez[entrez != ""]
entrez <- entrez[!is.na(entrez)]

length(entrez)


# GO enrichment: pvalueCutoff = 0.05, qvalueCutoff = 0.05
ego <- enrichGO(
  gene = entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(ego))

write.csv(
  as.data.frame(ego),
  "results/GO_BP_GSE23319_3h.csv",
  row.names = FALSE
)


# KEGG enrichment: pvalueCutoff = 0.05
ekegg <- enrichKEGG(
  gene = entrez,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(as.data.frame(ekegg))


write.csv(
  as.data.frame(ekegg),
  "results/KEGG_GSE23319_3h.csv",
  row.names = FALSE
)


# Dotplots
if (!dir.exists("figures")) {
  dir.create("figures")
}

library(enrichplot)

pdf("figures/GO_BP_dotplot_GSE23319.pdf", width = 10, height = 8)
dotplot(ego, showCategory = 20)
dev.off()

pdf("figures/KEGG_dotplot_GSE23319.pdf", width = 10, height = 8)
dotplot(ekegg, showCategory = 20)
dev.off()

saveRDS(
  ego,
  file = "results/GSE23319_GO_BP.rds"
)