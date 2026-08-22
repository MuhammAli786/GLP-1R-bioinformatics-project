# GO and KEGG enrichment for GSE23317 3h HI vs Sham
# Input: data/GSE23317.rds, results/GSE23317_3h_DEGs.csv -> Output: results/GO_BP_GSE23317_3h.csv, results/KEGG_GSE23317_3h.csv, results/GSE23317_GO_BP.rds, figures/GO_BP_dotplot_GSE23317.pdf, figures/KEGG_dotplot_23317.pdf

library(clusterProfiler)
library(org.Mm.eg.db)

gse <- readRDS("data/GSE23317.rds")
eset <- gse[[1]]

annot <- fData(eset)

deg <- read.csv(
  "results/GSE23317_3h_DEGs.csv",
  row.names = 1
)

deg$ID <- rownames(deg)

deg_annot <- merge(
  deg,
  annot[, c("ID", "Symbol", "Entrez_Gene_ID")],
  by = "ID"
)

deg_annot <- deg_annot[
  deg_annot$Symbol != "" &
    !is.na(deg_annot$Symbol),
]

deg_annot <- deg_annot[
  !duplicated(deg_annot$Symbol),
]

entrez <- unique(as.character(deg_annot$Entrez_Gene_ID))
entrez <- entrez[entrez != ""]
entrez <- entrez[!is.na(entrez)]

length(entrez)

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
  "results/GO_BP_GSE23317_3h.csv",
  row.names = FALSE
)

ekegg <- enrichKEGG(
  gene = entrez,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(as.data.frame(ekegg))

write.csv(
  as.data.frame(ekegg),
  "results/KEGG_GSE23317_3h.csv",
  row.names = FALSE
)

if (!dir.exists("figures")) {
  dir.create("figures")
}

library(enrichplot)

pdf("figures/GO_BP_dotplot_GSE23317.pdf", width = 10, height = 8)
dotplot(ego, showCategory = 20)
dev.off()

pdf("figures/KEGG_dotplot_23317.pdf", width = 10, height = 8)
dotplot(ekegg, showCategory = 20)
dev.off()

saveRDS(
  ego,
  file = "results/GSE23317_GO_BP.rds"
)