# Annotate and enrichment analysis for GSE36215 (HI vs Control)
# Input: GSE36215.top.table.tsv, GPL6246-18741.txt -> Output: GSE36215_HI_vs_Control_annotated.csv, GO/KEGG dotplots, bar, cnetplot
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)

deg <- read.table(
  "GSE36215.top.table.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

anno <- read.delim(
  "GPL6246-18741.txt",
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

anno_unique <- anno[, c(
  "ID",
  "gene_assignment"
)]

anno_unique <- anno_unique[
  anno_unique$gene_assignment != "",
]

anno_unique <- anno_unique[
  !duplicated(anno_unique$ID),
]

deg_annot <- merge(
  deg,
  anno_unique,
  by = "ID"
)

deg_annot$GENE_SYMBOL <- sapply(
  strsplit(
    deg_annot$gene_assignment,
    " // "
  ),
  function(x){
    if(length(x) >= 2){
      x[2]
    } else {
      NA
    }
  }
)

deg_annot <- deg_annot[
  !is.na(deg_annot$GENE_SYMBOL) &
    deg_annot$GENE_SYMBOL != "",
]

write.csv(
  deg_annot,
  "GSE36215_HI_vs_Control_annotated.csv",
  row.names = FALSE
)

deg_sig <- subset(
  deg_annot,
  adj.P.Val < 0.05
)

gene_list <- unique(
  deg_sig$GENE_SYMBOL
)

gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(
  ego,
  showCategory = 20
)

ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

dotplot(
  ekegg,
  showCategory = 20
)

top_up <- head(
  deg_sig[order(-deg_sig$logFC),],
  10
)

top_down <- head(
  deg_sig[order(deg_sig$logFC),],
  10
)

top_genes <- rbind(
  top_up,
  top_down
)

ggplot(
  top_genes,
  aes(
    x = reorder(GENE_SYMBOL, logFC),
    y = logFC
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top Differentially Expressed Genes",
    x = "Gene",
    y = "log2 Fold Change"
  )

gene_fc <- deg_sig[,c(
  "GENE_SYMBOL",
  "logFC"
)]

gene_fc <- gene_fc[
  !duplicated(gene_fc$GENE_SYMBOL),
]

fold_change <- gene_fc$logFC
names(fold_change) <- gene_fc$GENE_SYMBOL

cnetplot(
  ego,
  showCategory = 10,
  foldChange = fold_change,
  node_label = "all")