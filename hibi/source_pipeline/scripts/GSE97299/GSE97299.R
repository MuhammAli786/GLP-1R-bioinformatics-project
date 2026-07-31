BiocManager::install("org.Rn.eg.db")
BiocManager::install("GOplot")

library(GOplot)
library(dplyr)
library(org.Rn.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(enrichplot)

deg <- read.table(
  "GSE97299.top.table.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE)

anno <- read.delim(
  "GPL22388-48305.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE)

anno_clean <- anno[, c(
  "probeset_id",
  "gene_assignment")]


anno_clean <- anno_clean[
  anno_clean$gene_assignment != "",]

deg_annot <- merge(
  deg,
  anno_clean,
  by.x = "ID",
  by.y = "probeset_id")

deg_annot$GENE_SYMBOL <- sapply(
  strsplit(
    deg_annot$gene_assignment,
    " // "),
  function(x){if(length(x) >= 2){
      return(x[2])
    } else {
      return(NA)}})

deg_annot <- deg_annot[
  !is.na(deg_annot$GENE_SYMBOL) &
    deg_annot$GENE_SYMBOL != "",]

deg_sig <- subset(
  deg_annot,
  adj.P.Val < 0.05)

dim(deg_sig)

gene_list <- unique(
  deg_sig$GENE_SYMBOL)

gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db)

#GO Plot
ego <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Rn.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)

dotplot(
  ego,
  showCategory = 20)

up_genes <- unique(
  deg_annot$GENE_SYMBOL[
    deg_annot$logFC > 0])


down_genes <- unique(
  deg_annot$GENE_SYMBOL[
    deg_annot$logFC < 0])

write.csv(
  deg_annot,
  "GSE144455_RV_PBS_12h_DEG_annotated.csv",
  row.names = FALSE)

write.csv(
  gene_list,
  "GSE144455_RV_PBS_12h_gene_symbols.csv",
  row.names = FALSE)

gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db)

#GO Plot
ego <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)

dotplot(
  ego,
  showCategory = 20)

#KEGG
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "rno",
  pvalueCutoff = 0.05)

dotplot(
  ekegg,
  showCategory = 20)

# Bar plot of log2 fold change
top_up <- head(
  deg_sig[order(-deg_sig$logFC),],
  10)

top_down <- head(
  deg_sig[order(deg_sig$logFC),],
  10)

top_genes <- rbind(
  top_up,
  top_down)

ggplot(
  top_genes,
  aes(
    x = reorder(GENE_SYMBOL, logFC),
    y = logFC))+
  geom_bar(
    stat="identity")+
  coord_flip()+
  theme_classic()+
  labs(
    title="Top Differentially Expressed Genes HI SVZ P10 vs Sham SVZ P10",
    x="Gene",
    y="log2 Fold Change")

#Cnet Plot
gene_fc <- deg_sig[,c(
  "GENE_SYMBOL",
  "logFC")]

gene_fc <- gene_fc[
  !duplicated(gene_fc$GENE_SYMBOL),]

fold_change <- gene_fc$logFC
names(fold_change) <- gene_fc$GENE_SYMBOL

cnetplot(
  ego,
  showCategory = 10,
  foldChange = fold_change,
  node_label = "all")
