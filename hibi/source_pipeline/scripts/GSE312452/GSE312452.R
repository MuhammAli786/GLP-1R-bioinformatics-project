library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(EnhancedVolcano)

counts <- read.table(
  "GSE312452_read_count_all_groups.txt.gz",
  header = TRUE,
  sep = "\t",
  check.names = FALSE)

counts <- counts[
  !is.na(counts$`Gene ID`),
]

rownames(counts) <- make.unique(
  as.character(counts$`Gene ID`)
)

counts <- counts[, -1]

head(counts)

samples <- colnames(counts)

metadata <- data.frame(
  row.names = samples,
  genotype = c(
    "KO","KO","KO","KO",
    "KO","KO","KO","KO",
    "WT","WT","WT","WT",
    "WT","WT","WT"
  ),
  injury = c(
    "HI","HI","HI","HI",
    "Sham","Sham","Sham","Sham",
    "HI","HI","HI","HI",
    "Sham","Sham","Sham"
  )
)

metadata_WT <- metadata[
  metadata$genotype == "WT",
]


counts_WT <- counts[
  ,
  rownames(metadata_WT)
]


all(colnames(counts_WT)==rownames(metadata_WT))


sum(is.na(counts_WT))

counts_WT[is.na(counts_WT)] <- 0

sum(is.na(counts_WT))

counts_WT <- round(counts_WT)

metadata_WT$injury <- factor(
  metadata_WT$injury,
  levels = c("Sham","HI")
)

# DESeq2 analysis

dds <- DESeqDataSetFromMatrix(
  countData = counts_WT,
  colData = metadata_WT,
  design = ~ injury
)

dds <- dds[
  rowSums(counts(dds)) > 10,
]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "injury",
    "HI",
    "Sham"
  )
)

res <- as.data.frame(res)

deg_sig <- res[
  !is.na(res$padj) &
    res$padj < 0.05,]

dim(deg_sig)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  title =
    "GSE312452: WT HI vs WT Sham pericytes"
)

gene_list <- rownames(deg_sig)


gene_entrez <- bitr(
  gene_list,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Mm.eg.db
)

#GO Plot
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(
  ego,
  showCategory = 20,
  title = "GO Biological Processes: WT HI vs WT Sham"
)

#KEGG Plot
ekegg <- enrichKEGG(
  gene = gene_list,
  organism = "mmu",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05
)


dotplot(
  ekegg,
  showCategory = 20,
  title = "KEGG Pathways: WT HI vs WT Sham"
)


top_up <- head(
  deg_sig[order(-deg_sig$log2FoldChange),],
  10
)


top_down <- head(
  deg_sig[order(deg_sig$log2FoldChange),],
  10
)


top_genes <- rbind(
  top_up,
  top_down
)

#Bar Plot
ggplot(
  top_genes,
  aes(
    x = reorder(rownames(top_genes), log2FoldChange),
    y = log2FoldChange
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top Differentially Expressed Genes WT HI vs WT Sham",
    x = "Gene",
    y = "log2 Fold Change"
  )


gene_fc <- deg_sig$log2FoldChange

names(gene_fc) <- rownames(deg_sig)

#Cnet Plot
cnetplot(
  ego,
  showCategory = 10,
  foldChange = gene_fc,
  node_label = "all"
)