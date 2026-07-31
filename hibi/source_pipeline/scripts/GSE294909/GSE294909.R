library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(EnhancedVolcano)
library(dplyr)

data <- read.csv(
  "GSE294909_DESeq2_removedbatch_log_2025-04-15.csv.gz",
  check.names = FALSE
)

colnames(data)[1] <- "gene"

rownames(data) <- data$gene

data <- data[, -1]

data <- as.matrix(data)


group <- factor(
  c(
    rep("Sham",8),
    rep("LPS",6),
    rep("NaCl_HI",6),
    rep("LPS_HI",3)
  )
)


design <- model.matrix(~0 + group)

colnames(design) <- levels(group)


fit <- lmFit(
  data,
  design
)


contrast <- makeContrasts(
  LPS_HI_vs_NaCl_HI = LPS_HI - NaCl_HI,
  levels = design
)


fit2 <- contrasts.fit(
  fit,
  contrast
)


fit2 <- eBayes(fit2)


res <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)


head(res)


deg_sig <- res[
  res$adj.P.Val < 0.05,
]


dim(deg_sig)

#VolcanoPlot
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "logFC",
  y = "adj.P.Val",
  title = "GSE294909: LPS_HI vs NaCl_HI microglia"
)

#Bar Plot
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
    x=reorder(rownames(top_genes),logFC),
    y=logFC
  )
)+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_classic()+
  labs(
    title="Top DEGs LPS_HI vs NaCl_HI",
    x="Gene",
    y="log2 Fold Change"
  )

#GO Plot
gene_list <- rownames(deg_sig)


gene_entrez <- bitr(
  gene_list,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Rn.eg.db
)


ego <- enrichGO(
  gene=gene_entrez$ENTREZID,
  OrgDb=org.Rn.eg.db,
  ont="BP",
  pAdjustMethod="BH",
  pvalueCutoff=0.05,
  readable=TRUE
)


dotplot(
  ego,
  showCategory=20,
  title="GO Biological Processes LPS_HI vs NaCl_HI"
)

#KEGG Plot 
ekegg <- enrichKEGG(
  gene=gene_entrez$ENTREZID,
  organism="rno",
  pvalueCutoff=0.05
)


dotplot(
  ekegg,
  showCategory=20,
  title="KEGG Pathways LPS_HI vs NaCl_HI"
)

#Cnet Plot
gene_fc <- deg_sig$logFC

names(gene_fc) <- rownames(deg_sig)

cnetplot(
  ego,
  showCategory=10,
  foldChange=gene_fc,
  node_label="all")