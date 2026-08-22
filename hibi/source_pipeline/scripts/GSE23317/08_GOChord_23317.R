# Create GOChord plot for GSE23317 3h HI vs Sham
# Input: requires `ego` and `deg_annot` from script 06 -> Output: figures/GOChord.pdf

library(clusterProfiler)
library(GOplot)
library(dplyr)

go_df <- as.data.frame(ego)

terms <- go_df %>%
  mutate(category = "BP") %>%
  dplyr::select(
    category,
    ID,
    term = Description,
    adj_pval = p.adjust,
    genes = geneID
  )

circ <- deg_annot %>%
  dplyr::select(Symbol, logFC) %>%
  mutate(Symbol = toupper(Symbol)) %>%
  rename(ID = Symbol) %>%
  filter(!is.na(ID), ID != "") %>%
  distinct(ID, .keep_all = TRUE)

go_genes <- unique(unlist(strsplit(terms$genes, "/")))

cat("Genes in DEG table:", nrow(circ), "\n")
cat("Genes in GO terms:", length(go_genes), "\n")
cat("Matching genes:", length(intersect(circ$ID, go_genes)), "\n")

circ_data <- circle_dat(
  terms,
  circ
)

summary(circ_data$logFC)

if(all(is.na(circ_data$logFC))){
  stop("No logFC values were mapped. Gene IDs do not match between GO results and DEG table.")
}

pdf(
  "figures/GOChord.pdf",
  width = 12,
  height = 10
)

GOChord(
  circ_data,
  limit = c(5,5)
)

dev.off()

message("GOChord successfully created.")