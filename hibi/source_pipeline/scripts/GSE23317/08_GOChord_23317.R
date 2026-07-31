############################################################
## 08_GOChord.R
## GOChord plot for GSE23317 (3 h HI vs Sham)
############################################################

#=========================
# Load packages
#=========================

library(clusterProfiler)
library(GOplot)
library(dplyr)

#=========================
# Convert GO results
#=========================

go_df <- as.data.frame(ego)

#=========================
# Prepare GOplot terms table
#=========================

terms <- go_df %>%
  mutate(category = "BP") %>%
  dplyr::select(
    category,
    ID,
    term = Description,
    adj_pval = p.adjust,
    genes = geneID
  )

#=========================
# Prepare expression table
#=========================

circ <- deg_annot %>%
  dplyr::select(Symbol, logFC) %>%
  mutate(Symbol = toupper(Symbol)) %>%
  rename(ID = Symbol) %>%
  filter(!is.na(ID), ID != "") %>%
  distinct(ID, .keep_all = TRUE)

#=========================
# Check overlap
#=========================

go_genes <- unique(unlist(strsplit(terms$genes, "/")))

cat("Genes in DEG table:", nrow(circ), "\n")
cat("Genes in GO terms:", length(go_genes), "\n")
cat("Matching genes:", length(intersect(circ$ID, go_genes)), "\n")

#=========================
# Create GOplot object
#=========================

circ_data <- circle_dat(
  terms,
  circ
)

#=========================
# Check mapping
#=========================

summary(circ_data$logFC)

if(all(is.na(circ_data$logFC))){
  stop("No logFC values were mapped. Gene IDs do not match between GO results and DEG table.")
}

#=========================
# Plot GOChord
#=========================

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