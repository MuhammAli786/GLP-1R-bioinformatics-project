############################################################
# GSE236133
# GO Biological Process and KEGG enrichment
############################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(readr)
library(dplyr)
library(ggplot2)

#-----------------------------------------------------------
# 1. Paths
#-----------------------------------------------------------

de_directory <- "results/GSE236133/DE"

output_directory <- "results/GSE236133/enrichment/GO_KEGG"

figure_directory <- "figures/GSE236133/enrichment/GO_KEGG"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

#-----------------------------------------------------------
# 2. Find DEG files
#-----------------------------------------------------------

deg_files <- list.files(
  de_directory,
  pattern = "_DEGs\\.csv$",
  full.names = TRUE
)

if (length(deg_files) == 0) {
  stop(
    "No DEG files found. Run Script 04 first."
  )
}

#-----------------------------------------------------------
# 3. Helper for saving enrichment results
#-----------------------------------------------------------

save_enrichment_result <- function(
    enrichment_object,
    comparison_name,
    gene_direction,
    database_name
) {
  
  output_prefix <- paste(
    comparison_name,
    gene_direction,
    database_name,
    sep = "_"
  )
  
  saveRDS(
    enrichment_object,
    file.path(
      output_directory,
      paste0(output_prefix, ".rds")
    )
  )
  
  enrichment_table <- as.data.frame(
    enrichment_object
  )
  
  write_csv(
    enrichment_table,
    file.path(
      output_directory,
      paste0(output_prefix, ".csv")
    )
  )
  
  if (nrow(enrichment_table) == 0) {
    
    cat(
      "No significant",
      database_name,
      "terms for",
      comparison_name,
      gene_direction,
      "\n"
    )
    
    return(invisible(NULL))
  }
  
  number_to_show <- min(
    15,
    nrow(enrichment_table)
  )
  
  enrichment_plot <- dotplot(
    enrichment_object,
    showCategory = number_to_show
  ) +
    ggtitle(
      paste(
        comparison_name,
        gene_direction,
        database_name
      )
    )
  
  ggsave(
    file.path(
      figure_directory,
      paste0(
        output_prefix,
        "_dotplot.pdf"
      )
    ),
    enrichment_plot,
    width = 10,
    height = 7
  )
  
  ggsave(
    file.path(
      figure_directory,
      paste0(
        output_prefix,
        "_dotplot.png"
      )
    ),
    enrichment_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

#-----------------------------------------------------------
# 4. Run enrichment for one DEG file
#-----------------------------------------------------------

run_enrichment <- function(deg_file) {
  
  deg <- read_csv(
    deg_file,
    show_col_types = FALSE
  )
  
  comparison_name <- sub(
    "_DEGs\\.csv$",
    "",
    basename(deg_file)
  )
  
  cat("\n=================================================\n")
  cat("Enrichment:", comparison_name, "\n")
  cat("=================================================\n")
  
  full_file <- file.path(
    de_directory,
    paste0(
      comparison_name,
      "_all_genes.csv"
    )
  )
  
  if (!file.exists(full_file)) {
    stop(
      "Full result file missing for ",
      comparison_name
    )
  }
  
  full_result <- read_csv(
    full_file,
    show_col_types = FALSE
  )
  
  universe_entrez <- full_result %>%
    filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID) %>%
    as.character()
  
  gene_sets <- list(
    All = deg,
    Up = deg %>% filter(direction == "Up"),
    Down = deg %>% filter(direction == "Down")
  )
  
  all_results <- list()
  
  for (gene_direction in names(gene_sets)) {
    
    selected_degs <- gene_sets[[gene_direction]]
    
    selected_entrez <- selected_degs %>%
      filter(
        !is.na(ENTREZID),
        ENTREZID != ""
      ) %>%
      distinct(ENTREZID) %>%
      pull(ENTREZID) %>%
      as.character()
    
    cat(
      gene_direction,
      "mapped Entrez genes:",
      length(selected_entrez),
      "\n"
    )
    
    if (length(selected_entrez) < 5) {
      
      cat(
        "Skipping",
        gene_direction,
        "because fewer than five mapped genes are available.\n"
      )
      
      next
    }
    
    go_bp <- enrichGO(
      gene = selected_entrez,
      universe = universe_entrez,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 5,
      maxGSSize = 500,
      readable = TRUE
    )
    
    kegg <- enrichKEGG(
      gene = selected_entrez,
      universe = universe_entrez,
      organism = "mmu",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 5,
      maxGSSize = 500
    )
    
    save_enrichment_result(
      go_bp,
      comparison_name,
      gene_direction,
      "GO_BP"
    )
    
    save_enrichment_result(
      kegg,
      comparison_name,
      gene_direction,
      "KEGG"
    )
    
    all_results[[gene_direction]] <- list(
      GO_BP = go_bp,
      KEGG = kegg
    )
  }
  
  saveRDS(
    all_results,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_GO_KEGG_all_results.rds"
      )
    )
  )
}

#-----------------------------------------------------------
# 5. Run every comparison
#-----------------------------------------------------------

invisible(
  lapply(
    deg_files,
    run_enrichment
  )
)

cat("\nScript 06 complete.\n")