# GSE236133 Reactome pathway enrichment; input: results/GSE236133/DE/*_DEGs.csv -> output: results/GSE236133/enrichment/Reactome/

library(ReactomePA)
library(enrichplot)
library(readr)
library(dplyr)
library(ggplot2)

# Paths
de_directory <- "results/GSE236133/DE"

output_directory <- "results/GSE236133/enrichment/Reactome"

figure_directory <- "figures/GSE236133/enrichment/Reactome"

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

# Identify DEG files
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

# Save a Reactome result
save_reactome_result <- function(
    enrichment_object,
    comparison_name,
    gene_direction
) {
  
  output_prefix <- paste(
    comparison_name,
    gene_direction,
    "Reactome",
    sep = "_"
  )
  
  saveRDS(
    enrichment_object,
    file.path(
      output_directory,
      paste0(output_prefix, ".rds")
    )
  )
  
  result_table <- as.data.frame(
    enrichment_object
  )
  
  write_csv(
    result_table,
    file.path(
      output_directory,
      paste0(output_prefix, ".csv")
    )
  )
  
  if (nrow(result_table) == 0) {
    
    cat(
      "No significant Reactome pathways for",
      comparison_name,
      gene_direction,
      "\n"
    )
    
    return(invisible(NULL))
  }
  
  pathway_plot <- dotplot(
    enrichment_object,
    showCategory = min(
      15,
      nrow(result_table)
    )
  ) +
    ggtitle(
      paste(
        comparison_name,
        gene_direction,
        "Reactome"
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
    pathway_plot,
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
    pathway_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# Run one comparison
run_reactome <- function(deg_file) {
  
  comparison_name <- sub(
    "_DEGs\\.csv$",
    "",
    basename(deg_file)
  )
  
  deg <- read_csv(
    deg_file,
    show_col_types = FALSE
  )
  
  full_file <- file.path(
    de_directory,
    paste0(
      comparison_name,
      "_all_genes.csv"
    )
  )
  
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
  
  comparison_results <- list()
  
  for (gene_direction in names(gene_sets)) {
    
    selected_entrez <- gene_sets[[gene_direction]] %>%
      filter(
        !is.na(ENTREZID),
        ENTREZID != ""
      ) %>%
      distinct(ENTREZID) %>%
      pull(ENTREZID) %>%
      as.character()
    
    cat(
      "\n",
      comparison_name,
      gene_direction,
      "mapped genes:",
      length(selected_entrez),
      "\n"
    )
    
    if (length(selected_entrez) < 5) {
      
      cat(
        "Skipping because fewer than five mapped genes are available.\n"
      )
      
      next
    }
    
    # Reactome: pvalueCutoff = 0.05, qvalueCutoff = 0.20
    reactome_result <- enrichPathway(
      gene = selected_entrez,
      universe = universe_entrez,
      organism = "mouse",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 5,
      maxGSSize = 500,
      readable = TRUE
    )
    
    save_reactome_result(
      reactome_result,
      comparison_name,
      gene_direction
    )
    
    comparison_results[[gene_direction]] <- reactome_result
  }
  
  saveRDS(
    comparison_results,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_Reactome_all_results.rds"
      )
    )
  )
}

# Run all comparisons
invisible(
  lapply(
    deg_files,
    run_reactome
  )
)

cat("\nScript 07 complete.\n")