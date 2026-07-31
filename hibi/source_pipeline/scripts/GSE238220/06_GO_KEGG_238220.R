############################################################
# GSE238220
# GO Biological Process and KEGG enrichment
#
# Enrichment is run separately for:
#   All significant DEGs
#   Upregulated DEGs in HI
#   Downregulated DEGs in HI
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

de_directory <- "results/GSE238220/DE"

output_directory <- paste0(
  "results/GSE238220/enrichment/",
  "GO_KEGG"
)

figure_directory <- paste0(
  "figures/GSE238220/enrichment/",
  "GO_KEGG"
)

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
    "No DEG files were found in:\n",
    de_directory,
    "\nRun Script 04 first."
  )
}

cat("\nDEG files found:\n")
print(deg_files)

#-----------------------------------------------------------
# 3. Function to save enrichment output
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
  
  rds_file <- file.path(
    output_directory,
    paste0(output_prefix, ".rds")
  )
  
  csv_file <- file.path(
    output_directory,
    paste0(output_prefix, ".csv")
  )
  
  saveRDS(
    enrichment_object,
    rds_file
  )
  
  enrichment_table <- as.data.frame(
    enrichment_object
  )
  
  write_csv(
    enrichment_table,
    csv_file
  )
  
  cat(
    database_name,
    gene_direction,
    "significant terms:",
    nrow(enrichment_table),
    "\n"
  )
  
  if (nrow(enrichment_table) == 0) {
    return(invisible(NULL))
  }
  
  number_to_show <- min(
    15,
    nrow(enrichment_table)
  )
  
  dotplot_object <- enrichplot::dotplot(
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
  
  pdf_file <- file.path(
    figure_directory,
    paste0(
      output_prefix,
      "_dotplot.pdf"
    )
  )
  
  png_file <- file.path(
    figure_directory,
    paste0(
      output_prefix,
      "_dotplot.png"
    )
  )
  
  ggsave(
    filename = pdf_file,
    plot = dotplot_object,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = png_file,
    plot = dotplot_object,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  invisible(dotplot_object)
}

#-----------------------------------------------------------
# 4. Function to run enrichment for one comparison
#-----------------------------------------------------------

run_enrichment <- function(deg_file) {
  
  comparison_name <- sub(
    "_DEGs\\.csv$",
    "",
    basename(deg_file)
  )
  
  cat("\n\n=================================================\n")
  cat("Running enrichment:", comparison_name, "\n")
  cat("=================================================\n")
  
  full_result_file <- file.path(
    de_directory,
    paste0(
      comparison_name,
      "_all_genes.csv"
    )
  )
  
  if (!file.exists(full_result_file)) {
    stop(
      "Full result file was not found:\n",
      full_result_file
    )
  }
  
  deg <- read_csv(
    deg_file,
    show_col_types = FALSE
  )
  
  full_result <- read_csv(
    full_result_file,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "ENTREZID",
    "direction"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(deg)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      comparison_name,
      " DEG file is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  #---------------------------------------------------------
  # Background universe
  #---------------------------------------------------------
  
  universe_entrez <- full_result %>%
    dplyr::filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) %>%
    dplyr::distinct(ENTREZID) %>%
    dplyr::pull(ENTREZID) %>%
    as.character()
  
  cat(
    "\nMapped background genes:",
    length(universe_entrez),
    "\n"
  )
  
  #---------------------------------------------------------
  # Define gene sets
  #---------------------------------------------------------
  
  gene_sets <- list(
    All = deg,
    
    Up = deg %>%
      dplyr::filter(
        direction == "Up"
      ),
    
    Down = deg %>%
      dplyr::filter(
        direction == "Down"
      )
  )
  
  combined_results <- list()
  
  #---------------------------------------------------------
  # Run each direction
  #---------------------------------------------------------
  
  for (gene_direction in names(gene_sets)) {
    
    selected_degs <- gene_sets[[gene_direction]]
    
    selected_entrez <- selected_degs %>%
      dplyr::filter(
        !is.na(ENTREZID),
        ENTREZID != ""
      ) %>%
      dplyr::distinct(ENTREZID) %>%
      dplyr::pull(ENTREZID) %>%
      as.character()
    
    cat(
      "\n",
      gene_direction,
      "mapped DEG genes:",
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
    
    #-------------------------------------------------------
    # GO Biological Process
    #-------------------------------------------------------
    
    go_bp <- clusterProfiler::enrichGO(
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
    
    #-------------------------------------------------------
    # KEGG
    #-------------------------------------------------------
    
    kegg <- tryCatch(
      
      clusterProfiler::enrichKEGG(
        gene = selected_entrez,
        universe = universe_entrez,
        organism = "mmu",
        keyType = "ncbi-geneid",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        minGSSize = 5,
        maxGSSize = 500
      ),
      
      error = function(e) {
        
        warning(
          comparison_name,
          " ",
          gene_direction,
          ": KEGG enrichment failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    save_enrichment_result(
      enrichment_object = go_bp,
      comparison_name = comparison_name,
      gene_direction = gene_direction,
      database_name = "GO_BP"
    )
    
    if (!is.null(kegg)) {
      
      save_enrichment_result(
        enrichment_object = kegg,
        comparison_name = comparison_name,
        gene_direction = gene_direction,
        database_name = "KEGG"
      )
    }
    
    combined_results[[gene_direction]] <- list(
      GO_BP = go_bp,
      KEGG = kegg
    )
  }
  
  saveRDS(
    combined_results,
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