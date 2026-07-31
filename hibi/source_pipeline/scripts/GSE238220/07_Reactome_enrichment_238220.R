############################################################
# GSE238220
# Reactome pathway enrichment
#
# Enrichment is run separately for:
#   All significant DEGs
#   Upregulated DEGs in HI
#   Downregulated DEGs in HI
############################################################

library(ReactomePA)
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
  "Reactome"
)

figure_directory <- paste0(
  "figures/GSE238220/enrichment/",
  "Reactome"
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
    "No DEG files were found. Run Script 04 first."
  )
}

#-----------------------------------------------------------
# 3. Save one Reactome result
#-----------------------------------------------------------

save_reactome_result <- function(
    reactome_object,
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
    reactome_object,
    file.path(
      output_directory,
      paste0(
        output_prefix,
        ".rds"
      )
    )
  )
  
  result_table <- as.data.frame(
    reactome_object
  )
  
  write_csv(
    result_table,
    file.path(
      output_directory,
      paste0(
        output_prefix,
        ".csv"
      )
    )
  )
  
  cat(
    gene_direction,
    "significant Reactome pathways:",
    nrow(result_table),
    "\n"
  )
  
  if (nrow(result_table) == 0) {
    return(invisible(NULL))
  }
  
  dotplot_object <- enrichplot::dotplot(
    reactome_object,
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
    filename = file.path(
      figure_directory,
      paste0(
        output_prefix,
        "_dotplot.pdf"
      )
    ),
    plot = dotplot_object,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = file.path(
      figure_directory,
      paste0(
        output_prefix,
        "_dotplot.png"
      )
    ),
    plot = dotplot_object,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  invisible(dotplot_object)
}

#-----------------------------------------------------------
# 4. Run one comparison
#-----------------------------------------------------------

run_reactome <- function(deg_file) {
  
  comparison_name <- sub(
    "_DEGs\\.csv$",
    "",
    basename(deg_file)
  )
  
  cat("\n\n=================================================\n")
  cat("Reactome:", comparison_name, "\n")
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
  
  universe_entrez <- full_result %>%
    dplyr::filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) %>%
    dplyr::distinct(ENTREZID) %>%
    dplyr::pull(ENTREZID) %>%
    as.character()
  
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
  
  all_results <- list()
  
  for (gene_direction in names(gene_sets)) {
    
    selected_entrez <- gene_sets[[gene_direction]] %>%
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
    
    reactome_result <- tryCatch(
      
      ReactomePA::enrichPathway(
        gene = selected_entrez,
        universe = universe_entrez,
        organism = "mouse",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        minGSSize = 5,
        maxGSSize = 500,
        readable = TRUE
      ),
      
      error = function(e) {
        
        warning(
          comparison_name,
          " ",
          gene_direction,
          ": Reactome enrichment failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (is.null(reactome_result)) {
      next
    }
    
    save_reactome_result(
      reactome_object = reactome_result,
      comparison_name = comparison_name,
      gene_direction = gene_direction
    )
    
    all_results[[gene_direction]] <- reactome_result
  }
  
  saveRDS(
    all_results,
    file.path(
      output_directory,
      paste0(
        comparison_name,
        "_Reactome_all_results.rds"
      )
    )
  )
}

#-----------------------------------------------------------
# 5. Run all comparisons
#-----------------------------------------------------------

invisible(
  lapply(
    deg_files,
    run_reactome
  )
)

cat("\nScript 07 complete.\n")