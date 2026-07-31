############################################################
# GSE238220
# GO Biological Process and Reactome cnetplots
############################################################

library(enrichplot)
library(readr)
library(dplyr)
library(ggplot2)

#-----------------------------------------------------------
# 1. Paths
#-----------------------------------------------------------

de_directory <- "results/GSE238220/DE"

go_directory <- paste0(
  "results/GSE238220/enrichment/",
  "GO_KEGG"
)

reactome_directory <- paste0(
  "results/GSE238220/enrichment/",
  "Reactome"
)

figure_directory <- "figures/GSE238220/cnetplots"

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
# 3. Create named fold-change vector
#-----------------------------------------------------------

make_fold_change_vector <- function(deg_table) {
  
  fold_change_table <- deg_table %>%
    dplyr::filter(
      !is.na(ENTREZID),
      ENTREZID != "",
      !is.na(log2FoldChange)
    ) %>%
    dplyr::arrange(
      padj
    ) %>%
    dplyr::distinct(
      ENTREZID,
      .keep_all = TRUE
    )
  
  fold_change <- fold_change_table$log2FoldChange
  
  names(fold_change) <- as.character(
    fold_change_table$ENTREZID
  )
  
  fold_change
}

#-----------------------------------------------------------
# 4. Save one cnetplot
#-----------------------------------------------------------

save_cnetplot <- function(
    enrichment_object,
    fold_change,
    plot_title,
    output_prefix
) {
  
  result_table <- as.data.frame(
    enrichment_object
  )
  
  if (nrow(result_table) == 0) {
    
    cat(
      "Skipping empty result:",
      output_prefix,
      "\n"
    )
    
    return(invisible(NULL))
  }
  
  categories_to_show <- min(
    5,
    nrow(result_table)
  )
  
  cat(
    "\nCreating cnetplot:",
    output_prefix,
    "\n"
  )
  
  plot_object <- tryCatch(
    
    enrichplot::cnetplot(
      enrichment_object,
      showCategory = categories_to_show,
      foldChange = fold_change,
      layout = "kk"
    ) +
      ggtitle(plot_title),
    
    error = function(first_error) {
      
      message(
        "The requested layout failed. ",
        "Retrying with the default layout."
      )
      
      tryCatch(
        
        enrichplot::cnetplot(
          enrichment_object,
          showCategory = categories_to_show,
          foldChange = fold_change
        ) +
          ggtitle(plot_title),
        
        error = function(second_error) {
          
          warning(
            "Cnetplot failed for ",
            output_prefix,
            ": ",
            conditionMessage(second_error)
          )
          
          NULL
        }
      )
    }
  )
  
  if (is.null(plot_object)) {
    return(invisible(NULL))
  }
  
  pdf_file <- file.path(
    figure_directory,
    paste0(
      output_prefix,
      "_cnetplot.pdf"
    )
  )
  
  png_file <- file.path(
    figure_directory,
    paste0(
      output_prefix,
      "_cnetplot.png"
    )
  )
  
  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = 12,
    height = 10
  )
  
  ggsave(
    filename = png_file,
    plot = plot_object,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  cat("Saved:\n")
  cat(pdf_file, "\n")
  cat(png_file, "\n")
  
  invisible(plot_object)
}

#-----------------------------------------------------------
# 5. Process one comparison
#-----------------------------------------------------------

process_comparison <- function(deg_file) {
  
  comparison_name <- sub(
    "_DEGs\\.csv$",
    "",
    basename(deg_file)
  )
  
  cat("\n\n=================================================\n")
  cat("Cnetplots:", comparison_name, "\n")
  cat("=================================================\n")
  
  deg <- read_csv(
    deg_file,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "ENTREZID",
    "log2FoldChange",
    "padj",
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
  
  gene_directions <- c(
    "All",
    "Up",
    "Down"
  )
  
  for (gene_direction in gene_directions) {
    
    selected_degs <- switch(
      gene_direction,
      
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
    
    fold_change <- make_fold_change_vector(
      selected_degs
    )
    
    cat(
      "\n",
      gene_direction,
      "fold-change genes:",
      length(fold_change),
      "\n"
    )
    
    if (length(fold_change) == 0) {
      next
    }
    
    #-------------------------------------------------------
    # GO Biological Process
    #-------------------------------------------------------
    
    go_file <- file.path(
      go_directory,
      paste0(
        comparison_name,
        "_",
        gene_direction,
        "_GO_BP.rds"
      )
    )
    
    if (file.exists(go_file)) {
      
      go_result <- readRDS(
        go_file
      )
      
      save_cnetplot(
        enrichment_object = go_result,
        fold_change = fold_change,
        plot_title = paste(
          comparison_name,
          gene_direction,
          "GO Biological Process"
        ),
        output_prefix = paste0(
          comparison_name,
          "_",
          gene_direction,
          "_GO_BP"
        )
      )
      
    } else {
      
      cat(
        "GO result not found:",
        go_file,
        "\n"
      )
    }
    
    #-------------------------------------------------------
    # Reactome
    #-------------------------------------------------------
    
    reactome_file <- file.path(
      reactome_directory,
      paste0(
        comparison_name,
        "_",
        gene_direction,
        "_Reactome.rds"
      )
    )
    
    if (file.exists(reactome_file)) {
      
      reactome_result <- readRDS(
        reactome_file
      )
      
      save_cnetplot(
        enrichment_object = reactome_result,
        fold_change = fold_change,
        plot_title = paste(
          comparison_name,
          gene_direction,
          "Reactome"
        ),
        output_prefix = paste0(
          comparison_name,
          "_",
          gene_direction,
          "_Reactome"
        )
      )
      
    } else {
      
      cat(
        "Reactome result not found:",
        reactome_file,
        "\n"
      )
    }
  }
}

#-----------------------------------------------------------
# 6. Run all comparisons
#-----------------------------------------------------------

invisible(
  lapply(
    deg_files,
    process_comparison
  )
)

cat("\nScript 08 complete.\n")