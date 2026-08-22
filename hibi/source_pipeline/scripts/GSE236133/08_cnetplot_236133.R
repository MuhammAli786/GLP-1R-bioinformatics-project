# GSE236133 GO BP and Reactome cnetplots; input: results/GSE236133/DE/*_DEGs.csv, enrichment RDS -> output: figures/GSE236133/cnetplots/

library(enrichplot)
library(clusterProfiler)
library(ReactomePA)
library(readr)
library(dplyr)
library(ggplot2)

# Paths
de_directory <- "results/GSE236133/DE"

go_directory <- "results/GSE236133/enrichment/GO_KEGG"

reactome_directory <- "results/GSE236133/enrichment/Reactome"

figure_directory <- "figures/GSE236133/cnetplots"

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# Find DEG files
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

# Create named fold-change vector
make_fold_change_vector <- function(deg) {
  
  fold_change_table <- deg %>%
    filter(
      !is.na(ENTREZID),
      ENTREZID != "",
      !is.na(logFC)
    ) %>%
    arrange(adj.P.Val) %>%
    distinct(
      ENTREZID,
      .keep_all = TRUE
    )
  
  fold_change <- fold_change_table$logFC
  
  names(fold_change) <- as.character(
    fold_change_table$ENTREZID
  )
  
  fold_change
}

# Save one cnetplot safely with fallback layout
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
      "Skipping empty enrichment result:",
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
  
  network_plot <- tryCatch(
    
    enrichplot::cnetplot(
      enrichment_object,
      showCategory = categories_to_show,
      foldChange = fold_change,
      layout = "kk"
    ) +
      ggtitle(plot_title),
    
    error = function(e) {
      
      cat(
        "The requested layout failed. ",
        "Retrying with the default layout.\n"
      )
      
      enrichplot::cnetplot(
        enrichment_object,
        showCategory = categories_to_show,
        foldChange = fold_change
      ) +
        ggtitle(plot_title)
    }
  )
  
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
    plot = network_plot,
    width = 12,
    height = 10
  )
  
  ggsave(
    filename = png_file,
    plot = network_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  cat(
    "Saved:\n",
    pdf_file,
    "\n",
    png_file,
    "\n"
  )
  
  invisible(network_plot)
}

# Process one comparison
process_comparison <- function(deg_file) {
  
  comparison_name <- sub(
    "_DEGs\\.csv$",
    "",
    basename(deg_file)
  )
  
  deg <- read_csv(
    deg_file,
    show_col_types = FALSE
  )
  
  gene_directions <- c(
    "All",
    "Up",
    "Down"
  )
  
  for (gene_direction in gene_directions) {
    
    selected_degs <- switch(
      gene_direction,
      All = deg,
      Up = deg %>% filter(direction == "Up"),
      Down = deg %>% filter(direction == "Down")
    )
    
    fold_change <- make_fold_change_vector(
      selected_degs
    )
    
    if (length(fold_change) == 0) {
      
      cat(
        "No fold-change genes for",
        comparison_name,
        gene_direction,
        "\n"
      )
      
      next
    }
    
    # GO Biological Process
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
      
      go_result <- readRDS(go_file)
      
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
    
    # Reactome
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

# Run all comparisons
invisible(
  lapply(
    deg_files,
    process_comparison
  )
)

cat("\nScript 08 complete.\n")