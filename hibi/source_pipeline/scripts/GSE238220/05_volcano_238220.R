############################################################
# GSE238220
# Volcano plots for DESeq2 results
#
# Positive log2FoldChange = higher in HI
# Negative log2FoldChange = higher in Control
############################################################

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

#-----------------------------------------------------------
# 1. Paths and thresholds
#-----------------------------------------------------------

input_directory <- "results/GSE238220/DE"
figure_directory <- "figures/GSE238220/volcano"

fdr_threshold <- 0.05
logfc_threshold <- 0.2
number_of_labels <- 12

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

#-----------------------------------------------------------
# 2. Find all complete DESeq2 result files
#-----------------------------------------------------------

result_files <- list.files(
  input_directory,
  pattern = "_all_genes\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0) {
  stop(
    "No DESeq2 result files were found in:\n",
    input_directory,
    "\nRun Script 04 first."
  )
}

cat("\nResult files found:\n")
print(result_files)

#-----------------------------------------------------------
# 3. Function to make one volcano plot
#-----------------------------------------------------------

make_volcano <- function(result_file) {
  
  comparison_name <- sub(
    "_all_genes\\.csv$",
    "",
    basename(result_file)
  )
  
  cat("\n============================================\n")
  cat("Creating volcano plot:", comparison_name, "\n")
  cat("============================================\n")
  
  result <- read_csv(
    result_file,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "input_id",
    "SYMBOL",
    "log2FoldChange",
    "padj"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(result)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      comparison_name,
      " is missing these columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  #---------------------------------------------------------
  # Prepare plotting variables
  #---------------------------------------------------------
  
  result <- result %>%
    mutate(
      plot_status = case_when(
        !is.na(padj) &
          padj < fdr_threshold &
          log2FoldChange >= logfc_threshold ~ "Up in HI",
        
        !is.na(padj) &
          padj < fdr_threshold &
          log2FoldChange <= -logfc_threshold ~ "Down in HI",
        
        TRUE ~ "Not significant"
      ),
      
      minus_log10_fdr = case_when(
        is.na(padj) ~ NA_real_,
        padj == 0 ~ -log10(.Machine$double.xmin),
        TRUE ~ -log10(padj)
      ),
      
      gene_label = case_when(
        !is.na(SYMBOL) &
          SYMBOL != "" ~ SYMBOL,
        
        TRUE ~ input_id
      )
    )
  
  # Label the strongest significant genes
  label_data <- result %>%
    filter(
      plot_status != "Not significant",
      !is.na(minus_log10_fdr),
      !is.na(log2FoldChange)
    ) %>%
    arrange(
      padj,
      desc(abs(log2FoldChange))
    ) %>%
    slice_head(
      n = number_of_labels
    )
  
  up_count <- sum(
    result$plot_status == "Up in HI",
    na.rm = TRUE
  )
  
  down_count <- sum(
    result$plot_status == "Down in HI",
    na.rm = TRUE
  )
  
  subtitle_text <- paste0(
    "Up in HI: ",
    up_count,
    " | Down in HI: ",
    down_count,
    " | FDR < ",
    fdr_threshold,
    " and |log2FC| ≥ ",
    logfc_threshold
  )
  
  #---------------------------------------------------------
  # Build volcano plot
  #---------------------------------------------------------
  
  volcano_plot <- ggplot(
    result,
    aes(
      x = log2FoldChange,
      y = minus_log10_fdr
    )
  ) +
    geom_point(
      aes(shape = plot_status),
      alpha = 0.65,
      size = 1.8,
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = c(
        -logfc_threshold,
        logfc_threshold
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = gene_label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.25,
      min.segment.length = 0,
      na.rm = TRUE
    ) +
    labs(
      title = comparison_name,
      subtitle = subtitle_text,
      x = "log2 fold change: HI vs Control",
      y = expression(-log[10]("adjusted P value")),
      shape = "Expression status"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )
  
  #---------------------------------------------------------
  # Save PDF and PNG
  #---------------------------------------------------------
  
  pdf_file <- file.path(
    figure_directory,
    paste0(
      comparison_name,
      "_volcano.pdf"
    )
  )
  
  png_file <- file.path(
    figure_directory,
    paste0(
      comparison_name,
      "_volcano.png"
    )
  )
  
  ggsave(
    filename = pdf_file,
    plot = volcano_plot,
    width = 8,
    height = 6
  )
  
  ggsave(
    filename = png_file,
    plot = volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  cat("Saved:\n")
  cat(pdf_file, "\n")
  cat(png_file, "\n")
  
  #---------------------------------------------------------
  # Save plotted significant genes
  #---------------------------------------------------------
  
  write_csv(
    label_data %>%
      select(
        input_id,
        SYMBOL,
        log2FoldChange,
        padj,
        plot_status
      ),
    
    file.path(
      figure_directory,
      paste0(
        comparison_name,
        "_volcano_labeled_genes.csv"
      )
    )
  )
  
  invisible(volcano_plot)
}

#-----------------------------------------------------------
# 4. Generate all volcano plots
#-----------------------------------------------------------

invisible(
  lapply(
    result_files,
    make_volcano
  )
)

cat("\nScript 05 complete.\n")