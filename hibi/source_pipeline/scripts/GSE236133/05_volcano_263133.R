############################################################
# GSE236133
# Volcano plots for all differential-expression comparisons
############################################################

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

#-----------------------------------------------------------
# 1. Paths and thresholds
#-----------------------------------------------------------

results_directory <- "results/GSE236133/DE"
figure_directory <- "figures/GSE236133/volcano"

logfc_threshold <- 0.2
fdr_threshold <- 0.05

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

#-----------------------------------------------------------
# 2. Locate full DEG result files
#-----------------------------------------------------------

result_files <- list.files(
  results_directory,
  pattern = "_all_genes\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0) {
  stop(
    "No full differential-expression result files found.",
    "\nRun Script 04 first."
  )
}

#-----------------------------------------------------------
# 3. Create one volcano plot
#-----------------------------------------------------------

make_volcano <- function(result_file) {
  
  result <- read_csv(
    result_file,
    show_col_types = FALSE
  )
  
  comparison_name <- unique(result$comparison)
  
  if (length(comparison_name) != 1) {
    comparison_name <- sub(
      "_all_genes\\.csv$",
      "",
      basename(result_file)
    )
  }
  
  result <- result %>%
    mutate(
      plot_status = case_when(
        adj.P.Val < fdr_threshold &
          logFC >= logfc_threshold ~ "Up",
        
        adj.P.Val < fdr_threshold &
          logFC <= -logfc_threshold ~ "Down",
        
        TRUE ~ "Not significant"
      ),
      
      minus_log10_fdr = -log10(
        pmax(adj.P.Val, .Machine$double.xmin)
      ),
      
      label = case_when(
        plot_status != "Not significant" &
          !is.na(SYMBOL) &
          SYMBOL != "" ~ SYMBOL,
        
        plot_status != "Not significant" ~ ENSEMBL,
        
        TRUE ~ NA_character_
      )
    )
  
  label_data <- result %>%
    filter(plot_status != "Not significant") %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    slice_head(n = 12)
  
  subtitle_text <- paste0(
    "Up: ",
    sum(result$plot_status == "Up"),
    " | Down: ",
    sum(result$plot_status == "Down"),
    " | FDR < ",
    fdr_threshold,
    ", |logFC| ≥ ",
    logfc_threshold
  )
  
  volcano_plot <- ggplot(
    result,
    aes(
      x = logFC,
      y = minus_log10_fdr
    )
  ) +
    geom_point(
      aes(shape = plot_status),
      alpha = 0.65,
      size = 1.7
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
      aes(label = label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    labs(
      title = comparison_name,
      subtitle = subtitle_text,
      x = "log2 fold change: ipsilateral vs contralateral",
      y = expression(-log[10]("adjusted P value")),
      shape = "Status"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold"
      )
    )
  
  output_name <- file.path(
    figure_directory,
    paste0(
      comparison_name,
      "_volcano.pdf"
    )
  )
  
  ggsave(
    filename = output_name,
    plot = volcano_plot,
    width = 8,
    height = 6
  )
  
  ggsave(
    filename = sub(
      "\\.pdf$",
      ".png",
      output_name
    ),
    plot = volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  cat("Saved:", output_name, "\n")
}

#-----------------------------------------------------------
# 4. Generate all plots
#-----------------------------------------------------------

invisible(
  lapply(
    result_files,
    make_volcano
  )
)

cat("\nScript 05 complete.\n")