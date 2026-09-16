#!/usr/bin/env Rscript
# Figure 3.7: Body weight and blood glucose at 24 h post-injection.
# (A) Percent change in body weight, pre- to post-surgery, for Sham (n = 4),
#     HI-Vehicle (n = 19), HI-Semaglutide (n = 19). One-way ANOVA, with
#     baseline-adjusted confirmation by ANCOVA (covariate = pre-surgery weight)
#     and a litter-random-effect mixed model.
# (B) Blood glucose at the endpoint for the same three groups. One-way ANOVA.
# Data are mean +/- SEM; individual animals overlaid as points.
# Styled with ggprism to match GraphPad Prism's look (theme, palette, and
# significance brackets).
#
# Input:  ../../data/thesis_figures/figure_3_7_bodyweight_glucose.csv
# Output: ../../figures/thesis/Figure_3_7_bodyweight_glucose.{png,pdf}
#
# Litter is coded "First".."Fifth" per group and is shared across the
# HI-Vehicle and HI-Semaglutide arms (split-litter design: littermates
# assigned across treatment groups) -- so litter is used as-is as the
# mixed-model grouping factor, not nested within group.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(ggprism)
  library(patchwork)
  library(lme4)
  library(lmerTest)
})

# ---- paths (edit if the project moves) ----
PROJECT_ROOT <- "/Users/muhammadali/Documents/GLP-1R-bioinformatics-project"
DATA_CSV     <- file.path(PROJECT_ROOT, "data/thesis_figures/figure_3_7_bodyweight_glucose.csv")
OUT_DIR      <- file.path(PROJECT_ROOT, "figures/thesis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- font sizes -- bumped up throughout; change here to resize everything ----
FS_BASE        <- 20    # theme_prism base size (drives axis text/ticks)
FS_AXIS_TITLE  <- 24
FS_SIG_TEXT    <- 7.5   # ggprism::add_pvalue label size
FS_PANEL_LABEL <- 28    # "A" / "B" panel tags
FS_N_LABEL     <- 6     # "n = ..." label inside each bar

GROUP_LEVELS <- c("Sham", "HI-Vehicle", "HI-Semaglutide")

# ---- load + derive ----
df <- read.csv(DATA_CSV, stringsAsFactors = FALSE)
df$group  <- factor(df$group, levels = GROUP_LEVELS)
df$litter <- factor(df$litter, levels = c("First", "Second", "Third", "Fourth", "Fifth"))
df$pct_weight_change <- (df$weight_post - df$weight_pre) / df$weight_pre * 100

# ================= Panel A stats: body weight percent change =================
aov_weight <- aov(pct_weight_change ~ group, data = df)
cat("=== Body weight %% change: one-way ANOVA ===\n"); print(summary(aov_weight))

tukey_weight <- TukeyHSD(aov_weight)
cat("\n=== Body weight %% change: Tukey post-hoc ===\n"); print(tukey_weight)

ancova_weight <- aov(pct_weight_change ~ weight_pre + group, data = df)
cat("\n=== Body weight %% change: ANCOVA, baseline (pre-surgery weight) adjusted ===\n")
print(summary(ancova_weight))

mixed_weight <- lmer(pct_weight_change ~ group + (1 | litter), data = df)
cat("\n=== Body weight %% change: litter-random-effect mixed model ===\n")
print(anova(mixed_weight))

# ================= Panel B stats: endpoint blood glucose =================
aov_glucose <- aov(glucose ~ group, data = df)
cat("\n=== Blood glucose: one-way ANOVA ===\n"); print(summary(aov_glucose))

tukey_glucose <- TukeyHSD(aov_glucose)
cat("\n=== Blood glucose: Tukey post-hoc ===\n"); print(tukey_glucose)

# ---- helper: significance stars from an adjusted p-value ----
stars <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

# Build a Prism-style bracket table (group1, group2, label, y.position) from
# a TukeyHSD result, for ggprism::add_pvalue.
bracket_table <- function(tukey_tbl, y_top, step) {
  tt <- as.data.frame(tukey_tbl)
  tibble::tibble(
    group1 = c("Sham", "Sham", "HI-Vehicle"),
    group2 = c("HI-Vehicle", "HI-Semaglutide", "HI-Semaglutide"),
    label = c(
      stars(tt["HI-Vehicle-Sham", "p adj"]),
      stars(tt["HI-Semaglutide-Sham", "p adj"]),
      stars(tt["HI-Semaglutide-HI-Vehicle", "p adj"])
    ),
    y.position = y_top + c(0, step, step * 2)
  )
}

# ---- mean +/- SEM per group ----
summarise_group <- function(data, value_col) {
  data %>%
    group_by(group) %>%
    summarise(mean = mean(.data[[value_col]]),
              sem  = sd(.data[[value_col]]) / sqrt(dplyr::n()),
              n    = dplyr::n(),
              .groups = "drop")
}

weight_summary  <- summarise_group(df, "pct_weight_change")
glucose_summary <- summarise_group(df, "glucose")

# n-label y position: just inside each bar, near the baseline (bottom of the
# bar/chart) rather than mid-bar. Weight bars can be negative (semaglutide
# loses weight), so the offset is signed to stay inside the bar either way.
N_LABEL_OFFSET_WEIGHT  <- 3
N_LABEL_OFFSET_GLUCOSE <- 0.4
weight_summary$n_label_y  <- sign(weight_summary$mean) * N_LABEL_OFFSET_WEIGHT
glucose_summary$n_label_y <- N_LABEL_OFFSET_GLUCOSE

prism_theme <- theme_prism(base_size = FS_BASE) +
  theme(
    axis.title      = element_text(size = FS_AXIS_TITLE, face = "bold"),
    legend.position = "none",
    plot.title      = element_text(size = FS_PANEL_LABEL, face = "bold", hjust = -0.15)
  )

# ================= Panel A plot =================
weight_brackets <- bracket_table(
  tukey_weight$group,
  y_top = max(df$pct_weight_change) + 4,
  step  = 6
)

panel_a <- ggplot(df, aes(x = group, y = pct_weight_change)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.6) +
  geom_bar(data = weight_summary, aes(x = group, y = mean, fill = group),
           stat = "identity", width = 0.6, colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  geom_errorbar(data = weight_summary, aes(x = group, ymin = mean - sem, ymax = mean + sem),
                width = 0.2, linewidth = 0.8, inherit.aes = FALSE) +
  geom_jitter(width = 0.12, size = 2.5, shape = 21, colour = "black", fill = "white", stroke = 0.6, alpha = 0.9) +
  geom_text(data = weight_summary, aes(x = group, y = n_label_y, label = paste0("n = ", n)),
            inherit.aes = FALSE, size = FS_N_LABEL, fontface = "bold", colour = "black") +
  scale_fill_prism(palette = "colors") +
  add_pvalue(weight_brackets, label = "label", tip.length = 0.01,
             label.size = FS_SIG_TEXT, bracket.size = 0.9) +
  scale_y_continuous(limits = c(-20, 60), breaks = seq(-20, 60, by = 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Body weight change (%)", title = "A") +
  prism_theme

# ================= Panel B plot =================
glucose_brackets <- bracket_table(
  tukey_glucose$group,
  y_top = max(df$glucose) + 0.6,
  step  = 0.9
)

panel_b <- ggplot(df, aes(x = group, y = glucose)) +
  geom_bar(data = glucose_summary, aes(x = group, y = mean, fill = group),
           stat = "identity", width = 0.6, colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  geom_errorbar(data = glucose_summary, aes(x = group, ymin = mean - sem, ymax = mean + sem),
                width = 0.2, linewidth = 0.8, inherit.aes = FALSE) +
  geom_jitter(width = 0.12, size = 2.5, shape = 21, colour = "black", fill = "white", stroke = 0.6, alpha = 0.9) +
  geom_text(data = glucose_summary, aes(x = group, y = n_label_y, label = paste0("n = ", n)),
            inherit.aes = FALSE, size = FS_N_LABEL, fontface = "bold", colour = "black") +
  scale_fill_prism(palette = "colors") +
  add_pvalue(glucose_brackets, label = "label", tip.length = 0.01,
             label.size = FS_SIG_TEXT, bracket.size = 0.9) +
  labs(x = NULL, y = "Blood glucose", title = "B") +
  prism_theme

combined <- panel_a + panel_b

ggsave(file.path(OUT_DIR, "Figure_3_7_bodyweight_glucose.png"), combined, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Figure_3_7_bodyweight_glucose.pdf"), combined, width = 14, height = 7)

cat("\nSaved Figure 3.7 to", OUT_DIR, "\n")
