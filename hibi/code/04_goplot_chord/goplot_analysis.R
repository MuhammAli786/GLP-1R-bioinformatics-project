#!/usr/bin/env Rscript
# GOplot visualisations of the HIBI (hypoxic-ischemic brain injury) consensus data;
# a port of the GLP-1R 04_goplot_chord/goplot_analysis.R with the same calls,
# parameters and plot set.
# Inputs: hibi_data/goplot/<name>_terms.csv and _genes.csv
# Outputs: hibi_goplot/<Name>/<PNG|PDF>/
#
# Plots produced: GOBubble (top 20 GO terms per category for consensus, top 5 for
# pathways), GOCircle overview for consensus, gene-to-term GOChord for consensus
# and restricted pathways across GO, KEGG, Reactome and AllEnrichments, and a
# genes x terms heatmap. GOplot::GOHeat's own x-axis labels break under modern
# ggplot2 (it calls scale_x_discrete on a numeric axis), so that heatmap is
# reproduced by the go_heat helper.
suppressMessages(library(GOplot))
suppressMessages(library(ggplot2))

BASE <- Sys.getenv("GLP1R_BASE")
if (!nzchar(BASE)) stop("Set GLP1R_BASE to the directory holding the analysis data tree, e.g. export GLP1R_BASE=/path/to/workspace")
DATA  <- file.path(BASE, "mnt", "outputs", "hibi_data", "goplot")
PLOTS <- file.path(BASE, "mnt", "outputs", "hibi_goplot")

save_plot <- function(fn, name, fname, w = 11, h = 9) {
  for (d in list(c("PNG", "png"), c("PDF", "pdf"))) {
    dir.create(file.path(PLOTS, name, d[1]), recursive = TRUE, showWarnings = FALSE)
    f <- file.path(PLOTS, name, d[1], paste0(fname, ".", d[2]))
    if (d[2] == "png") png(f, width = w, height = h, units = "in", res = 200, bg = "white") else pdf(f, width = w, height = h)
    p <- tryCatch(fn(), error = function(e) { message("  ERR ", fname, ": ", conditionMessage(e)); NULL })
    if (inherits(p, "ggplot")) print(p) else if (inherits(p, c("grob", "gtable"))) grid::grid.draw(p)
    dev.off()
  }
  cat("  saved", name, "/", fname, "\n")
}

load_circ <- function(name) {
  terms <- read.csv(file.path(DATA, paste0(name, "_terms.csv")), stringsAsFactors = FALSE)
  genes <- read.csv(file.path(DATA, paste0(name, "_genes.csv")), stringsAsFactors = FALSE)
  circle_dat(terms, genes)
}
top_per_cat <- function(circ, n = 5) {
  u <- circ[!duplicated(circ$term), ]
  sel <- do.call(rbind, lapply(split(u, u$category), function(d) head(d[order(d$adj_pval), ], n)))
  circ[circ$term %in% sel$term, ]
}
top_terms <- function(circ, n = 6) {
  u <- circ[!duplicated(circ$term), ]
  head(u[order(u$adj_pval), "term"], n)
}

build_mem <- function(circ, terms, genes_vec = NULL) {
  if (is.null(genes_vec)) genes_vec <- unique(circ$genes[circ$term %in% terms])
  m <- sapply(terms, function(t) as.integer(genes_vec %in% circ$genes[circ$term == t]))
  if (is.null(dim(m))) m <- matrix(m, nrow = length(genes_vec))
  rownames(m) <- genes_vec; colnames(m) <- terms
  m[rowSums(m) > 0, , drop = FALSE]
}

go_heat <- function(mem, mode = "count", lfc = NULL, title = "") {
  cnt <- rowSums(mem)
  ord <- if (nrow(mem) > 2) hclust(dist(mem))$order else seq_len(nrow(mem))
  genes <- rownames(mem)[ord]; terms <- colnames(mem)
  long <- expand.grid(gene = genes, term = terms, stringsAsFactors = FALSE)
  long$member <- mapply(function(g, t) mem[g, t], long$gene, long$term)
  if (mode == "count") {
    long$val <- ifelse(long$member == 1, cnt[long$gene], 0)
    fill <- scale_fill_gradientn(colours = c("#2E80E0", "#8E6FB0", "#B0182C"), name = "Count")
  } else {
    cap <- 2
    long$val <- ifelse(long$member == 1, lfc[long$gene], 0)
    long$val <- pmax(pmin(long$val, cap), -cap)
    fill <- scale_fill_gradient2(low = "green", mid = "yellow", high = "red",
                                 midpoint = 0, name = "logFC", limits = c(-cap, cap))
  }
  long$gene <- factor(long$gene, levels = genes)
  long$term <- factor(long$term, levels = rev(terms))
  ggplot(long, aes(gene, term, fill = val)) +
    geom_tile() + fill +
    labs(title = title) + coord_cartesian(expand = FALSE) + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.text.y = element_text(size = 12), axis.title = element_blank(),
          plot.title = element_text(size = 13, face = "bold"), panel.grid = element_blank())
}

make_cat_chord <- function(prefix, folder, suffix, label, n_terms = 6) {
  tf <- file.path(DATA, paste0(prefix, "_", suffix, "_terms.csv"))
  if (!file.exists(tf)) return(invisible())
  terms <- read.csv(tf, stringsAsFactors = FALSE)
  if (nrow(terms) < 2) { cat("  skip", folder, label, "(<2 terms)\n"); return(invisible()) }
  genes <- read.csv(file.path(DATA, paste0(prefix, "_genes.csv")), stringsAsFactors = FALSE)
  genes$ID <- toupper(genes$ID)
  circ2 <- circle_dat(terms, genes)
  proc <- top_terms(circ2, n_terms)
  if (prefix == "consensus") {
    gv <- unique(unlist(lapply(proc, function(t) {
      s <- circ2[circ2$term == t, ]; head(s[order(-abs(s$logFC)), "genes"], 5) })))
    gdf <- genes[genes$ID %in% gv, ]
  } else gdf <- genes
  save_plot(function() {
    ch <- chord_dat(data = circ2, genes = gdf, process = proc)
    GOChord(ch, space = 0.02, gene.order = "logFC", gene.space = 0.25,
            gene.size = 4, process.label = 8)
  }, folder, paste0("GOChord_", folder, "_", label), 12, 11)
}

make_combined_chord <- function(prefix, folder, per_cat = 2) {
  files <- file.path(DATA, paste0(prefix, c("_terms.csv", "_kegg_terms.csv", "_reactome_terms.csv")))
  files <- files[file.exists(files)]
  terms <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  genes <- read.csv(file.path(DATA, paste0(prefix, "_genes.csv")), stringsAsFactors = FALSE)
  genes$ID <- toupper(genes$ID)
  circ2 <- circle_dat(terms, genes)
  u <- circ2[!duplicated(circ2$term), ]
  sel <- do.call(rbind, lapply(split(u, u$category), function(d) head(d[order(d$adj_pval), ], per_cat)))
  proc <- sel$term
  if (prefix == "consensus") {
    gv <- unique(unlist(lapply(proc, function(t) {
      s <- circ2[circ2$term == t, ]; head(s[order(-abs(s$logFC)), "genes"], 5) })))
    gdf <- genes[genes$ID %in% gv, ]
  } else gdf <- genes
  save_plot(function() {
    ch <- chord_dat(data = circ2, genes = gdf, process = proc)
    GOChord(ch, space = 0.02, gene.order = "logFC", gene.space = 0.25, gene.size = 4, process.label = 7)
  }, folder, paste0("GOChord_", folder, "_AllEnrichments"), 14, 13)
}

# Consensus set
cat("Consensus...\n")
circ <- load_circ("consensus")
save_plot(function() GOBubble(top_per_cat(circ, 20), display = "multiple", bg.col = TRUE,
          labels = 6, table.legend = FALSE, ID = FALSE,
          title = "HIBI consensus DEGs - top 20 enriched GO terms per category"),
          "Consensus", "GOBubble_consensus_top20percat", 15, 10)
save_plot(function() GOCircle(circ, nsub = top_terms(circ, 10),
          title = "HIBI consensus DEGs - GO enrichment overview"),
          "Consensus", "GOCircle_consensus_top10", 12, 10)
heat_consensus <- function() {
  cg <- read.csv(file.path(DATA, "consensus_genes.csv"), stringsAsFactors = FALSE)
  lfc <- setNames(cg$logFC, toupper(cg$ID))
  rc <- reduce_overlap(circ, overlap = 0.75)
  tt <- top_terms(rc, 7)
  rep_genes <- unique(unlist(lapply(tt, function(t) {
    s <- circ[circ$term == t, ]; head(s[order(-abs(s$logFC)), "genes"], 6) })))
  mem <- build_mem(circ, tt, rep_genes)
  go_heat(mem, mode = "logfc", lfc = lfc, title = "HIBI consensus DEGs - genes x GO terms (logFC)")
}
save_plot(heat_consensus, "Consensus", "GOHeat_consensus", 15, 8)
# Consensus GO chord, capping genes to the top 5 per term by |logFC|
save_plot(function() {
  proc <- top_terms(circ, 6)
  gv <- unique(unlist(lapply(proc, function(t) {
    s <- circ[circ$term == t, ]; head(s[order(-abs(s$logFC)), "genes"], 5) })))
  cg <- read.csv(file.path(DATA, "consensus_genes.csv"), stringsAsFactors = FALSE)
  cg$ID <- toupper(cg$ID)
  ch <- chord_dat(data = circ, genes = cg[cg$ID %in% gv, ], process = proc)
  GOChord(ch, space = 0.02, gene.order = "logFC", gene.space = 0.25, gene.size = 4, process.label = 8)
}, "Consensus", "GOChord_Consensus", 12, 11)

# Restricted pathway gene sets
for (name in c("BBB", "Inflammatory", "Survival", "IonChannel")) {
  if (!file.exists(file.path(DATA, paste0(name, "_terms.csv")))) { cat(name, "skipped\n"); next }
  cat(name, "...\n")
  cc <- load_circ(name)
  genes <- read.csv(file.path(DATA, paste0(name, "_genes.csv")), stringsAsFactors = FALSE)
  genes$ID <- toupper(genes$ID)
  lfc <- setNames(genes$logFC, genes$ID)
  save_plot(function() GOBubble(top_per_cat(cc, 5), display = "multiple", bg.col = TRUE,
            labels = 0.5, table.legend = FALSE, ID = FALSE,
            title = paste0("HIBI ", name, " - top 5 enriched GO terms per category")),
            name, paste0("GOBubble_", name, "_top5percat"), 13, 9)
  proc <- top_terms(cc, 7)
  save_plot(function() {
    ch <- chord_dat(data = cc, genes = genes, process = proc)
    GOChord(ch, space = 0.02, gene.order = "logFC", gene.space = 0.25, gene.size = 4, process.label = 8)
  }, name, paste0("GOChord_", name), 12, 11)
  save_plot(function() {
    mem <- build_mem(cc, proc, intersect(genes$ID, unique(cc$genes)))
    go_heat(mem, mode = "logfc", lfc = lfc, title = paste0("HIBI ", name, " - genes x GO terms (logFC)"))
  }, name, paste0("GOHeat_", name), 12, 8)
}

# KEGG and Reactome chords
cat("KEGG / Reactome chords...\n")
targets <- list(c("consensus", "Consensus"), c("BBB", "BBB"),
                c("Inflammatory", "Inflammatory"), c("Survival", "Survival"),
                c("IonChannel", "IonChannel"))
for (tg in targets) {
  make_cat_chord(tg[1], tg[2], "kegg", "KEGG")
  make_cat_chord(tg[1], tg[2], "reactome", "Reactome")
}
for (tg in targets) make_combined_chord(tg[1], tg[2], per_cat = 2)
cat("DONE\n")
