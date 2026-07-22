================================================================
 GOPLOT ANALYSIS
 GLP-1R agonist CNS meta-analysis — GOplot visualizations
================================================================

Reference: GOplot package, https://wencke.github.io/
  Walter, Sanchez-Cabo & Ricote (2015), Bioinformatics, btv300.
GOplot is a VISUALIZATION package: it combines expression (logFC)
with functional-enrichment results. It does not perform the
enrichment itself — that is done upstream (Enrichr/gseapy) and fed in.

INPUT FORMAT (circle_dat)
-------------------------
circle_dat(terms, genes) needs two tables (built by
scripts/prepare_goplot_data.py into data/):
  <name>_terms.csv : Category(BP/MF/CC), ID(GO:...), Term, Genes, adj_pval
  <name>_genes.csv : ID (gene symbol), logFC
It returns the "circ" object (category, ID, term, count, gene, logFC,
adj_pval, zscore) used by every plotting function.

DATA SOURCE
-----------
All inputs come from the consensus overlapping DEGs (>=2 groups) at
the LFC0.2 threshold:
  * Consensus  : the existing consensus GO enrichment (BP/MF/CC) +
                 consensus gene mean log2FC.
  * Restricted : for BBB, Inflammatory (JAK-STAT3), Survival (PI3K/Akt)
                 and IonChannel, the consensus genes are intersected
                 with the curated pathway gene list (cnet_gene_lists.py)
                 and re-enriched against GO (BP/MF/CC).

PLOTS (plots/<Name>/<PNG|PDF>/)
-------------------------------
  Consensus/
    GOBubble_consensus_top20percat  z-score vs -log10(padj); bubble area =
                                    gene count; faceted BP/MF/CC; top 20 per category
    GOCircle_consensus_top10        circular overview of the top 10 terms
                                    (outer ring = per-gene logFC scatter)
    GOHeat_consensus                genes x GO terms; GOHeat nlfc=1 (logFC) style:
                                    tiles = gene logFC, red/yellow/green, no borders
  BBB/ Inflammatory/ Survival/ IonChannel/
    GOBubble_<name>_top5percat      top 5 enriched GO terms per category
    GOChord_<name>                  gene<->term chord (GO terms); ribbons =
                                    membership, gene squares coloured by logFC
    GOHeat_<name>                   genes x terms; logFC (red/yellow/green), no borders

  KEGG & Reactome chords (GOplot is database-agnostic, not GO-only):
    GOChord_<name>_KEGG             gene<->pathway chord using KEGG terms
    GOChord_<name>_Reactome         gene<->pathway chord using Reactome terms
    (produced for Consensus + BBB/Inflammatory/Survival/IonChannel; the
     consensus chords cap genes to the top 5 per term by |logFC|. Note the
     chord legend is titled "GO Terms" by GOplot even for KEGG/Reactome.)
    GOChord_<name>_AllEnrichments   ALL databases in one chord: top 2 terms
                                    from each of GO BP/MF/CC + KEGG + Reactome
                                    (10 terms) linked to the genes by logFC.

  BBB presence (>=1 group) chords  (scripts/prepare_bbb_freq1.py):
    GOChord_BBB_Freq1               GO chord using BBB genes in >=1 group (36
    GOChord_BBB_KEGG_Freq1          genes vs 26 consensus); re-enriched on the
    GOChord_BBB_Reactome_Freq1      wider set, so it surfaces the MMPs and the
    GOChord_BBB_AllEnrichments_Freq1   metallopeptidase / focal-adhesion terms.

HOW TO REPRODUCE
----------------
  # 1. build the GOplot input tables (needs pandas + gseapy)
  python scripts/prepare_goplot_data.py
  # 2. render all plots (needs R + GOplot)
  Rscript scripts/goplot_analysis.R

ENVIRONMENT
-----------
  R 4.5.3 with GOplot 1.0.2 (+ ggplot2, ggdendro, gridExtra, RColorBrewer),
  installed into a user-space micromamba env. PNG (200 dpi, white
  background so labels are legible) + vector PDF for every figure.

NOTES
-----
  * Gene symbols are upper-cased for GOplot's internal matching (so the
    plots show e.g. CAV1 rather than Cav1).
  * GOBubble, GOChord and GOCircle use the real GOplot functions. GOplot's
    own GOHeat drops its x-axis gene labels under modern ggplot2 (it calls
    scale_x_discrete on a numeric axis), so the GOHeat panels are
    reproduced with a small helper (go_heat) that renders the identical
    GOHeat nlfc=1 logFC style WITH gene labels: tiles coloured by gene
    logFC (red high / yellow 0 / green low), unassigned tiles yellow, NO
    borders between squares; values clamped to +/- 2.
  * The consensus GOHeat shows representative genes (top 6 by |logFC| per
    term) across 7 non-redundant (reduce_overlap) GO terms.
