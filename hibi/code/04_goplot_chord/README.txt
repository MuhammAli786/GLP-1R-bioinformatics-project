================================================================
 GOPLOT ANALYSIS — HIBI
 Neonatal hypoxic-ischemic brain injury meta-analysis
 (port of GLP-1R-bioinformatics-project/code/04_goplot_chord)
================================================================

Reference: GOplot package, https://wencke.github.io/
  Walter, Sanchez-Cabo & Ricote (2015), Bioinformatics, btv300.
GOplot is a VISUALIZATION package: it combines expression (logFC)
with functional-enrichment results. It does not perform the
enrichment itself — that is done upstream (Enrichr/gseapy) and fed in.

INPUT FORMAT (circle_dat)
-------------------------
circle_dat(terms, genes) needs two tables (built by
prepare_goplot_data.py into ../../data/goplot/):
  <name>_terms.csv : Category(BP/MF/CC), ID(GO:...), Term, Genes, adj_pval
  <name>_genes.csv : ID (gene symbol), logFC
It returns the "circ" object (category, ID, term, count, gene, logFC,
adj_pval, zscore) used by every plotting function.

DATA SOURCE
-----------
All inputs come from the HIBI consensus overlapping DEGs (>=2 groups)
at the LFC0.2 threshold (12 groups, 6 GEO datasets):
  * Consensus  : the top 1500 consensus genes by recurrence (n_groups,
                 tie-broken by |mean log2FC|), enriched against
                 GO BP/MF/CC + KEGG + Reactome.
  * Restricted : for BBB, Inflammatory (JAK-STAT3), Survival (PI3K/Akt)
                 and IonChannel, the consensus genes are intersected
                 with the curated pathway gene list (cnet_gene_lists.py)
                 and re-enriched. Gene counts present in HIBI consensus:
                 BBB 35/44, Inflammatory 22/33, Survival 41/68,
                 IonChannel 18/30.

PLOTS (../../figures_goplot/<Name>/<PNG|PDF>/)  — 31 figures, PNG + PDF
----------------------------------------------------------------------
  Consensus/
    GOBubble_consensus_top20percat  z-score vs -log10(padj); bubble area =
                                    gene count; faceted BP/MF/CC; top 20 per category
    GOCircle_consensus_top10        circular overview of the top 10 terms
                                    (outer ring = per-gene logFC scatter)
    GOHeat_consensus                genes x GO terms; GOHeat nlfc=1 (logFC) style
    GOChord_Consensus               gene<->GO-term chord
    GOChord_Consensus_KEGG          gene<->KEGG pathway chord
    GOChord_Consensus_Reactome      gene<->Reactome pathway chord
    GOChord_Consensus_AllEnrichments  top 2 terms from each of BP/MF/CC/KEGG/Reactome

  BBB/ Inflammatory/ Survival/ IonChannel/   (each)
    GOBubble_<name>_top5percat      top 5 enriched GO terms per category
    GOChord_<name>                  gene<->term chord (GO terms)
    GOHeat_<name>                   genes x terms; logFC (red/yellow/green)
    GOChord_<name>_KEGG             gene<->pathway chord using KEGG terms
    GOChord_<name>_Reactome         gene<->pathway chord using Reactome terms
    GOChord_<name>_AllEnrichments   all databases in one chord

  (The consensus chords cap genes to the top 5 per term by |logFC|. Note
   the chord legend is titled "GO Terms" by GOplot even for KEGG/Reactome.)

HOW TO REPRODUCE
----------------
  # 1. build the GOplot input tables (needs pandas + gseapy)
  python prepare_goplot_data.py
  # 2. render all plots (needs R + GOplot)
  Rscript goplot_analysis.R

ENVIRONMENT
-----------
  R 4.3.3 with GOplot 1.0.2 (+ ggplot2), installed into a user-space
  micromamba env — the same approach used for the GLP-1R project.
  PNG (200 dpi, white background so labels are legible) + vector PDF
  for every figure.

NOTES
-----
  * Gene symbols are upper-cased for GOplot's internal matching (so the
    plots show e.g. CAV1 rather than Cav1).
  * GOBubble, GOChord and GOCircle use the real GOplot functions. GOplot's
    own GOHeat drops its x-axis gene labels under modern ggplot2 (it calls
    scale_x_discrete on a numeric axis), so the GOHeat panels are
    reproduced with the same go_heat helper used in the GLP-1R project:
    tiles coloured by gene logFC (red high / yellow 0 / green low),
    unassigned tiles yellow, NO borders between squares; values clamped
    to +/- 2.
  * consensus_terms.csv is capped to the top 60 terms per GO category by
    adjusted p-value. GOplot's reduce_overlap() is O(n^2) and does not
    complete in reasonable time on the full 959-term consensus table; the
    plots only ever draw the top 20/10/7 terms, so the cap does not change
    any figure. The uncapped table is kept as consensus_terms_FULL.csv.
