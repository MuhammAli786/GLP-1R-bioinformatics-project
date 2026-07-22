Cnet Plot Generation Scripts
============================

These scripts generate all concept network (Cnet) plots for the
GLP-1R agonist CNS meta-analysis.

Files
-----
  cnet_style.py          Style configuration (colours, sizes, fonts, layout)
  cnet_gene_lists.py     Curated pathway gene lists and keyword filters
  cnet_engine.py         Core plot builder (data loading, graph, rendering)
  generate_all_cnets.py  Master script — generates every Cnet plot

Quick Start
-----------
  cd "Final analysis/Scripts"
  python generate_all_cnets.py \
      --data-dir /path/to/csv/data \
      --out-dir  /path/to/Final_analysis \
      --workspace /path/to/Bulk_RNA_sequencing

Dependencies
------------
  Python 3.8+
  numpy, matplotlib, networkx

  Install: pip install numpy matplotlib networkx

Style Reference
---------------
Plot style is based on 02_build_cnet_plots.py (v2 style):
  - Gene nodes: small scatter dots, coloured by mean log2FC (RdBu_r)
  - Gene labels: white text on dark box (#333333), positioned below nodes
  - Term nodes: coloured squares (by enrichment database)
  - Term labels: white text on coloured box, positioned above nodes,
    wrapped at 24 characters, GO/Reactome IDs stripped
  - Edges: thin (0.7 width), subtle (0.2 alpha)
  - Colorbar: symmetric, capped at +/- 2.0
  - Layout: spring layout, k=3.5, 150 iterations, seed=42
  - Output: 200 dpi PNG + PDF, transparent backgrounds

Thresholds
----------
  LFC0   padj < 0.05 only
  LFC02  |log2FC| >= 0.2 + padj < 0.05
  LFC05  |log2FC| >= 0.5 + padj < 0.05
  LFC1   |log2FC| >= 1.0 + padj < 0.05

Plots Generated
---------------
  Consensus-Based (4 plots):
    Cnet_Consensus_{LFC0, LFC02, LFC05, LFC1}

  Pathway-Based (12 plots):
    Cnet_BBB_MMP_{LFC0, LFC02, LFC05, LFC1}
    Cnet_JAKSTAT3_Inflammatory_{LFC0, LFC02, LFC05, LFC1}
    Cnet_AktPI3K_ProSurvival_{LFC0, LFC02, LFC05, LFC1}

  Specific Pathway (up to 16 plots):
    Cnet_BBB_{LFC0, LFC02, LFC05, LFC1}
    Cnet_Inflammatory_{LFC0, LFC02, LFC05, LFC1}
    Cnet_Survival_{LFC0, LFC02, LFC05, LFC1}
    Cnet_IonChannel_{LFC0, LFC02}  (LFC05/LFC1 skipped: too few genes)
