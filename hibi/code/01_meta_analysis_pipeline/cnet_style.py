"""
cnet_style.py — Shared style configuration for all Cnet plots
=============================================================
Centralises colours, sizes, fonts, layout parameters, and library
metadata so every plot in the pipeline looks identical.

Import with:
    from cnet_style import *
"""
# Colorbar / LFC 
VMAX = 2.0                          

# Figure 
FIGSIZE = (30, 26)                   # Inches
DPI = 200                            # PNG resolution

# Edges 
EDGE_WIDTH = 2.2
EDGE_ALPHA = 0.45
EDGE_COLOR = '#555555'

# Gene nodes 
GENE_MIN_SIZE = 320                   # Scatter size for degree-1 genes
GENE_DEGREE_SCALE = 110               # Added per degree
GENE_EDGE_COLOR = 'white'
GENE_EDGE_WIDTH = 0.5
GENE_LABEL_FONTSIZE = 20
GENE_LABEL_OFFSET = 0.03            # Below node (subtracted from y)
GENE_BOX_COLOR = '#333333'
GENE_BOX_ALPHA = 0.85

# Term (pathway) nodes 
TERM_SQUARE_SIZE = 350               # Marker size for square node
TERM_MARKER = 's'
TERM_EDGE_COLOR = 'white'
TERM_EDGE_WIDTH = 1.5
TERM_LABEL_FONTSIZE = 20
TERM_LABEL_OFFSET = 0.035           # Above node (added to y)
TERM_BOX_ALPHA = 0.9
TERM_WRAP_WIDTH = 24                 # Characters before wrapping

# Layout 
K_SPRING = 3.5                       # Spring constant (higher = more spacing)
ITERATIONS = 150                     # Spring layout iterations
SEED = 42                            # Reproducible layout

# Term selection 
MAX_TERMS = 15                       # Max enrichment terms per plot
MAX_GENES_PER_TERM = 20              # Cap genes per term edge list

# Library colours (enrichment databases)
LIB_COLORS = {
    'GO_Biological_Process_2023': '#2196F3',   # Blue
    'GO_Molecular_Function_2023': '#9C27B0',   # Purple
    'GO_Cellular_Component_2023': '#FF9800',   # Orange
    'KEGG_2021_Human':            '#4CAF50',   # Green
    'Reactome_2022':              '#F44336',   # Red
}

LIB_LABELS = {
    'GO_Biological_Process_2023': 'GO Biological Process',
    'GO_Molecular_Function_2023': 'GO Molecular Function',
    'GO_Cellular_Component_2023': 'GO Cellular Component',
    'KEGG_2021_Human':            'KEGG',
    'Reactome_2022':              'Reactome',
}

# Threshold display names 
THRESHOLD_LABELS = {
    'LFC0':  'LFC 0',
    'LFC02': 'LFC 0.2',
    'LFC05': 'LFC 0.5',
    'LFC1':  'LFC 1.0',
}
