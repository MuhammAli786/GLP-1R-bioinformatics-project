#!/usr/bin/env python3
"""
06b_consensus_top5_cnet.py
-------------------------------------------------------------
A SEPARATE consensus concept-network for the LFC0.2 threshold that
includes ONLY the top 5 most-enriched terms (smallest adjusted
p-value across all databases). Same cnet_style.py styling as the
main cnets (enlarged gene dots, thick edges).

Output:
  Plots/Cnet plots/Consensus/LFC0.2/<PDF|PNG>/Cnet_Consensus_Top5Terms_LFC0.2.*
"""
import os, re, csv, textwrap, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnet_style import *

BASE = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
OUTDIR = os.path.join(BASE, "Plots", "Cnet plots", "Consensus", "LFC0.2")
THR = "LFC02"
N_TERMS = 5


def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", t)
    t = re.sub(r"\s*R-HSA-\d+", "", t)
    return t[:52] + "..." if len(t) > 55 else t

def wrap_label(t, width=TERM_WRAP_WIDTH):
    return "\n".join(textwrap.wrap(t, width=width))


# Load data
cons_map = {r["gene_symbol"].upper(): r["gene_symbol"]
            for r in csv.DictReader(open(os.path.join(DATA, f"consensus_{THR}.csv")))}
lfc_map = {r["gene_symbol"].upper(): float(r["mean_log2FC"])
           for r in csv.DictReader(open(os.path.join(DATA, f"gene_lfc_comprehensive_{THR}.csv")))}
enr = [r for r in csv.DictReader(open(os.path.join(DATA, f"enrichment_{THR}.csv")))
       if r.get("Adjusted P-value") and float(r["Adjusted P-value"]) < 0.05]

# Top 5 enriched terms by adjusted p-value (most significant)
enr.sort(key=lambda r: float(r["Adjusted P-value"]))
selected = enr[:N_TERMS]
input_upper = set(cons_map.keys())

# Build graph
G = nx.Graph()
for r in selected:
    t = clean_term(r["Term"])
    G.add_node(t, ntype="term", lib=r["Gene_set"])
    tg = set(r["Genes"].upper().split(";")) & input_upper
    top_genes = sorted(tg, key=lambda g: -abs(lfc_map.get(g, 0)))
    for gu in top_genes[:MAX_GENES_PER_TERM]:
        gn = cons_map.get(gu, gu)
        if gn not in G:
            G.add_node(gn, ntype="gene", lfc=lfc_map.get(gu, 0))
        G.add_edge(t, gn)

G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
term_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "term"]
gene_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "gene"]
print(f"Top-{N_TERMS} terms: {len(term_nodes)} terms, {len(gene_nodes)} genes, {G.number_of_edges()} edges")
for r in selected:
    print(f"  - {clean_term(r['Term'])}  (padj={float(r['Adjusted P-value']):.1e}, {r['Gene_set']})")

pos = nx.spring_layout(G, k=K_SPRING, iterations=ITERATIONS, seed=SEED)
fig, ax = plt.subplots(figsize=FIGSIZE)
fig.patch.set_alpha(0); ax.set_facecolor("none")
nx.draw_networkx_edges(G, pos, ax=ax, alpha=EDGE_ALPHA, width=EDGE_WIDTH, edge_color=EDGE_COLOR)
gene_lfcs = [G.nodes[g].get("lfc", 0) for g in gene_nodes]
gene_sizes = [max(GENE_MIN_SIZE, G.degree(g) * GENE_DEGREE_SCALE) for g in gene_nodes]
gene_xy = np.array([pos[g] for g in gene_nodes])
sc = ax.scatter(gene_xy[:, 0], gene_xy[:, 1],
                c=[np.clip(l, -VMAX, VMAX) for l in gene_lfcs],
                cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, s=gene_sizes, zorder=3,
                edgecolors=GENE_EDGE_COLOR, linewidths=GENE_EDGE_WIDTH)
cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02)
cbar.set_label("mean log2FC", fontsize=20); cbar.ax.tick_params(labelsize=14)
for t in term_nodes:
    color = LIB_COLORS.get(G.nodes[t].get("lib", ""), "#888888"); xy = pos[t]
    ax.scatter(xy[0], xy[1], s=TERM_SQUARE_SIZE, c=color, zorder=4,
               edgecolors=TERM_EDGE_COLOR, linewidths=TERM_EDGE_WIDTH, marker=TERM_MARKER)
    ax.text(xy[0], xy[1] + TERM_LABEL_OFFSET, wrap_label(t),
            fontsize=TERM_LABEL_FONTSIZE, fontweight="bold", color="white",
            ha="center", va="bottom", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=TERM_BOX_ALPHA, edgecolor="none"))
for g in gene_nodes:
    xy = pos[g]
    ax.text(xy[0], xy[1] - GENE_LABEL_OFFSET, g, fontsize=GENE_LABEL_FONTSIZE,
            fontweight="bold", color="white", ha="center", va="top", zorder=5,
            bbox=dict(boxstyle="round,pad=0.2", facecolor=GENE_BOX_COLOR, alpha=GENE_BOX_ALPHA, edgecolor="none"))
libs_used = {G.nodes[t].get("lib") for t in term_nodes}
legend_els = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=12,
                     label=LIB_LABELS.get(lib, lib)) for lib, c in LIB_COLORS.items() if lib in libs_used]
if legend_els:
    ax.legend(handles=legend_els, loc="lower left", fontsize=20, framealpha=0.8, facecolor="white")
ax.set_title("Consensus Cnet — Top 5 Enriched Terms  (LFC 0.2)", fontsize=24, fontweight="bold", pad=20)
ax.axis("off")

for fmt in ("PNG", "PDF"):
    d = os.path.join(OUTDIR, fmt); os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, f"Cnet_Consensus_Top5Terms_LFC0.2.{fmt.lower()}"),
                dpi=DPI, bbox_inches="tight", transparent=True)
plt.close(fig)
print("saved Cnet_Consensus_Top5Terms_LFC0.2 (png + pdf)")
