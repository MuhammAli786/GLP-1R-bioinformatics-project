#!/usr/bin/env python3
"""
06_make_cnets.py
-------------------------------------------------------------
Concept-network (Cnet) plots in the cnet_style.py reference style
(with enlarged gene dots and thicker edges, per request).

Produces, for both thresholds (LFC0.2, LFC1):
  * Consensus    : all consensus genes vs all enrichment terms
  * BBB          : blood-brain-barrier / MMP gene+term set
  * Inflammatory : JAK-STAT3 inflammatory gene+term set
  * Survival     : PI3K/Akt pro-survival gene+term set

Saved as transparent PNG + PDF into:
  Plots/Cnet plots/<sub>/<LFCfolder>/<PDF|PNG>/
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
from cnet_gene_lists import (BBB_GENES, JAK_STAT3_GENES, PI3K_AKT_GENES,
                             ION_CHANNEL_BASE_GENES, KEYWORD_FILTERS)

BASE = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
PLOTS = os.path.join(BASE, "Plots", "Cnet plots")
THR_FOLDER = {"LFC02": "LFC0.2", "LFC05": "LFC0.5", "LFC1": "LFC1"}


def load_consensus(thr):
    p = os.path.join(DATA, f"consensus_{thr}.csv")
    return {r["gene_symbol"].upper(): r["gene_symbol"] for r in csv.DictReader(open(p))}

def load_enrichment(thr):
    p = os.path.join(DATA, f"enrichment_{thr}.csv")
    return [r for r in csv.DictReader(open(p)) if r.get("Adjusted P-value") and float(r["Adjusted P-value"]) < 0.05]

def load_lfc(thr):
    p = os.path.join(DATA, f"gene_lfc_comprehensive_{thr}.csv")
    return {r["gene_symbol"].upper(): float(r["mean_log2FC"]) for r in csv.DictReader(open(p))}

def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", t)
    t = re.sub(r"\s*R-HSA-\d+", "", t)
    return t[:52] + "..." if len(t) > 55 else t

def wrap_label(t, width=TERM_WRAP_WIDTH):
    return "\n".join(textwrap.wrap(t, width=width))


def build_cnet(ref_genes, term_filter, title, sub, thr, lfc_map, cons_map, enr_rows):
    input_genes = [cons_map[g.upper()] for g in ref_genes if g.upper() in cons_map] if ref_genes else list(cons_map.values())
    input_upper = {g.upper() for g in input_genes}
    if len(input_genes) < 3:
        print(f"  SKIP {sub}/{thr}: only {len(input_genes)} consensus genes"); return False

    if term_filter:
        pat = re.compile(term_filter)
        ft = [r for r in enr_rows if pat.search(r["Term"]) or
              len(set(r["Genes"].upper().split(";")) & input_upper) >= 2]
    else:
        ft = list(enr_rows)
    if not ft:
        print(f"  SKIP {sub}/{thr}: no matching terms"); return False
    ft.sort(key=lambda r: float(r["Adjusted P-value"]))

    covered, selected, remaining = set(), [], list(range(len(ft)))
    while len(selected) < MAX_TERMS and remaining:
        best = max(remaining, key=lambda i: (
            len(set(ft[i]["Genes"].upper().split(";")) & input_upper - covered),
            float(ft[i].get("Combined Score", 0))))
        tg = set(ft[best]["Genes"].upper().split(";")) & input_upper
        if not (tg - covered) and len(selected) >= 8:
            for idx in sorted(remaining, key=lambda i: float(ft[i].get("Combined Score", 0)), reverse=True):
                if len(selected) >= MAX_TERMS: break
                selected.append(idx); covered.update(set(ft[idx]["Genes"].upper().split(";")) & input_upper)
            break
        selected.append(best); covered.update(tg); remaining.remove(best)
    if not selected:
        print(f"  SKIP {sub}/{thr}: no terms after set-cover"); return False

    G = nx.Graph()
    for idx in selected:
        r = ft[idx]; t = clean_term(r["Term"])
        G.add_node(t, ntype="term", lib=r["Gene_set"])
        tg = set(r["Genes"].upper().split(";"))
        top = sorted(tg & input_upper, key=lambda g: -abs(lfc_map.get(g, 0)))
        for gu in top[:MAX_GENES_PER_TERM]:
            gn = cons_map.get(gu, gu)
            if gn not in G:
                G.add_node(gn, ntype="gene", lfc=lfc_map.get(gu, 0))
            G.add_edge(t, gn)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
    term_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "term"]
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "gene"]
    if len(gene_nodes) < 2:
        print(f"  SKIP {sub}/{thr}: <2 gene nodes"); return False
    print(f"  {sub}/{thr}: {len(term_nodes)} terms, {len(gene_nodes)} genes, {G.number_of_edges()} edges")

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
    ax.set_title(title, fontsize=24, fontweight="bold", pad=20); ax.axis("off")

    for fmt in ("PNG", "PDF"):
        d = os.path.join(PLOTS, sub, THR_FOLDER[thr], fmt); os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, f"Cnet_{sub}_{THR_FOLDER[thr]}.{fmt.lower()}"),
                    dpi=DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return True


def main():
    SPECS = [
        ("Consensus", None, None, "Consensus genes — GLP-1R agonist CNS meta-analysis"),
        ("BBB", BBB_GENES, KEYWORD_FILTERS["BBB"], "Blood-Brain Barrier / MMP Cnet"),
        ("Inflammatory", JAK_STAT3_GENES, KEYWORD_FILTERS["INFLAMMATORY"], "JAK-STAT3 Inflammatory Cnet"),
        ("Survival", PI3K_AKT_GENES, KEYWORD_FILTERS["SURVIVAL"], "PI3K/Akt Pro-Survival Cnet"),
        ("IonChannel", ION_CHANNEL_BASE_GENES, KEYWORD_FILTERS["ION_CHANNEL"], "Ion Channel Cnet"),
    ]
    thrs = [sys.argv[1]] if len(sys.argv) > 1 else ["LFC02", "LFC05", "LFC1"]
    for thr in thrs:
        cons_map = load_consensus(thr); enr = load_enrichment(thr); lfc = load_lfc(thr)
        print(f"=== {thr}: {len(cons_map)} consensus genes, {len(enr)} sig terms ===")
        for sub, genes, filt, title in SPECS:
            build_cnet(genes, filt, f"{title}  ({THR_FOLDER[thr]})", sub, thr, lfc, cons_map, enr)


if __name__ == "__main__":
    main()
