#!/usr/bin/env python3
"""
06c_bbb_freq1_cnet.py
-------------------------------------------------------------
A SEPARATE Blood-Brain-Barrier consensus concept-network that relaxes
the consensus gene-frequency requirement from >= 2 groups to >= 1
group (i.e. show every BBB gene PRESENT/significant in any group, not
only recurrent ones). Same cnet_style.py styling as the main cnets.

Generated for each threshold (LFC0.2, LFC0.5, LFC1). Saved alongside
the standard BBB cnet with a _Freq1 suffix:
  Plots/Cnet plots/BBB/<LFC>/<PDF|PNG>/Cnet_BBB_Freq1_<LFC>.*
"""
import os, re, csv, textwrap, sys, time
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnet_style import *
from cnet_gene_lists import BBB_GENES, KEYWORD_FILTERS

DATABASES = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']

BASE = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
PLOTS = os.path.join(BASE, "Plots", "Cnet plots", "BBB")
THR_FOLDER = {"LFC02": "LFC0.2", "LFC05": "LFC0.5", "LFC1": "LFC1"}


def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", t)
    t = re.sub(r"\s*R-HSA-\d+", "", t)
    return t[:52] + "..." if len(t) > 55 else t

def wrap_label(t, width=TERM_WRAP_WIDTH):
    return "\n".join(textwrap.wrap(t, width=width))


def presence_maps(thr):
    """genes present (significant) in >= 1 group, with mean log2FC."""
    m = pd.read_csv(os.path.join(DATA, f"master_deg_{thr}.csv"))
    g = m.groupby("symbol").agg(mean_lfc=("log2FC", "mean")).reset_index()
    cons_map = {s.upper(): s for s in g["symbol"]}
    lfc_map = {r.symbol.upper(): float(r.mean_lfc) for r in g.itertuples()}
    return cons_map, lfc_map


def enrich_present(genes):
    """Run Enrichr on the present BBB gene list so term gene-lists include
    the present-only genes (e.g. Mmp2, Mmp14)."""
    rows = []
    for db in DATABASES:
        for attempt in range(4):
            try:
                res = gp.enrichr(gene_list=genes, gene_sets=db, organism="mouse", no_plot=True).results
                res["Gene_set"] = db
                rows.append(res); break
            except Exception:
                time.sleep(3 * (attempt + 1))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def build(thr):
    cons_map, lfc_map = presence_maps(thr)
    input_genes = [cons_map[g.upper()] for g in BBB_GENES if g.upper() in cons_map]
    input_upper = {g.upper() for g in input_genes}
    if len(input_genes) < 3:
        print(f"  {thr}: only {len(input_genes)} BBB genes present (>=1); skip"); return

    edf = enrich_present(input_genes)
    if edf.empty:
        print(f"  {thr}: enrichment failed"); return
    enr = edf.to_dict("records")
    pat = re.compile(KEYWORD_FILTERS["BBB"])
    # keep BBB-keyword terms OR significant terms OR terms overlapping >=2 input genes
    ft = [r for r in enr if pat.search(str(r["Term"])) or float(r["Adjusted P-value"]) < 0.05
          or len(set(str(r["Genes"]).upper().split(";")) & input_upper) >= 2]
    if not ft:
        print(f"  {thr}: no matching BBB terms"); return
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

    G = nx.Graph()
    for idx in selected:
        r = ft[idx]; t = clean_term(r["Term"])
        G.add_node(t, ntype="term", lib=r["Gene_set"])
        tg = set(r["Genes"].upper().split(";")) & input_upper
        for gu in sorted(tg, key=lambda g: -abs(lfc_map.get(g, 0)))[:MAX_GENES_PER_TERM]:
            gn = cons_map.get(gu, gu)
            if gn not in G:
                G.add_node(gn, ntype="gene", lfc=lfc_map.get(gu, 0))
            G.add_edge(t, gn)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
    term_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "term"]
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "gene"]
    if len(gene_nodes) < 2:
        print(f"  {thr}: <2 BBB gene nodes; skip"); return
    print(f"  {thr}: {len(term_nodes)} terms, {len(gene_nodes)} BBB genes (>=1 group), {G.number_of_edges()} edges")

    pos = nx.spring_layout(G, k=K_SPRING, iterations=ITERATIONS, seed=SEED)
    fig, ax = plt.subplots(figsize=FIGSIZE); fig.patch.set_alpha(0); ax.set_facecolor("none")
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=EDGE_ALPHA, width=EDGE_WIDTH, edge_color=EDGE_COLOR)
    gene_sizes = [max(GENE_MIN_SIZE, G.degree(g) * GENE_DEGREE_SCALE) for g in gene_nodes]
    gene_xy = np.array([pos[g] for g in gene_nodes])
    sc = ax.scatter(gene_xy[:, 0], gene_xy[:, 1],
                    c=[np.clip(G.nodes[g].get("lfc", 0), -VMAX, VMAX) for g in gene_nodes],
                    cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, s=gene_sizes, zorder=3,
                    edgecolors=GENE_EDGE_COLOR, linewidths=GENE_EDGE_WIDTH)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02)
    cbar.set_label("mean log2FC", fontsize=20); cbar.ax.tick_params(labelsize=14)
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get("lib", ""), "#888888"); xy = pos[t]
        ax.scatter(xy[0], xy[1], s=TERM_SQUARE_SIZE, c=color, zorder=4,
                   edgecolors=TERM_EDGE_COLOR, linewidths=TERM_EDGE_WIDTH, marker=TERM_MARKER)
        ax.text(xy[0], xy[1] + TERM_LABEL_OFFSET, wrap_label(t), fontsize=TERM_LABEL_FONTSIZE,
                fontweight="bold", color="white", ha="center", va="bottom", zorder=5,
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
    ax.set_title(f"Blood-Brain Barrier / MMP Cnet — gene presence >=1 group  ({THR_FOLDER[thr]})",
                 fontsize=24, fontweight="bold", pad=20)
    ax.axis("off")
    for fmt in ("PNG", "PDF"):
        d = os.path.join(PLOTS, THR_FOLDER[thr], fmt); os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, f"Cnet_BBB_Freq1_{THR_FOLDER[thr]}.{fmt.lower()}"),
                    dpi=DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)


THRS = [sys.argv[1]] if len(sys.argv) > 1 else ["LFC02", "LFC05", "LFC1"]
for thr in THRS:
    build(thr)
print("done")
