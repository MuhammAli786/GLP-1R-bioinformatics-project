#!/usr/bin/env python3
"""Pathway-restricted HIBI Cnets, a port of the GLP-1R 06d_focused_pathway_cnet.py.

Terms come from enriching only that pathway's HIBI consensus genes, as the
GOplot chord does, rather than keyword-filtering a whole-consensus enrichment.
Uses cnet_style.py styling and the four curated gene sets in cnet_gene_lists.py
(BBB/MMP, JAK-STAT3 inflammatory, PI3K-Akt pro-survival, ion channel); these are
the HIBI equivalents of GLP-1R Fig08A-D.
Outputs: figures/Cnet_<Pathway>_FocusedTerms_LFC0.2.{png,pdf}
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

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
from cnet_gene_lists import BBB_GENES, JAK_STAT3_GENES, PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES

DATA = BASE + "/mnt/outputs/hibi_data/meta_analysis"
OUTDIR = BASE + "/mnt/outputs/hibi_figures"
THR_FOLDER = {"LFC02": "LFC0.2", "LFC05": "LFC0.5", "LFC1": "LFC1"}
DATABASES = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']

PATHWAYS = {"BBB": BBB_GENES, "Inflammatory": JAK_STAT3_GENES,
            "Survival": PI3K_AKT_GENES, "IonChannel": ION_CHANNEL_BASE_GENES}
TITLES = {"BBB": "Blood-Brain Barrier / MMP", "Inflammatory": "JAK-STAT3 Inflammatory",
          "Survival": "PI3K/Akt Pro-Survival", "IonChannel": "Ion Channel"}


def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", t)
    t = re.sub(r"\s*R-HSA-\d+", "", t)
    return t[:52] + "..." if len(t) > 55 else t

def wrap_label(t, width=TERM_WRAP_WIDTH):
    return "\n".join(textwrap.wrap(t, width=width))


def build(pathway, thr="LFC02"):
    cons = pd.read_csv(os.path.join(DATA, f"consensus_{thr}.csv"))
    cons_map = {r.gene_symbol.upper(): r.gene_symbol for r in cons.itertuples()}
    lfc_map = {r.gene_symbol.upper(): float(r.mean_log2FC) for r in cons.itertuples()}
    genes = [cons_map[g.upper()] for g in PATHWAYS[pathway] if g.upper() in cons_map]
    input_upper = {g.upper() for g in genes}
    if len(genes) < 3:
        print(f"  {pathway}/{thr}: only {len(genes)} consensus genes; skip"); return
    print(f"  {pathway}/{thr}: {len(genes)} pathway genes in HIBI consensus")

    # Enrich only this pathway's consensus genes
    rows = []
    for db in DATABASES:
        for attempt in range(4):
            try:
                res = gp.enrichr(gene_list=genes, gene_sets=db, organism="mouse", no_plot=True).results
                res["Gene_set"] = db
                rows.append(res); break
            except Exception:
                time.sleep(3 * (attempt + 1))
    edf = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    enr = [r for r in edf.to_dict("records") if float(r["Adjusted P-value"]) < 0.05]
    if not enr:
        print(f"  {pathway}/{thr}: no significant focused terms"); return
    enr.sort(key=lambda r: float(r["Adjusted P-value"]))

    # Select top terms by adjusted p-value, as the GOplot chord does, dropping
    # near-duplicate terms (pathway-gene Jaccard >= 0.8 with one already chosen)
    # so the network stays readable.
    MAXT = 12
    selected, sel_sets = [], []
    for i in range(len(enr)):
        tg = set(enr[i]["Genes"].upper().split(";")) & input_upper
        if len(tg) < 2:
            continue
        if any(len(tg & s) / max(1, len(tg | s)) >= 0.8 for s in sel_sets):
            continue
        selected.append(i); sel_sets.append(tg)
        if len(selected) >= MAXT:
            break

    G = nx.Graph()
    for idx in selected:
        r = enr[idx]; t = clean_term(r["Term"])
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
        print(f"  {pathway}/{thr}: <2 gene nodes; skip"); return
    print(f"  {pathway}/{thr}: {len(term_nodes)} focused terms, {len(gene_nodes)} genes, {G.number_of_edges()} edges")

    pos = nx.spring_layout(G, k=K_SPRING, iterations=ITERATIONS, seed=SEED)
    fig, ax = plt.subplots(figsize=FIGSIZE); fig.patch.set_alpha(0); ax.set_facecolor("none")
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=EDGE_ALPHA, width=EDGE_WIDTH, edge_color=EDGE_COLOR)
    gene_sizes = [max(GENE_MIN_SIZE, G.degree(g) * GENE_DEGREE_SCALE) for g in gene_nodes]
    gxy = np.array([pos[g] for g in gene_nodes])
    sc = ax.scatter(gxy[:, 0], gxy[:, 1], c=[np.clip(G.nodes[g].get("lfc", 0), -VMAX, VMAX) for g in gene_nodes],
                    cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, s=gene_sizes, zorder=3,
                    edgecolors=GENE_EDGE_COLOR, linewidths=GENE_EDGE_WIDTH)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02); cbar.set_label("mean log2FC", fontsize=20)
    cbar.ax.tick_params(labelsize=14)
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get("lib", ""), "#888888"); xy = pos[t]
        ax.scatter(xy[0], xy[1], s=TERM_SQUARE_SIZE, c=color, zorder=4,
                   edgecolors=TERM_EDGE_COLOR, linewidths=TERM_EDGE_WIDTH, marker=TERM_MARKER)
        ax.text(xy[0], xy[1] + TERM_LABEL_OFFSET, wrap_label(t), fontsize=TERM_LABEL_FONTSIZE,
                fontweight="bold", color="white", ha="center", va="bottom", zorder=5,
                bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=TERM_BOX_ALPHA, edgecolor="none"))
    for g in gene_nodes:
        xy = pos[g]
        ax.text(xy[0], xy[1] - GENE_LABEL_OFFSET, g, fontsize=GENE_LABEL_FONTSIZE, fontweight="bold",
                color="white", ha="center", va="top", zorder=5,
                bbox=dict(boxstyle="round,pad=0.2", facecolor=GENE_BOX_COLOR, alpha=GENE_BOX_ALPHA, edgecolor="none"))
    libs_used = {G.nodes[t].get("lib") for t in term_nodes}
    leg = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=12, label=LIB_LABELS.get(lib, lib))
           for lib, c in LIB_COLORS.items() if lib in libs_used]
    if leg:
        ax.legend(handles=leg, loc="lower left", fontsize=20, framealpha=0.8, facecolor="white")
    ax.set_title(f"HIBI - {TITLES[pathway]} Cnet - focused enrichment  ({THR_FOLDER[thr]})",
                 fontsize=24, fontweight="bold", pad=20); ax.axis("off")
    os.makedirs(OUTDIR, exist_ok=True)
    for fmt in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"Cnet_{pathway}_FocusedTerms_{THR_FOLDER[thr]}.{fmt}"),
                    dpi=DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  saved Cnet_{pathway}_FocusedTerms_{THR_FOLDER[thr]} (png + pdf)")


if __name__ == "__main__":
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(PATHWAYS)
    for p in targets:
        build(p, "LFC02")
        time.sleep(2)
    print("done")
