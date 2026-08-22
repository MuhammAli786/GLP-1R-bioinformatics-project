#!/usr/bin/env python3
"""Consensus-signature Cnet built from the genes recurring in the most groups (ranked by n_groups, not by log2FC or enrichment strength), over the GLP-1R-only groups at threshold LFC0.2.

Data/master_deg_LFC02.csv plus a fresh Enrichr run -> Plots/Cnet plots/Consensus/LFC0.2/<PDF|PNG>/Cnet_Consensus_TopFrequentGenes_LFC0.2.* and GOplot input tables GOPLOT analyis/data/freqsig_{terms,kegg_terms,reactome_terms,genes}.csv.
Gene nodes are sized by recurrence and coloured by mean log2FC; terms come from enriching only these frequent genes.
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, re, textwrap, sys, time
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

BASE = BASE + "/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
OUTDIR = os.path.join(BASE, "Plots", "Cnet plots", "Consensus", "LFC0.2")
GO_DATA = BASE + "/mnt/Bulk RNA sequencing/GOPLOT analyis/data"
GROUPS = BASE + "/work/gsea_groups.csv"   # 31 GLP-1R-only groups
DBS = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
       'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']
CATMAP = {'GO_Biological_Process_2023': 'BP', 'GO_Molecular_Function_2023': 'MF',
          'GO_Cellular_Component_2023': 'CC', 'KEGG_2021_Human': 'KEGG', 'Reactome_2022': 'Reactome'}
TOP_N = 30
MAXT = 12


def clean(t):
    """Strip GO and Reactome accessions from a term."""
    t = re.sub(r"\s*\(GO:\d+\)", "", str(t)); return re.sub(r"\s*R-HSA-\d+", "", t)
def wrap(t):
    """Clean, truncate and hard-wrap a term label for use as node text."""
    return "\n".join(textwrap.wrap(clean(t)[:52], TERM_WRAP_WIDTH))


def freq_genes():
    """Return the TOP_N genes recurring in the most GLP-1R-only groups, with mean log2FC."""
    keep = set(pd.read_csv(GROUPS)["group"])
    m = pd.read_csv(os.path.join(DATA, "master_deg_LFC02.csv"))
    m = m[m["group"].isin(keep)]
    agg = (m.groupby("symbol").agg(n_groups=("group", "nunique"), mean_lfc=("log2FC", "mean"))
           .reset_index().sort_values("n_groups", ascending=False))
    return agg.head(TOP_N)


def parse_id(term):
    """Split an Enrichr term into its (GO or Reactome accession, bare name)."""
    mm = re.search(r'\(GO:(\d+)\)', term)
    if mm: return f"GO:{mm.group(1)}", re.sub(r'\s*\(GO:\d+\)', '', term).strip()
    mm = re.search(r'\bR-HSA-(\d+)\b', term)
    if mm: return f"R-HSA-{mm.group(1)}", re.sub(r'\s*R-HSA-\d+', '', term).strip()
    return "", term.strip()


def main():
    """Write the GOplot tables and draw the frequency-ranked consensus Cnet."""
    fg = freq_genes()
    genes = fg["symbol"].tolist()
    lfc = {g.upper(): l for g, l in zip(fg["symbol"], fg["mean_lfc"])}
    nfreq = {g.upper(): int(n) for g, n in zip(fg["symbol"], fg["n_groups"])}
    disp = {g.upper(): g for g in genes}
    input_upper = set(disp)
    print(f"Top {TOP_N} most frequent genes (n_groups {fg['n_groups'].min()}-{fg['n_groups'].max()})")

    rows = []
    for db in DBS:
        for _ in range(4):
            try:
                r = gp.enrichr(gene_list=genes, gene_sets=db, organism="mouse", no_plot=True).results
                r["Gene_set"] = db; rows.append(r); break
            except Exception:
                time.sleep(3)
    enr = pd.concat(rows, ignore_index=True)

    # GOplot input tables for the matching GOChord. Per category, keep the terms with
    # padj < 0.05, or the 6 most significant if fewer than 2 reach significance.
    os.makedirs(GO_DATA, exist_ok=True)
    trows = []
    for _, r in enr.iterrows():
        tid, name = parse_id(str(r["Term"]))
        trows.append({"Category": CATMAP[r["Gene_set"]], "ID": tid, "Term": name,
                      "Genes": str(r["Genes"]).replace(";", ", "), "adj_pval": float(r["Adjusted P-value"])})
    tdf = pd.DataFrame(trows)
    keep = []
    for cat in ["BP", "MF", "CC", "KEGG", "Reactome"]:
        c = tdf[tdf["Category"] == cat]; s = c[c["adj_pval"] < 0.05]
        keep.append(s if len(s) >= 2 else c.nsmallest(6, "adj_pval"))
    tdf = pd.concat(keep, ignore_index=True)
    tdf[tdf["Category"].isin(["BP", "MF", "CC"])].to_csv(f"{GO_DATA}/freqsig_terms.csv", index=False)
    tdf[tdf["Category"] == "KEGG"].to_csv(f"{GO_DATA}/freqsig_kegg_terms.csv", index=False)
    tdf[tdf["Category"] == "Reactome"].to_csv(f"{GO_DATA}/freqsig_reactome_terms.csv", index=False)
    pd.DataFrame({"ID": genes, "logFC": [lfc[g.upper()] for g in genes]}).to_csv(f"{GO_DATA}/freqsig_genes.csv", index=False)
    print("wrote GOplot freqsig_* tables")

    # Select terms by greedy set cover over the frequent genes, breaking ties on adjusted p-value.
    er = [r for r in enr.to_dict("records") if float(r["Adjusted P-value"]) < 0.05]
    if len(er) < 3:
        er = sorted(enr.to_dict("records"), key=lambda r: float(r["Adjusted P-value"]))[:MAXT]
    er.sort(key=lambda r: float(r["Adjusted P-value"]))
    covered, selected, remaining = set(), [], list(range(len(er)))
    while len(selected) < MAXT and remaining:
        best = max(remaining, key=lambda i: (len(set(er[i]["Genes"].upper().split(";")) & input_upper - covered),
                                             -float(er[i]["Adjusted P-value"])))
        tg = set(er[best]["Genes"].upper().split(";")) & input_upper
        if not (tg - covered) and len(selected) >= 6: break
        selected.append(best); covered.update(tg); remaining.remove(best)

    G = nx.Graph()
    for idx in selected:
        r = er[idx]; t = clean(r["Term"])[:52]
        G.add_node(("T", t), ntype="term", lib=r["Gene_set"])
        for gu in (set(r["Genes"].upper().split(";")) & input_upper):
            if ("G", gu) not in G: G.add_node(("G", gu), ntype="gene")
            G.add_edge(("T", t), ("G", gu))
    # Every frequent gene gets a node, even one that appears in no selected term.
    for gu in input_upper:
        if ("G", gu) not in G:
            G.add_node(("G", gu), ntype="gene")
    # Drop only term nodes with no genes; keep every frequent gene node.
    G.remove_nodes_from([n for n, d in G.nodes(data=True) if d["ntype"] == "term" and G.degree(n) == 0])
    term_nodes = [n for n, d in G.nodes(data=True) if d["ntype"] == "term"]
    gene_nodes = [n for n, d in G.nodes(data=True) if d["ntype"] == "gene"]
    print(f"  cnet: {len(term_nodes)} terms, {len(gene_nodes)} frequent genes, {G.number_of_edges()} edges")

    pos = nx.spring_layout(G, k=K_SPRING, iterations=ITERATIONS, seed=SEED)
    fig, ax = plt.subplots(figsize=FIGSIZE); fig.patch.set_alpha(0); ax.set_facecolor("none")
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=EDGE_ALPHA, width=EDGE_WIDTH, edge_color=EDGE_COLOR)
    # Gene node size encodes recurrence (n_groups).
    gsz = [GENE_MIN_SIZE + nfreq.get(g[1], 1) * 90 for g in gene_nodes]
    gxy = np.array([pos[g] for g in gene_nodes])
    sc = ax.scatter(gxy[:, 0], gxy[:, 1], c=[np.clip(lfc.get(g[1], 0), -VMAX, VMAX) for g in gene_nodes],
                    cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, s=gsz, zorder=3,
                    edgecolors=GENE_EDGE_COLOR, linewidths=GENE_EDGE_WIDTH)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02); cbar.set_label("mean log2FC", fontsize=20)
    cbar.ax.tick_params(labelsize=14)
    for t in term_nodes:
        col = LIB_COLORS.get(G.nodes[t]["lib"], "#888"); xy = pos[t]
        ax.scatter(xy[0], xy[1], s=TERM_SQUARE_SIZE, c=col, zorder=4, edgecolors=TERM_EDGE_COLOR,
                   linewidths=TERM_EDGE_WIDTH, marker="s")
        ax.text(xy[0], xy[1] + TERM_LABEL_OFFSET, wrap(t[1]), fontsize=TERM_LABEL_FONTSIZE, fontweight="bold",
                color="white", ha="center", va="bottom", zorder=5,
                bbox=dict(boxstyle="round,pad=0.3", facecolor=col, alpha=TERM_BOX_ALPHA, edgecolor="none"))
    for g in gene_nodes:
        xy = pos[g]
        ax.text(xy[0], xy[1] - GENE_LABEL_OFFSET, f"{disp.get(g[1],g[1])}", fontsize=GENE_LABEL_FONTSIZE,
                fontweight="bold", color="white", ha="center", va="top", zorder=5,
                bbox=dict(boxstyle="round,pad=0.2", facecolor=GENE_BOX_COLOR, alpha=GENE_BOX_ALPHA, edgecolor="none"))
    libs = {G.nodes[t]["lib"] for t in term_nodes}
    leg = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=12, label=LIB_LABELS.get(l, l))
           for l, c in LIB_COLORS.items() if l in libs]
    if leg: ax.legend(handles=leg, loc="lower left", fontsize=20, framealpha=0.8, facecolor="white")
    ax.set_title("Consensus signature Cnet - most frequently recurring genes  (LFC 0.2)\n"
                 "(node size = number of groups the gene recurs in)", fontsize=22, fontweight="bold", pad=20)
    ax.axis("off")
    for fmt in ("PNG", "PDF"):
        d = os.path.join(OUTDIR, fmt); os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, f"Cnet_Consensus_TopFrequentGenes_LFC0.2.{fmt.lower()}"),
                    dpi=DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("saved Cnet_Consensus_TopFrequentGenes_LFC0.2")


if __name__ == "__main__":
    main()
