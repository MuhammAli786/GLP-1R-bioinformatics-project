"""
cnet_engine.py — Core Cnet plot builder
=======================================
Provides:
  - Data loading (consensus genes, enrichment, LFC map)
  - Term cleaning (strip GO/Reactome IDs)
  - Greedy set-cover term selection
  - NetworkX graph construction
  - Matplotlib figure rendering in the reference style
"""
import csv, re, os, textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import networkx as nx

from cnet_style import *


def load_consensus(data_dir, threshold):
    """Return {GENE_UPPER: original_case_name} for a threshold."""
    path = os.path.join(data_dir, f'final_consensus_{threshold}.csv')
    return {r['gene_symbol'].upper(): r['gene_symbol']
            for r in csv.DictReader(open(path))}

def load_enrichment(data_dir, threshold):
    """Return list of significant (padj < 0.05) enrichment rows."""
    path = os.path.join(data_dir, f'final_enrichment_{threshold}.csv')
    return [r for r in csv.DictReader(open(path))
            if float(r['Adjusted P-value']) < 0.05]

def load_lfc(data_dir):
    """Return {GENE_UPPER: mean_log2FC}."""
    path = os.path.join(data_dir, 'final_gene_lfc_COMPREHENSIVE.csv')
    return {r['gene_symbol'].upper(): float(r['mean_log2FC'])
            for r in csv.DictReader(open(path))}

# Helpers
def clean_term(t):
    """Remove GO / Reactome IDs and truncate long names."""
    t = re.sub(r'\s*\(GO:\d+\)', '', t)
    t = re.sub(r'\s*R-HSA-\d+', '', t)
    return t[:52] + '...' if len(t) > 55 else t

def wrap_label(t, width=TERM_WRAP_WIDTH):
    return '\n'.join(textwrap.wrap(t, width=width))
# Core builder

def build_cnet(ref_genes, term_filter, title, fname, out_dir,
               lfc_map, cons_map, enr_rows):

    input_genes = [cons_map[g.upper()] for g in ref_genes
                   if g.upper() in cons_map]
    input_upper = {g.upper() for g in input_genes}

    if len(input_genes) < 3:
        print(f"  SKIP {fname}: only {len(input_genes)} consensus genes")
        return False

    #  Filter enrichment terms 
    if term_filter:
        pat = re.compile(term_filter)
        ft = [r for r in enr_rows
              if pat.search(r['Term'])
              or len(set(r['Genes'].upper().split(';')) & input_upper) >= 2]
    else:
        ft = list(enr_rows)

    if not ft:
        print(f"  SKIP {fname}: no matching enrichment terms")
        return False

    ft.sort(key=lambda r: float(r['Adjusted P-value']))

    # Greedy set-cover
    covered, selected = set(), []
    remaining = list(range(len(ft)))

    while len(selected) < MAX_TERMS and remaining:
        best = max(remaining, key=lambda i: (
            len(set(ft[i]['Genes'].upper().split(';'))
                & input_upper - covered),
            float(ft[i].get('Combined Score', 0))))
        tg = set(ft[best]['Genes'].upper().split(';')) & input_upper
        new_genes = tg - covered

        if not new_genes and len(selected) >= 8:
            # Fill remaining slots by Combined Score
            rest = sorted(remaining,
                          key=lambda i: float(ft[i].get('Combined Score', 0)),
                          reverse=True)
            for idx in rest:
                if len(selected) >= MAX_TERMS:
                    break
                selected.append(idx)
                covered.update(
                    set(ft[idx]['Genes'].upper().split(';')) & input_upper)
            break

        selected.append(best)
        covered.update(tg)
        remaining.remove(best)

    if not selected:
        print(f"  SKIP {fname}: no terms after set-cover")
        return False

    # Build NetworkX graph 
    G = nx.Graph()
    for idx in selected:
        r = ft[idx]
        t = clean_term(r['Term'])
        G.add_node(t, ntype='term', lib=r['Gene_set'])
        tg = set(r['Genes'].upper().split(';'))
        top_genes = sorted(tg & input_upper,
                           key=lambda g: -abs(lfc_map.get(g, 0)))
        for gu in top_genes[:MAX_GENES_PER_TERM]:
            gn = cons_map.get(gu, gu)
            if gn not in G:
                G.add_node(gn, ntype='gene', lfc=lfc_map.get(gu, 0))
            G.add_edge(t, gn)

    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
    term_nodes = [n for n, d in G.nodes(data=True) if d.get('ntype') == 'term']
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get('ntype') == 'gene']

    if len(gene_nodes) < 2:
        print(f"  SKIP {fname}: fewer than 2 gene nodes")
        return False

    n_edges = G.number_of_edges()
    print(f"  {fname}: {len(term_nodes)} terms, "
          f"{len(gene_nodes)} genes, {n_edges} edges")

    # ── Spring layout ────────────────────────────────────────
    pos = nx.spring_layout(G, k=K_SPRING, iterations=ITERATIONS, seed=SEED)

    # ── Figure ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=FIGSIZE)
    fig.patch.set_alpha(0)
    ax.set_facecolor('none')

    # Edges
    nx.draw_networkx_edges(G, pos, ax=ax,
                           alpha=EDGE_ALPHA, width=EDGE_WIDTH,
                           edge_color=EDGE_COLOR)

    # Gene nodes (scatter, coloured by LFC)
    gene_lfcs = [G.nodes[g].get('lfc', 0) for g in gene_nodes]
    gene_sizes = [max(GENE_MIN_SIZE, G.degree(g) * GENE_DEGREE_SCALE)
                  for g in gene_nodes]
    gene_xy = np.array([pos[g] for g in gene_nodes])

    sc = ax.scatter(
        gene_xy[:, 0], gene_xy[:, 1],
        c=[np.clip(l, -VMAX, VMAX) for l in gene_lfcs],
        cmap='RdBu_r', vmin=-VMAX, vmax=VMAX,
        s=gene_sizes, zorder=3,
        edgecolors=GENE_EDGE_COLOR, linewidths=GENE_EDGE_WIDTH)

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02)
    cbar.set_label('mean log₂FC', fontsize=20)
    cbar.ax.tick_params(labelsize=14)

    # Term nodes (coloured squares)
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get('lib', ''), '#888888')
        xy = pos[t]
        ax.scatter(xy[0], xy[1], s=TERM_SQUARE_SIZE, c=color,
                   zorder=4, edgecolors=TERM_EDGE_COLOR,
                   linewidths=TERM_EDGE_WIDTH, marker=TERM_MARKER)

    # Term labels (white on coloured box, ABOVE node)
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get('lib', ''), '#888888')
        xy = pos[t]
        ax.text(xy[0], xy[1] + TERM_LABEL_OFFSET,
                wrap_label(t, TERM_WRAP_WIDTH),
                fontsize=TERM_LABEL_FONTSIZE, fontweight='bold',
                color='white', ha='center', va='bottom', zorder=5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=color,
                          alpha=TERM_BOX_ALPHA, edgecolor='none'))

    # Gene labels (white on dark box, BELOW node)
    for g in gene_nodes:
        xy = pos[g]
        ax.text(xy[0], xy[1] - GENE_LABEL_OFFSET, g,
                fontsize=GENE_LABEL_FONTSIZE, fontweight='bold',
                color='white', ha='center', va='top', zorder=5,
                bbox=dict(boxstyle='round,pad=0.2',
                          facecolor=GENE_BOX_COLOR,
                          alpha=GENE_BOX_ALPHA, edgecolor='none'))

    # Legend
    libs_used = {G.nodes[t].get('lib') for t in term_nodes}
    legend_els = [
        Line2D([0], [0], marker='s', color='w',
               markerfacecolor=c, markersize=12,
               label=LIB_LABELS.get(lib, lib))
        for lib, c in LIB_COLORS.items() if lib in libs_used
    ]
    if legend_els:
        ax.legend(handles=legend_els, loc='lower left',
                  fontsize=20, framealpha=0.8, facecolor='white')

    ax.set_title(title, fontsize=24, fontweight='bold', pad=20)
    ax.axis('off')

    # ── Save ─────────────────────────────────────────────────
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f'{fname}.png'),
                dpi=DPI, bbox_inches='tight', transparent=True)
    fig.savefig(os.path.join(out_dir, f'{fname}.pdf'),
                bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f"    saved {fname}.png + .pdf")
    return True
