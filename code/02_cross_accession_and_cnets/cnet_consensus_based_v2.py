import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import textwrap, re, os

DATA = '/sessions/practical-ecstatic-mendel/mnt/outputs'
OUT_CONS = '/sessions/practical-ecstatic-mendel/mnt/Bulk RNA sequencing/Final analysis/Plots/Cnet_Plots_Consensus_Based'

LIB_COLORS = {
    'GO_Biological_Process_2023': '#2196F3',
    'GO_Molecular_Function_2023': '#9C27B0',
    'GO_Cellular_Component_2023': '#FF9800',
    'KEGG_2021_Human':            '#4CAF50',
    'Reactome_2022':              '#F44336',
}
LIB_LABELS = {
    'GO_Biological_Process_2023': 'GO Biological Process',
    'GO_Molecular_Function_2023': 'GO Molecular Function',
    'GO_Cellular_Component_2023': 'GO Cellular Component',
    'KEGG_2021_Human':            'KEGG',
    'Reactome_2022':              'Reactome',
}

lfc_df = pd.read_csv(f'{DATA}/final_gene_lfc_COMPREHENSIVE.csv')
lfc_upper = {k.upper(): v for k, v in zip(lfc_df['symbol'], lfc_df['mean_lfc'])}
symbol_upper = {g.upper(): g for g in lfc_df['symbol']}

def clean_term(t):
    t = re.sub(r'\s*\(GO:\d+\)', '', t)
    t = re.sub(r'\s*R-HSA-\d+', '', t)
    return t[:52] + '...' if len(t) > 55 else t

def wrap_label(t, width=24):
    return '\n'.join(textwrap.wrap(t, width=width))

def build_cnet(enr_sig, input_genes_upper, title, outpath, max_terms=15):
    sig = enr_sig.copy()
    
    # Pre-process
    term_data = []
    for _, row in sig.iterrows():
        genes_raw = [g.strip() for g in str(row['Genes']).split(';') if g.strip()]
        genes_in = [g.upper() for g in genes_raw if g.upper() in input_genes_upper]
        if genes_in:
            term_data.append({
                'term': clean_term(row['Term']),
                'lib': row.get('Gene_set', 'Unknown'),
                'combined_score': row.get('Combined Score', 0),
                'input_genes': genes_in,
            })
    
    if len(term_data) < 2:
        print(f"  SKIPPED: too few terms with gene matches")
        return None
    
    # Greedy set-cover
    selected, covered, remaining = [], set(), list(range(len(term_data)))
    while len(selected) < max_terms and remaining:
        best = max(remaining, key=lambda i: (
            len(set(term_data[i]['input_genes']) - covered),
            term_data[i]['combined_score']))
        new_genes = set(term_data[best]['input_genes']) - covered
        if not new_genes and len(selected) >= 8:
            rest = sorted(remaining, key=lambda i: term_data[i]['combined_score'], reverse=True)
            for idx in rest:
                if len(selected) >= max_terms: break
                if term_data[idx]['term'] not in {s['term'] for s in selected}:
                    selected.append(term_data[idx])
                    covered.update(term_data[idx]['input_genes'])
            break
        selected.append(term_data[best])
        covered.update(term_data[best]['input_genes'])
        remaining.remove(best)
    
    print(f"  Selected {len(selected)} terms, {len(covered)} genes covered")
    for s in selected:
        print(f"    [{LIB_LABELS.get(s['lib'], s['lib'])[:10]}] {s['term']}")
    
    # Build graph
    G = nx.Graph()
    for td in selected:
        G.add_node(td['term'], ntype='term', lib=td['lib'])
        for g in td['input_genes']:
            G.add_node(g, ntype='gene', lfc=lfc_upper.get(g, 0))
            G.add_edge(td['term'], g)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
    
    term_nodes = [n for n, d in G.nodes(data=True) if d.get('ntype') == 'term']
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get('ntype') == 'gene']
    
    # Layout
    pos = nx.spring_layout(G, k=3.5, iterations=150, seed=42)
    fig, ax = plt.subplots(figsize=(30, 26))
    fig.patch.set_alpha(0); ax.set_facecolor('none')
    
    # Edges
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2, width=0.7, edge_color='#888888')
    
    # Gene nodes
    gene_lfcs = [G.nodes[g].get('lfc', 0) for g in gene_nodes]
    gene_sizes = [max(80, G.degree(g) * 40) for g in gene_nodes]
    vmax = max(abs(min(gene_lfcs)), abs(max(gene_lfcs)), 0.5)
    gene_xy = np.array([pos[g] for g in gene_nodes])
    sc = ax.scatter(gene_xy[:, 0], gene_xy[:, 1], c=gene_lfcs,
                    cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                    s=gene_sizes, zorder=3, edgecolors='white', linewidths=0.5)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02)
    cbar.set_label('mean log₂FC', fontsize=20)
    
    # Term nodes
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get('lib'), '#888888')
        xy = pos[t]
        ax.scatter(xy[0], xy[1], s=350, c=color, zorder=4,
                   edgecolors='white', linewidths=1.5, marker='s')
    
    # Term labels (white on colored box)
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get('lib'), '#888888')
        xy = pos[t]
        ax.text(xy[0], xy[1] + 0.035, wrap_label(t, 24),
                fontsize=20, fontweight='bold', color='white',
                ha='center', va='bottom', zorder=5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=color,
                          alpha=0.9, edgecolor='none'))
    
    # Gene labels (white on dark box)
    for g in gene_nodes:
        xy = pos[g]
        ax.text(xy[0], xy[1] - 0.03, g,
                fontsize=20, fontweight='bold', color='white',
                ha='center', va='top', zorder=5,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='#333333',
                          alpha=0.85, edgecolor='none'))
    
    # Legend
    legend_els = []
    libs_used = {G.nodes[t].get('lib') for t in term_nodes}
    for lib, color in LIB_COLORS.items():
        if lib in libs_used:
            legend_els.append(
                Line2D([0], [0], marker='s', color='w',
                       markerfacecolor=color, markersize=12,
                       label=LIB_LABELS.get(lib, lib)))
    if legend_els:
        ax.legend(handles=legend_els, loc='lower left', fontsize=20,
                  framealpha=0.8, facecolor='white')
    
    ax.set_title(title, fontsize=24, fontweight='bold', pad=20)
    ax.axis('off')
    
    fig.savefig(outpath, dpi=200, bbox_inches='tight', transparent=True)
    plt.close(fig)
    non_zero = sum(1 for v in gene_lfcs if abs(v) > 0.001)
    print(f"  SAVED: {os.path.basename(outpath)} ({non_zero}/{len(gene_nodes)} genes with LFC≠0)")


# ============================================================
# ALL CONSENSUS AT LFC≥1
# ============================================================
print("=" * 60)
print("CONSENSUS-BASED CNET (LFC≥1)")
print("=" * 60)

cons1 = pd.read_csv(f'{DATA}/final_consensus_LFC1.csv')
genes1 = cons1['symbol'].tolist()
input_upper = {g.upper() for g in genes1}
print(f"Consensus genes at |LFC|≥1: {len(genes1)}")

enr1 = pd.read_csv(f'{DATA}/final_enrichment_LFC1.csv')
# Use nominal p<0.01 since few reach padj<0.05 with 491 genes
sig1 = enr1[enr1['P-value'] < 0.01].copy()
print(f"Terms at p<0.01: {len(sig1)}")

build_cnet(sig1, input_upper,
           f'Concept Network – All Consensus DEGs (|log₂FC| ≥ 1.0)\n{len(genes1)} genes × All Databases',
           f'{OUT_CONS}/Cnet_All_Consensus_LFC1.png',
           max_terms=15)

# BBB-related at LFC≥1
print("\n--- BBB-related consensus LFC≥1 ---")
bbb_kw = r'(?i)blood.brain|barrier|tight.junction|endotheli|vascular|permeab|transport|efflux|cell.junction|angiogen|MMP|extracellular.matrix|cell.migrat'
bbb_sig = sig1[sig1['Term'].str.contains(bbb_kw, na=False)]
print(f"BBB terms: {len(bbb_sig)}")
if len(bbb_sig) >= 2:
    build_cnet(bbb_sig, input_upper,
               f'BBB-Related Consensus DEGs (|log₂FC| ≥ 1.0)\nConcept Network',
               f'{OUT_CONS}/Cnet_BBB_Consensus_LFC1.png',
               max_terms=12)

# Inflammation+Survival at LFC≥1
print("\n--- Inflammation+Survival consensus LFC≥1 ---")
is_kw = r'(?i)inflam|immune|cytokine|interleukin|NF.?kB|apoptosis|survival|cell.death|necrosis|caspase|TNF|toll|TLR|innate|chemokine'
is_sig = sig1[sig1['Term'].str.contains(is_kw, na=False)]
print(f"InflamSurv terms: {len(is_sig)}")
if len(is_sig) >= 2:
    build_cnet(is_sig, input_upper,
               f'Inflammation & Survival Consensus (|log₂FC| ≥ 1.0)\nConcept Network',
               f'{OUT_CONS}/Cnet_InflamSurv_Consensus_LFC1.png',
               max_terms=12)

print("\n\nAll done!")
