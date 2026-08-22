# Rebuilds the consensus cnet plots with a per-database quota so terms are drawn
# from several enrichment libraries rather than one.
# final_consensus_*.csv + final_enrichment_*.csv -> final_plots/cnet_all_*.png

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import os

# Mean log2FC per gene, from the comprehensive lookup table.
lfc_all = pd.read_csv('final_gene_lfc_COMPREHENSIVE.csv').set_index('symbol')['mean_lfc'].to_dict()
lfc_lower = {k.lower(): v for k, v in lfc_all.items()}
print(f"LFC lookup: {len(lfc_all)} genes")

# Look up a gene's mean log2FC, falling back to a case-insensitive match, then 0.
def get_lfc(gene):
    v = lfc_all.get(gene)
    if v is not None: return v
    v = lfc_lower.get(gene.lower())
    if v is not None: return v
    return 0.0

def diverse_greedy_select(sig_df, input_genes, max_terms=10, max_per_db=4):
    """Greedy set cover over enrichment terms, capped at max_per_db terms per database."""
    cons_lower = {g.lower(): g for g in input_genes}
    
    sig = sig_df.copy()
    sig['parsed'] = sig['Genes'].apply(lambda x: [g.strip() for g in str(x).split(';') if g.strip()])
    sig['overlap'] = sig['parsed'].apply(lambda gs: set(cons_lower[g.lower()] for g in gs if g.lower() in cons_lower))
    sig['n_overlap'] = sig['overlap'].apply(len)
    sig = sig[sig['n_overlap'] >= 3].copy()
    
    if sig.empty:
        return []
    
    covered = set()
    selected = []
    db_counts = {}
    remaining = sig.copy()
    
    for _ in range(max_terms):
        if remaining.empty:
            break
        remaining['new_genes'] = remaining['overlap'].apply(lambda gs: len(gs - covered))
        # Drop databases that have hit their quota
        eligible = remaining[remaining['Gene_set'].apply(lambda db: db_counts.get(db, 0) < max_per_db)]
        if eligible.empty:
            eligible = remaining  # every database is at quota; ignore it
        
        eligible = eligible[eligible['new_genes'] > 0]
        if eligible.empty:
            break
            
        best_idx = eligible['new_genes'].idxmax()
        best = remaining.loc[best_idx]
        
        selected.append({'term': best['Term'], 'genes': best['overlap'], 'db': best['Gene_set']})
        covered |= best['overlap']
        db_counts[best['Gene_set']] = db_counts.get(best['Gene_set'], 0) + 1
        remaining = remaining.drop(best_idx)
    
    return selected

def build_cnet(selected_terms, title, outpath, input_set=None):
    """Render one cnet figure from a set-cover term selection, colouring genes by mean log2FC."""
    G = nx.Graph()
    
    for item in selected_terms:
        term = item['term']
        genes = item['genes']
        short = term.split('(GO:')[0].split('R-HSA')[0].strip()
        if len(short) > 40:
            short = short[:37] + '...'
        G.add_node(short, node_type='term')
        for g in genes:
            G.add_node(g, node_type='gene')
            G.add_edge(short, g)
    
    if len(G.nodes()) == 0:
        print(f"  No nodes for {title}, skipping")
        return
    
    terms = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'term']
    genes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
    
    gene_lfc = {g: get_lfc(g) for g in genes}
    
    pos = nx.spring_layout(G, k=4.0, iterations=80, seed=42)
    
    fig, ax = plt.subplots(figsize=(30, 26))
    fig.patch.set_alpha(0)
    ax.set_facecolor('none')
    
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.35, width=2.5, edge_color='gray')
    
    # Term nodes, coloured by source library
    nx.draw_networkx_nodes(G, pos, nodelist=terms, ax=ax, 
                           node_color='lightgray', node_size=3000, 
                           edgecolors='black', linewidths=2)
    
    # Gene nodes, coloured by mean log2FC
    lfc_vals = [gene_lfc[g] for g in genes]
    vmin, vmax = -2, 2
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.RdBu_r
    colors = [cmap(norm(v)) for v in lfc_vals]
    
    nx.draw_networkx_nodes(G, pos, nodelist=genes, ax=ax,
                           node_color=colors, node_size=800,
                           edgecolors='black', linewidths=1)
    
    nx.draw_networkx_labels(G, pos, labels={n: n for n in terms}, ax=ax,
                            font_size=20, font_weight='bold')
    nx.draw_networkx_labels(G, pos, labels={n: n for n in genes}, ax=ax,
                            font_size=12)
    
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label('Mean log2FC', fontsize=16)
    
    ax.set_title(title, fontsize=24, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight', transparent=True)
    plt.close()
    
    non_zero = sum(1 for v in lfc_vals if abs(v) > 0.001)
    print(f"  Saved {outpath} ({len(terms)} terms, {len(genes)} genes, {non_zero}/{len(genes)} with LFC≠0)")

os.makedirs('final_plots', exist_ok=True)

# Consensus-based cnets: consensus genes -> enrichment -> database-diverse term selection.
print("\n=== CONSENSUS-BASED CNETS ===")
for thresh in ['LFC0', 'LFC05']:
    cons = pd.read_csv(f'final_consensus_{thresh}.csv')
    genes = cons['symbol'].tolist()
    enrich = pd.read_csv(f'final_enrichment_{thresh}.csv')
    sig = enrich[enrich['Adjusted P-value'] < 0.05]
    
    print(f"\n{thresh}: {len(genes)} consensus genes, {len(sig)} sig terms")
    selected = diverse_greedy_select(sig, genes, max_terms=10, max_per_db=3)
    
    dbs_used = set(s['db'] for s in selected)
    print(f"  Databases used: {dbs_used}")
    
    build_cnet(selected, f'Consensus Gene Network ({thresh.replace("LFC","LFC≥")})\nAll Databases',
               f'final_plots/cnet_all_{thresh}.png', set(genes))

print("\nDone with consensus cnets!")
