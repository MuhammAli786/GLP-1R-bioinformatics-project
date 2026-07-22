import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import os

lfc_all = pd.read_csv('final_gene_lfc_COMPREHENSIVE.csv').set_index('symbol')['mean_lfc'].to_dict()
lfc_lower = {k.lower(): v for k, v in lfc_all.items()}

def get_lfc(gene):
    v = lfc_all.get(gene)
    if v is not None: return v
    v = lfc_lower.get(gene.lower())
    if v is not None: return v
    return 0.0

def diverse_greedy_select(sig_df, input_genes, max_terms=10, max_per_db=3):
    cons_lower = {g.lower(): g for g in input_genes}
    sig = sig_df.copy()
    sig['parsed'] = sig['Genes'].apply(lambda x: [g.strip() for g in str(x).split(';') if g.strip()])
    sig['overlap'] = sig['parsed'].apply(lambda gs: set(cons_lower[g.lower()] for g in gs if g.lower() in cons_lower))
    sig['n_overlap'] = sig['overlap'].apply(len)
    sig = sig[sig['n_overlap'] >= 3].copy()
    if sig.empty: return []
    covered = set()
    selected = []
    db_counts = {}
    remaining = sig.copy()
    for _ in range(max_terms):
        if remaining.empty: break
        remaining['new_genes'] = remaining['overlap'].apply(lambda gs: len(gs - covered))
        eligible = remaining[remaining['Gene_set'].apply(lambda db: db_counts.get(db, 0) < max_per_db)]
        if eligible.empty: eligible = remaining
        eligible = eligible[eligible['new_genes'] > 0]
        if eligible.empty: break
        best_idx = eligible['new_genes'].idxmax()
        best = remaining.loc[best_idx]
        selected.append({'term': best['Term'], 'genes': best['overlap'], 'db': best['Gene_set']})
        covered |= best['overlap']
        db_counts[best['Gene_set']] = db_counts.get(best['Gene_set'], 0) + 1
        remaining = remaining.drop(best_idx)
    return selected

def build_cnet(selected_terms, title, outpath):
    G = nx.Graph()
    for item in selected_terms:
        term = item['term']
        genes = item['genes']
        short = term.split('(GO:')[0].split('R-HSA')[0].strip()
        if len(short) > 40: short = short[:37] + '...'
        G.add_node(short, node_type='term')
        for g in genes:
            G.add_node(g, node_type='gene')
            G.add_edge(short, g)
    if len(G.nodes()) == 0:
        print(f"  No nodes, skipping"); return
    terms = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'term']
    genes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
    gene_lfc = {g: get_lfc(g) for g in genes}
    pos = nx.spring_layout(G, k=4.0, iterations=80, seed=42)
    fig, ax = plt.subplots(figsize=(30, 26))
    fig.patch.set_alpha(0); ax.set_facecolor('none')
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.35, width=2.5, edge_color='gray')
    nx.draw_networkx_nodes(G, pos, nodelist=terms, ax=ax, node_color='lightgray', node_size=3000, edgecolors='black', linewidths=2)
    lfc_vals = [gene_lfc[g] for g in genes]
    norm = Normalize(vmin=-2, vmax=2)
    cmap = plt.cm.RdBu_r
    colors = [cmap(norm(v)) for v in lfc_vals]
    nx.draw_networkx_nodes(G, pos, nodelist=genes, ax=ax, node_color=colors, node_size=800, edgecolors='black', linewidths=1)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in terms}, ax=ax, font_size=20, font_weight='bold')
    nx.draw_networkx_labels(G, pos, labels={n: n for n in genes}, ax=ax, font_size=12)
    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label('Mean log2FC', fontsize=16)
    ax.set_title(title, fontsize=24, fontweight='bold'); ax.axis('off')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight', transparent=True)
    plt.close()
    non_zero = sum(1 for v in lfc_vals if abs(v) > 0.001)
    print(f"  Saved {outpath} ({len(terms)} terms, {len(genes)} genes, {non_zero}/{len(genes)} LFC≠0)")

# ============================================================
# CONSENSUS-BASED BBB & InflamSurv CNETS
# ============================================================
print("=== CONSENSUS-BASED PATHWAY CNETS ===")

# BBB-related keywords
bbb_terms = ['blood-brain','barrier','tight junction','claudin','occludin','endotheli',
             'vascular','permeab','transport','efflux','abc trans']
inflam_surv_terms = ['inflam','immune','cytokine','interleukin','nf-kb','apoptosis',
                     'survival','cell death','necrosis','caspase','bcl','tnf']

for thresh in ['LFC0', 'LFC05']:
    cons = pd.read_csv(f'final_consensus_{thresh}.csv')
    genes = cons['symbol'].tolist()
    enrich = pd.read_csv(f'final_enrichment_{thresh}.csv')
    sig = enrich[enrich['Adjusted P-value'] < 0.05].copy()
    
    # BBB consensus cnet
    bbb_sig = sig[sig['Term'].str.lower().apply(lambda t: any(kw in t for kw in bbb_terms))]
    if len(bbb_sig) >= 1:
        selected = diverse_greedy_select(bbb_sig, genes, max_terms=8, max_per_db=3)
        if selected:
            print(f"  BBB consensus {thresh}: {len(selected)} terms")
            build_cnet(selected, f'BBB-Related Consensus Genes ({thresh.replace("LFC","LFC≥")})\nAll Databases',
                       f'final_plots/cnet_NEW_BBB_consensus_{thresh}.png')
    
    # Inflammation+Survival consensus cnet
    is_sig = sig[sig['Term'].str.lower().apply(lambda t: any(kw in t for kw in inflam_surv_terms))]
    if len(is_sig) >= 1:
        selected = diverse_greedy_select(is_sig, genes, max_terms=8, max_per_db=3)
        if selected:
            print(f"  InflamSurv consensus {thresh}: {len(selected)} terms")
            build_cnet(selected, f'Inflammation & Survival Consensus ({thresh.replace("LFC","LFC≥")})\nAll Databases',
                       f'final_plots/cnet_NEW_InflamSurv_consensus_{thresh}.png')

# ============================================================
# LFC≥1 DOTPLOTS
# ============================================================
print("\n=== LFC≥1 DOTPLOTS ===")
enrich1 = pd.read_csv('final_enrichment_LFC1.csv')

for db in enrich1['Gene_set'].unique():
    subset = enrich1[enrich1['Gene_set'] == db].copy()
    # Use top 15 by P-value regardless of significance
    top = subset.nsmallest(15, 'P-value')
    if len(top) == 0: continue
    
    fig, ax = plt.subplots(figsize=(12, max(6, len(top)*0.5)))
    fig.patch.set_alpha(0); ax.set_facecolor('none')
    
    top['neg_log_p'] = -np.log10(top['Adjusted P-value'].clip(lower=1e-50))
    top['gene_count'] = top['Genes'].apply(lambda x: len(str(x).split(';')))
    top = top.sort_values('neg_log_p')
    
    # Shorten term names
    def shorten(t):
        t = t.split('(GO:')[0].split('R-HSA')[0].strip()
        return t[:60] + '...' if len(t) > 60 else t
    top['short_term'] = top['Term'].apply(shorten)
    
    scatter = ax.scatter(top['neg_log_p'], range(len(top)), 
                        s=top['gene_count']*15, c=top['neg_log_p'],
                        cmap='RdYlBu_r', edgecolors='black', linewidths=0.5, alpha=0.8)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['short_term'], fontsize=10)
    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=12)
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='padj=0.05')
    ax.legend(fontsize=10)
    ax.set_title(f'{db} (LFC≥1.0 Consensus)', fontsize=14, fontweight='bold')
    
    plt.colorbar(scatter, ax=ax, label='-log10(padj)', fraction=0.03)
    plt.tight_layout()
    plt.savefig(f'final_plots/dotplot_{db}_LFC1.png', dpi=150, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"  Saved dotplot_{db}_LFC1.png")

print("\nAll done!")
