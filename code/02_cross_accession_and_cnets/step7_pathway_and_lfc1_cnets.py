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
        print(f"  No nodes for {title}, skipping"); return
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
# LFC≥1 CNET
# ============================================================
print("=== LFC≥1 CNET ===")
cons1 = pd.read_csv('final_consensus_LFC1.csv')
genes1 = cons1['symbol'].tolist()
enrich1 = pd.read_csv('final_enrichment_LFC1.csv')
sig1 = enrich1[enrich1['Adjusted P-value'] < 0.05]
print(f"LFC≥1: {len(genes1)} genes, {len(sig1)} sig terms (all Reactome)")

# Since only Reactome is significant, use relaxed threshold for other DBs
# Use padj < 0.1 for non-Reactome to get some diversity
sig_relaxed = enrich1[(enrich1['Adjusted P-value'] < 0.1) | (enrich1['Gene_set'] == 'Reactome_2022') & (enrich1['Adjusted P-value'] < 0.05)]
sig_strict = enrich1[enrich1['Adjusted P-value'] < 0.05]
# Try with nominal p-value < 0.01 for diversity
sig_nominal = enrich1[enrich1['P-value'] < 0.01]
print(f"  Nominal p<0.01: {len(sig_nominal)} terms")
for db in sig_nominal['Gene_set'].unique():
    print(f"    {db}: {len(sig_nominal[sig_nominal['Gene_set']==db])}")

# Use nominal p<0.01 for diverse cnet, but note in title
selected1 = diverse_greedy_select(sig_nominal, genes1, max_terms=10, max_per_db=3)
dbs = set(s['db'] for s in selected1)
print(f"  Databases used: {dbs}")
build_cnet(selected1, 'Consensus Gene Network (LFC≥1.0)\nAll Databases (p<0.01)', 
           'final_plots/cnet_NEW_all_LFC1.png')

# ============================================================
# PATHWAY-BASED CNETS
# ============================================================
print("\n=== PATHWAY-BASED CNETS ===")

# Load pathway gene lists
pathway_genes = {
    'JAKSTAT3': ['Jak1','Jak2','Jak3','Tyk2','Stat1','Stat2','Stat3','Stat4','Stat5a','Stat5b','Stat6',
                 'Socs1','Socs2','Socs3','Socs4','Socs5','Socs6','Socs7','Cish','Pias1','Pias2','Pias3','Pias4',
                 'Il6','Il6ra','Il6st','Lif','Lifr','Cntf','Cntfr','Osm','Osmr','Il10','Il10ra','Il10rb',
                 'Il21','Il21r','Il23a','Il12rb1','Ifng','Ifngr1','Ifngr2','Epo','Epor','Ghr','Prl','Prlr',
                 'Csf2','Csf2ra','Csf2rb','Il2','Il2ra','Il2rb','Il2rg','Il7','Il7r','Il15','Il15ra',
                 'Myc','Ccnd1','Bcl2','Bcl2l1','Mcl1','Pim1','Pim2','Vegfa','Hif1a','Mmp2','Mmp9','Birc5'],
    'BBB': ['Cldn1','Cldn3','Cldn5','Cldn11','Cldn12','Ocln','Tjp1','Tjp2','Tjp3','F11r','Jam2','Jam3',
            'Cdh5','Pecam1','Esam','Nectin1','Nectin2','Nectin3','Pvrl4',
            'Abcb1a','Abcb1b','Abcg2','Abcc1','Abcc4','Slc2a1','Slc7a5','Slc16a1','Slc1a1',
            'Slc7a1','Slc3a2','Slc6a6','Mfsd2a',
            'Pdgfrb','Pdgfb','Tgfb1','Tgfbr1','Tgfbr2','Angpt1','Angpt2','Tek','Tie1',
            'Vegfa','Vegfr1','Vegfr2','Kdr','Flt1','Nrp1','Nrp2',
            'Aqp4','Gja1','Sox9','S100b','Gfap','Aldh1l1',
            'Mmp2','Mmp9','Mmp3','Timp1','Timp2','Timp3',
            'Icam1','Vcam1','Sele','Selp','Itgal','Itgb2','Cxcl12','Cxcr4'],
    'AktPI3K': ['Pik3ca','Pik3cb','Pik3cd','Pik3cg','Pik3r1','Pik3r2','Pik3r3',
                'Akt1','Akt2','Akt3','Pdk1','Pdpk1',
                'Mtor','Rptor','Rictor','Mlst8','Tsc1','Tsc2','Rheb',
                'Pten','Inpp5d','Inpp4b',
                'Gsk3a','Gsk3b','Foxo1','Foxo3','Foxo4',
                'Bad','Bcl2','Bcl2l1','Mcl1','Casp9','Casp3',
                'Cdkn1a','Cdkn1b','Ccnd1','Ccnd2','Ccne1','Cdk2','Cdk4','Cdk6','Rb1','Mdm2','Tp53',
                'Vegfa','Nos3','Hif1a','Rps6kb1','Eif4e','Eif4ebp1',
                'Egfr','Erbb2','Igf1r','Insr','Pdgfra','Pdgfrb','Fgfr1',
                'Grb2','Sos1','Hras','Kras','Nras','Raf1','Map2k1','Map2k2','Mapk1','Mapk3'],
}

# Combined = all pathway genes
combined_genes = set()
for v in pathway_genes.values():
    combined_genes.update(v)
pathway_genes['Combined'] = list(combined_genes)

for pname, pgenes in pathway_genes.items():
    for thresh in ['LFC0', 'LFC05']:
        enrich = pd.read_csv(f'final_enrichment_{thresh}.csv')
        sig = enrich[enrich['Adjusted P-value'] < 0.05]
        selected = diverse_greedy_select(sig, pgenes, max_terms=8, max_per_db=3)
        if not selected:
            print(f"  {pname} {thresh}: no terms found, skipping")
            continue
        dbs = set(s['db'] for s in selected)
        print(f"  {pname} {thresh}: {len(selected)} terms from {dbs}")
        build_cnet(selected, f'{pname} Pathway Network ({thresh.replace("LFC","LFC≥")})\nAll Databases',
                   f'final_plots/cnet_{pname}_{thresh}.png')

print("\nDone!")
