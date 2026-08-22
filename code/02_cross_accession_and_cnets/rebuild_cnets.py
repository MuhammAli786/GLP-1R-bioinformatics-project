"""
Rebuilds the cnet plots, colouring genes by the mean log2FC lookup derived from the
full master DEG table rather than from the consensus table.
final_gene_lfc_lookup*.csv + final_enrichment_*.csv -> final_plots/cnet_*.png
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import os

OUT = BASE + "/mnt/outputs"
PLOT_DIR = os.path.join(OUT, "final_plots")

# Mean log2FC per gene, averaged over all master records rather than consensus genes only.
lfc_lookup_05 = pd.read_csv(os.path.join(OUT, 'final_gene_lfc_lookup.csv')).set_index('symbol')['mean_lfc'].to_dict()
lfc_lookup_0 = pd.read_csv(os.path.join(OUT, 'final_gene_lfc_lookup_LFC0.csv')).set_index('symbol')['mean_lfc'].to_dict()

enr05 = pd.read_csv(os.path.join(OUT, 'final_enrichment_LFC05.csv'))
enr0 = pd.read_csv(os.path.join(OUT, 'final_enrichment_LFC0.csv'))
for df in [enr05, enr0]:
    if 'Combined Score' in df.columns and 'Combined_Score' not in df.columns:
        df['Combined_Score'] = df['Combined Score']

# Pathway gene lists
JAK_GENES = ['Gfap','Il1r1','Il6ra','Il6st','Lifr','Osmr','Mapk1','Fyn','Akt1','Jak1','Nfkbia','Nfkb2','Ikbkg','Pias2','Pias1','Ptpn11','Ptprd','Ptprt','Egfr','Erbb2','Fgfr1','Socs5','Socs3','Socs4','Socs6','Stat1','Cdkn1a','Bcl2l1','Vegfa','Hif1a','Mmp2','Mcl1','Ccnd1']
BBB_GENES = ['Adam10','Adam9','Adam15','Cdh5','Ctnnb1','Ctnnd1','Vegfa','Hif1a','Aqp4','Gfap','Slc2a1','Lrp1','Slc7a5','Slc16a1','Tfrc','Abcb1b','Nid1','Lamc1','Fn1','Spp1','Kdr','Cav1','Vwf','Pecam1','Eng','Nos3','Mmp2','Mmp16','Mmp14','Mmp25','Nfkbia','Nfkb2','Mapk1','Rock2','Fyn','Rock1','Akt1','Cldn5','Tjp1','Jam2','F11r','Cldn12','Ocln','Tjp2']
AKT_GENES = ['Irs1','Irs2','Sos1','Gab1','Grb2','Pdpk1','Akt3','Akt1','Bcl2l1','Bad','Mcl1','Bax','Xiap','Birc2','Ulk2','Map1lc3a','Becn1','Ulk1','Creb1','Creb5','Creb3','Cdkn1a','Ccnd1','Foxo1','Foxo3','Foxo4','Foxo6','Gsk3a','Gsk3b','Vegfa','Kdr','Insr','Igf1r','Egfr','Erbb2','Erbb3','Fgfr1','Fgfr2','Glp1r','Nfkbia','Nfkb2','Phlpp1','Phlpp2','Pik3ca','Pik3r1','Pik3r3','Pik3c2a','Pik3cb','Pik3cg','Pik3c2b','Pten','Hras','Kras','Mapk1','Map2k2','Nras','Raf1','Fyn','Mtor','Tsc2','Rptor','Rictor','Rheb','Tsc1','Eif4g1','Eif4ebp1','Rps6kb1','Rps6kb2']
COMBINED_GENES = list(set(JAK_GENES + AKT_GENES))

# Render one cnet figure: greedy set-cover term selection, then a spring-layout network.
def build_cnet(enr_df, lfc_dict, title, save_name, max_terms=15, input_set=None):
    sig = enr_df[enr_df['Adjusted P-value'] < 0.05].copy()
    if sig.empty: print(f"  SKIP {save_name}: no sig terms"); return
    score_col = 'Combined_Score'
    sig['gene_list'] = sig['Genes'].apply(lambda x: [g.strip() for g in str(x).split(';') if g.strip()])
    
    if input_set:
        input_lower = set(g.lower() for g in input_set)
        sig['relevant_genes'] = sig['gene_list'].apply(lambda gl: [g for g in gl if g.lower() in input_lower or g in input_set])
        sig = sig[sig['relevant_genes'].apply(len) > 0]
        if sig.empty: print(f"  SKIP {save_name}: no terms with pathway genes"); return
    else:
        sig['relevant_genes'] = sig['gene_list']
    
    # Greedy set-cover
    selected_terms = []; covered_genes = set(); remaining = sig.copy()
    for _ in range(max_terms):
        if remaining.empty: break
        remaining['new_genes'] = remaining['relevant_genes'].apply(lambda gl: len(set(g.lower() for g in gl) - set(g.lower() for g in covered_genes)))
        remaining = remaining.sort_values(['new_genes', score_col], ascending=[False, False])
        best = remaining.iloc[0]
        if best['new_genes'] == 0:
            rest = remaining.nlargest(max_terms - len(selected_terms), score_col)
            for _, row in rest.iterrows(): selected_terms.append(row)
            break
        selected_terms.append(best)
        covered_genes.update(g.lower() for g in best['relevant_genes'])
        remaining = remaining.iloc[1:]
    
    if not selected_terms: print(f"  SKIP {save_name}"); return
    
    # Build the term-gene graph
    G = nx.Graph()
    for term_row in selected_terms:
        term_name = term_row['Term']; genes = term_row['relevant_genes']
        G.add_node(term_name, node_type='term')
        for g in genes: G.add_node(g, node_type='gene'); G.add_edge(term_name, g)
    
    if len(G.nodes) < 3: print(f"  SKIP {save_name}: too few nodes"); return
    
    gene_nodes = [n for n in G.nodes if G.nodes[n].get('node_type')=='gene']
    term_nodes = [n for n in G.nodes if G.nodes[n].get('node_type')=='term']
    
    # Colour genes from the master log2FC lookup, symmetric about zero.
    lfc_vals = [lfc_dict.get(g, 0) for g in gene_nodes]
    vmax = max(abs(v) for v in lfc_vals) if lfc_vals and max(abs(v) for v in lfc_vals) > 0 else 1
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap = plt.cm.RdBu_r
    gene_colors = [cmap(norm(lfc_dict.get(g, 0))) for g in gene_nodes]
    
    pos = nx.spring_layout(G, k=4.0, iterations=150, seed=42)
    fig, ax = plt.subplots(figsize=(30, 26)); ax.set_facecolor('none'); fig.patch.set_alpha(0)
    nx.draw_networkx_edges(G, pos, ax=ax, width=2.5, alpha=0.35, edge_color='#666666')
    
    gene_sizes = [max(120, G.degree(g)*50) for g in gene_nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_color=gene_colors, node_size=gene_sizes, ax=ax, edgecolors='black', linewidths=0.5)
    nx.draw_networkx_nodes(G, pos, nodelist=term_nodes, node_color='none', node_size=1, ax=ax)
    
    colors_palette = plt.cm.Set3(np.linspace(0,1,max(len(term_nodes),12)))
    for i, term in enumerate(term_nodes):
        x,y = pos[term]
        ax.text(x,y,term, fontsize=20, fontweight='bold', ha='center', va='center',
               bbox=dict(boxstyle='round,pad=0.3', facecolor=colors_palette[i%len(colors_palette)], alpha=0.85, edgecolor='none'), color='white', zorder=5)
    for g in gene_nodes:
        x,y = pos[g]
        lfc_val = lfc_dict.get(g, 0)
        ax.text(x,y-0.03,f"{g}\n({lfc_val:+.2f})", fontsize=18, fontweight='bold', ha='center', va='top',
               bbox=dict(boxstyle='round,pad=0.15', facecolor='#333333', alpha=0.7, edgecolor='none'), color='white', zorder=4)
    
    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.3, pad=0.02); cbar.set_label('Mean log2FC', fontsize=20); cbar.ax.tick_params(labelsize=16)
    ax.set_title(title, fontsize=24, fontweight='bold', pad=20); ax.axis('off')
    ax.text(0.02, 0.02, f"Genes: {len(gene_nodes)} | Terms: {len(term_nodes)}", transform=ax.transAxes, fontsize=20, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, save_name), dpi=150, transparent=True, bbox_inches='tight')
    plt.close()
    print(f"  Created: {save_name} ({len(gene_nodes)} genes, {len(term_nodes)} terms)")

print("=== Rebuilding ALL cnet plots with corrected LFC ===\n")

# |LFC| >= 0.5 consensus
print("--- LFC≥0.5 ---")
build_cnet(enr05, lfc_lookup_05, 'GLP-1R Agonist Consensus Genes\n(|LFC|≥0.5, ≥2 Groups)', 'cnet_all_LFC05.png')
build_cnet(enr05, lfc_lookup_05, 'JAK-STAT3/Inflammatory Pathway\n(|LFC|≥0.5)', 'cnet_JAKSTAT3_LFC05.png', input_set=JAK_GENES)
build_cnet(enr05, lfc_lookup_05, 'Blood-Brain Barrier / MMP Pathway\n(|LFC|≥0.5)', 'cnet_BBB_LFC05.png', input_set=BBB_GENES)
build_cnet(enr05, lfc_lookup_05, 'PI3K/Akt/mTOR Survival Pathway\n(|LFC|≥0.5)', 'cnet_AktPI3K_LFC05.png', input_set=AKT_GENES)
build_cnet(enr05, lfc_lookup_05, 'Combined JAK-STAT3 + PI3K/Akt Pathways\n(|LFC|≥0.5)', 'cnet_Combined_LFC05.png', input_set=COMBINED_GENES)

# |LFC| >= 0 consensus (padj < 0.05 only)
print("\n--- LFC≥0 ---")
build_cnet(enr0, lfc_lookup_0, 'GLP-1R Agonist Consensus Genes\n(padj<0.05 only, ≥2 Groups)', 'cnet_all_LFC0.png')
build_cnet(enr0, lfc_lookup_0, 'JAK-STAT3/Inflammatory Pathway\n(padj<0.05 only)', 'cnet_JAKSTAT3_LFC0.png', input_set=JAK_GENES)
build_cnet(enr0, lfc_lookup_0, 'Blood-Brain Barrier / MMP Pathway\n(padj<0.05 only)', 'cnet_BBB_LFC0.png', input_set=BBB_GENES)
build_cnet(enr0, lfc_lookup_0, 'PI3K/Akt/mTOR Survival Pathway\n(padj<0.05 only)', 'cnet_AktPI3K_LFC0.png', input_set=AKT_GENES)
build_cnet(enr0, lfc_lookup_0, 'Combined JAK-STAT3 + PI3K/Akt Pathways\n(padj<0.05 only)', 'cnet_Combined_LFC0.png', input_set=COMBINED_GENES)

print("\nDone!")
