"""
Pathway-based concept network (cnet) plots for JAK-STAT3, BBB/MMP, and PI3K/Akt.
pathway_specific_enrichment.csv + final_gene_lfc_COMPREHENSIVE.csv -> Cnet_*_v2.png
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import textwrap, re, os

DATA = BASE + '/mnt/outputs'
OUT_PW = BASE + '/mnt/Bulk RNA sequencing/Final analysis/Plots/Cnet_Plots_Pathway_Based'
OUT_CONS = BASE + '/mnt/Bulk RNA sequencing/Final analysis/Plots/Cnet_Plots_Consensus_Based'

# Library colours
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

# Term and gene regex filters applied per pathway.
FILTERS = {
    'JAK_STAT3': {
        'term': r'(?i)JAK|STAT|interleukin|cytokine|inflam|interferon|NF.?kB|immune|innate|'
                r'toll.like|TLR|chemokine|TNF|apoptotic|cell.death|MAPK|signal.transduct|'
                r'growth.factor.recep|phosphatase|kinase.activ|receptor.signal|Src.family',
        'gene': r'JAK1|JAK2|STAT1|STAT3|IL6|IL6ST|NFKBIA|SOCS|MAPK1|FYN|BCL2L1|VEGFA|HIF1A|MMP2|PIAS|CDKN1A|MCL1',
    },
    'BBB_MMP': {
        'term': r'(?i)blood.brain|BBB|tight.junction|endotheli|barrier|claudin|occludin|'
                r'cell.junction|adherens|basement.membrane|extracellular.matrix|focal.adhesion|'
                r'leukocyte.transendotheli|VEGF|vascular|pericyte|astrocyte|transport.*barrier|'
                r'permeab|cell.cell|intercellular|junction|cell.adhesion|anchoring|gap.junction|'
                r'wound.heal|angiogen|MMP|matrix.metallo|collagen|laminin|transporter|'
                r'efflux|ABC.transporter|aquaporin|membrane.transport|caveol|'
                r'lipid.raft|cell.migrat|tube.morpho',
        'gene': r'MMP2|MMP9|MMP3|CLDN|OCLN|TJP1|TJP2|AQP4|ABCB|ABCG|SLC2A1|PECAM1|CDH5|'
                r'KDR|CAV1|F11R|JAM|ADAM|ROCK|LRP1|TFRC|SLC7A5|SLC16A1|HIF1A|'
                r'VEGFA|FYN|NFKBIA|ICAM|VCAM|PDGF|TGFB|ANGPT|NRP',
    },
    'Akt_PI3K': {
        'term': r'(?i)PI3K|AKT|mTOR|phosphoinositide|insulin|IGF|growth.factor|FOXO|GSK|'
                r'autophagy|apoptosis|pro.?surviv|MAPK|RAS|RAF|MEK|ERK|PTEN|phosphoryl|'
                r'kinase.signal|RTK|receptor.tyrosine|CREB|translation|EIF4|ribosom|cell.cycle|'
                r'proliferat|survival|ubiquitin|proteolys',
        'gene': r'AKT|PIK3|MTOR|PTEN|FOXO|GSK3|IRS1|IRS2|PDPK1|HRAS|KRAS|BAD|BCL2L1|MAPK1|'
                r'MAP2K|TSC2|EIF4|INSR|SOS1|CREB|RHEB|RPTOR|RICTOR|'
                r'NFKBIA|EGFR|ERBB|IGF1R|GRB2|NRAS|RAF1|CDK',
    },
}

lfc_df = pd.read_csv(f'{DATA}/final_gene_lfc_COMPREHENSIVE.csv')
lfc_dict = dict(zip(lfc_df['symbol'], lfc_df['mean_lfc']))
lfc_upper = {k.upper(): v for k, v in lfc_dict.items()}
symbol_upper = {g.upper(): g for g in lfc_dict.keys()}

# Strip GO / Reactome IDs from a term name and truncate long ones.
def clean_term(t):
    t = re.sub(r'\s*\(GO:\d+\)', '', t)
    t = re.sub(r'\s*R-HSA-\d+', '', t)
    return t[:52] + '...' if len(t) > 55 else t

# Wrap a label to 24 characters per line.
def wrap_label(t, width=24):
    return '\n'.join(textwrap.wrap(t, width=width))

def build_cnet(enr_sig, input_genes_upper, lfc_lookup_upper, title, outpath, 
               max_terms=15, filter_cfg=None):
    """Render one cnet figure: keyword-filtered terms, greedy set cover, spring layout."""
    sig = enr_sig.copy()
    
    # Keyword filters, when configured for this pathway
    if filter_cfg:
        term_kw = filter_cfg.get('term')
        gene_kw = filter_cfg.get('gene')
        if term_kw or gene_kw:
            mask = pd.Series(False, index=sig.index)
            if term_kw:
                mask |= sig['Term'].str.contains(term_kw, na=False)
            if gene_kw:
                mask |= sig['Genes'].str.upper().str.contains(gene_kw, na=False)
            sig = sig[mask]
    
    if len(sig) < 2:
        print(f"  SKIPPED {title}: too few terms after filter ({len(sig)})")
        return None
    
    # Genes of each term that are also in the input list
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
        print(f"  SKIPPED {title}: too few terms with gene matches ({len(term_data)})")
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
    
    print(f"  {title.split('–')[0].strip()}: {len(selected)} terms, {len(covered)} genes")
    
    # Build the term-gene graph
    G = nx.Graph()
    for td in selected:
        G.add_node(td['term'], ntype='term', lib=td['lib'])
        for g in td['input_genes']:
            G.add_node(g, ntype='gene', lfc=lfc_lookup_upper.get(g, 0))
            G.add_edge(td['term'], g)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
    
    term_nodes = [n for n, d in G.nodes(data=True) if d.get('ntype') == 'term']
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get('ntype') == 'gene']
    if len(gene_nodes) < 2:
        print(f"  SKIPPED: <2 genes")
        return None
    
    pos = nx.spring_layout(G, k=3.5, iterations=150, seed=42)
    fig, ax = plt.subplots(figsize=(30, 26))
    fig.patch.set_alpha(0)
    ax.set_facecolor('none')
    
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2, width=0.7, edge_color='#888888')
    
    # Gene nodes, coloured by mean log2FC
    gene_lfcs = [G.nodes[g].get('lfc', 0) for g in gene_nodes]
    gene_sizes = [max(80, G.degree(g) * 40) for g in gene_nodes]
    vmax = max(abs(min(gene_lfcs)), abs(max(gene_lfcs)), 0.5)
    gene_xy = np.array([pos[g] for g in gene_nodes])
    sc = ax.scatter(gene_xy[:, 0], gene_xy[:, 1], c=gene_lfcs,
                    cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                    s=gene_sizes, zorder=3, edgecolors='white', linewidths=0.5)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.35, pad=0.02)
    cbar.set_label('mean log₂FC', fontsize=20)
    
    # Term nodes, coloured by source library
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get('lib'), '#888888')
        xy = pos[t]
        ax.scatter(xy[0], xy[1], s=350, c=color, zorder=4,
                   edgecolors='white', linewidths=1.5, marker='s')
    
    # Term labels sit above their node
    for t in term_nodes:
        color = LIB_COLORS.get(G.nodes[t].get('lib'), '#888888')
        xy = pos[t]
        ax.text(xy[0], xy[1] + 0.035, wrap_label(t, 24),
                fontsize=20, fontweight='bold', color='white',
                ha='center', va='bottom', zorder=5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=color,
                          alpha=0.9, edgecolor='none'))
    
    # Gene labels sit below their node
    for g in gene_nodes:
        xy = pos[g]
        ax.text(xy[0], xy[1] - 0.03, g,
                fontsize=20, fontweight='bold', color='white',
                ha='center', va='top', zorder=5,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='#333333',
                          alpha=0.85, edgecolor='none'))
    
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
    print(f"    SAVED: {os.path.basename(outpath)} ({non_zero}/{len(gene_nodes)} genes with LFC≠0)")
    return True


# Pathway-based cnets, one per pathway, from the pathway-specific enrichment table.
print("=" * 60)
print("PATHWAY-BASED CNETS (from pathway-specific enrichment)")
print("=" * 60)

enr = pd.read_csv(f'{DATA}/pathway_specific_enrichment.csv')

pathway_genes_lists = {
    'JAK_STAT3': ['Jak1','Jak2','Jak3','Tyk2','Stat1','Stat2','Stat3','Stat4','Stat5a','Stat5b','Stat6',
                 'Socs1','Socs2','Socs3','Socs4','Socs5','Socs6','Socs7','Cish','Pias1','Pias2','Pias3','Pias4',
                 'Il6','Il6ra','Il6st','Lif','Lifr','Cntf','Cntfr','Osm','Osmr','Il10','Il10ra','Il10rb',
                 'Il21','Il21r','Il23a','Il12rb1','Ifng','Ifngr1','Ifngr2','Epo','Epor','Ghr','Prl','Prlr',
                 'Csf2','Csf2ra','Csf2rb','Il2','Il2ra','Il2rb','Il2rg','Il7','Il7r','Il15','Il15ra',
                 'Myc','Ccnd1','Bcl2','Bcl2l1','Mcl1','Pim1','Pim2','Vegfa','Hif1a','Mmp2','Mmp9','Birc5'],
    'BBB_MMP': ['Cldn1','Cldn3','Cldn5','Cldn11','Cldn12','Ocln','Tjp1','Tjp2','Tjp3','F11r','Jam2','Jam3',
            'Cdh5','Pecam1','Esam','Nectin1','Nectin2','Nectin3',
            'Abcb1a','Abcb1b','Abcg2','Abcc1','Abcc4','Slc2a1','Slc7a5','Slc16a1','Slc1a1',
            'Slc7a1','Slc3a2','Slc6a6','Mfsd2a',
            'Pdgfrb','Pdgfb','Tgfb1','Tgfbr1','Tgfbr2','Angpt1','Angpt2','Tek','Tie1',
            'Vegfa','Kdr','Flt1','Nrp1','Nrp2',
            'Aqp4','Gja1','Sox9','S100b','Gfap','Aldh1l1',
            'Mmp2','Mmp9','Mmp3','Timp1','Timp2','Timp3',
            'Icam1','Vcam1','Sele','Selp','Itgal','Itgb2','Cxcl12','Cxcr4'],
    'Akt_PI3K': ['Pik3ca','Pik3cb','Pik3cd','Pik3cg','Pik3r1','Pik3r2','Pik3r3',
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

titles = {
    'JAK_STAT3': 'JAK-STAT3 / Inflammatory Pathway – Concept Network',
    'BBB_MMP': 'MMP-2/9 & Blood-Brain Barrier – Concept Network',
    'Akt_PI3K': 'PI3K/Akt Pro-Survival Pathway – Concept Network',
}

for pw_name in ['JAK_STAT3', 'BBB_MMP', 'Akt_PI3K']:
    gene_list = pathway_genes_lists[pw_name]
    input_upper = {g.upper() for g in gene_list if g.upper() in symbol_upper}
    
    pw_enr = enr[enr['Pathway'] == pw_name].copy()
    sig = pw_enr[pw_enr['Adjusted P-value'] < 0.05]
    
    print(f"\n--- {pw_name} ---")
    print(f"  Input genes in data: {len(input_upper)}")
    print(f"  Significant terms: {len(sig)}")
    
    build_cnet(sig, input_upper, lfc_upper,
               titles[pw_name],
               f'{OUT_PW}/Cnet_{pw_name}_v2.png',
               max_terms=15, filter_cfg=FILTERS.get(pw_name))

print("\n\nPathway-based cnets done!")
