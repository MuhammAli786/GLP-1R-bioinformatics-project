import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

OUT = "/sessions/practical-ecstatic-mendel/mnt/outputs"
PLOT_DIR = os.path.join(OUT, "final_plots")
os.makedirs(PLOT_DIR, exist_ok=True)


# ENRICHMENT DOT PLOTS
def make_dotplot(enr_df, db_name, threshold_label, top_n=20):
    """Create enrichment dot plot for a single database"""
    subset = enr_df[(enr_df['Database'] == db_name) & (enr_df['Adjusted P-value'] < 0.05)]
    if subset.empty:
        return
    
    # Get Combined_Score column
    score_col = 'Combined_Score' if 'Combined_Score' in subset.columns else 'Combined Score'
    if score_col not in subset.columns:
        return
    
    subset = subset.nlargest(top_n, score_col)
    subset = subset.sort_values(score_col, ascending=True)
    
    fig, ax = plt.subplots(figsize=(12, max(6, len(subset)*0.4)))
    
    # Parse overlap for dot size
    def get_overlap_count(x):
        try:
            return int(str(x).split('/')[0])
        except:
            return 5
    
    subset['n_genes'] = subset['Overlap'].apply(get_overlap_count)
    
    scatter = ax.scatter(subset[score_col], range(len(subset)),
                        s=subset['n_genes']*8, 
                        c=-np.log10(subset['Adjusted P-value']),
                        cmap='RdYlBu_r', alpha=0.8, edgecolors='black', linewidths=0.5)
    
    ax.set_yticks(range(len(subset)))
    ax.set_yticklabels(subset['Term'].values, fontsize=9)
    ax.set_xlabel('Combined Score', fontsize=12)
    ax.set_title(f'{db_name.replace("_", " ")} - {threshold_label}\nTop {len(subset)} Enriched Terms', fontsize=13)
    
    plt.colorbar(scatter, label='-log10(adj p-value)', shrink=0.7)
    plt.tight_layout()
    
    safe_db = db_name.replace(' ', '_')
    fname = f"dotplot_{safe_db}_{threshold_label}.png"
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, transparent=True, bbox_inches='tight')
    plt.close()
    return fname

databases = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']

for thr_label, enr_file in [('LFC05', 'final_enrichment_LFC05.csv'), ('LFC0', 'final_enrichment_LFC0.csv')]:
    enr = pd.read_csv(os.path.join(OUT, enr_file))
    if 'Combined Score' in enr.columns and 'Combined_Score' not in enr.columns:
        enr['Combined_Score'] = enr['Combined Score']
    
    for db in databases:
        fname = make_dotplot(enr, db, thr_label)
        if fname:
            print(f"  Created: {fname}")


# TOP 20 GENE BARPLOT (by group frequency)
def make_top_gene_barplot(cons_df, threshold_label):
    top = cons_df.head(20).copy()
    
    fig, ax = plt.subplots(figsize=(12, 8))
    colors = ['#d73027' if d == 'UP' else '#4575b4' for d in top['predominant_direction']]
    
    bars = ax.barh(range(len(top)), top['n_groups'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['symbol'], fontsize=11)
    ax.set_xlabel('Number of Treatment×Region Groups', fontsize=12)
    ax.set_title(f'Top 20 Consensus Genes by Group Frequency ({threshold_label})', fontsize=14)
    ax.invert_yaxis()
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#d73027', label='Predominantly UP'),
                       Patch(facecolor='#4575b4', label='Predominantly DOWN')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=11)
    
    plt.tight_layout()
    fname = f"top20_genes_barplot_{threshold_label}.png"
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, transparent=True, bbox_inches='tight')
    plt.close()
    print(f"  Created: {fname}")

for thr_label, cons_file in [('LFC05', 'final_consensus_LFC05.csv'), ('LFC0', 'final_consensus_LFC0.csv')]:
    cons = pd.read_csv(os.path.join(OUT, cons_file))
    make_top_gene_barplot(cons, thr_label)


# HEATMAP: Top genes × groups (average LFC)
def make_heatmap(master_file, cons_file, threshold_label, top_n=30):
    master = pd.read_csv(os.path.join(OUT, master_file))
    cons = pd.read_csv(os.path.join(OUT, cons_file))
    top_genes = cons.head(top_n)['symbol'].tolist()
    
    master['group'] = master['treatment_class'] + '|' + master['region'].fillna('unknown') if 'treatment_class' in master.columns else master['group']
    
    sub = master[master['symbol'].isin(top_genes)]
    pivot = sub.pivot_table(index='symbol', columns='group', values='log2FC', aggfunc='mean')
    
    # Reorder by consensus frequency
    pivot = pivot.reindex(top_genes)
    # Drop columns with all NaN
    pivot = pivot.dropna(axis=1, how='all')
    
    fig, ax = plt.subplots(figsize=(max(14, len(pivot.columns)*0.8), max(8, len(pivot)*0.4)))
    sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax, 
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Mean log2FC', 'shrink': 0.7},
                xticklabels=True, yticklabels=True)
    
    ax.set_title(f'Top {top_n} Consensus Genes: Mean LFC by Group ({threshold_label})', fontsize=13)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    
    fname = f"heatmap_top{top_n}_{threshold_label}.png"
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, transparent=True, bbox_inches='tight')
    plt.close()
    print(f"  Created: {fname}")

make_heatmap('final_master_deg.csv', 'final_consensus_LFC05.csv', 'LFC05')
make_heatmap('final_master_deg_LFC0.csv', 'final_consensus_LFC0.csv', 'LFC0')

print(f"\nAll plots saved to: {PLOT_DIR}")
print(f"Total files: {len(os.listdir(PLOT_DIR))}")
