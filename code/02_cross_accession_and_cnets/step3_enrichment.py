import pandas as pd
import gseapy as gp
import os

OUT = "/sessions/practical-ecstatic-mendel/mnt/outputs"
databases = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']

def run_enrichment(gene_list, label):
    print(f"\n=== Running enrichment: {label} ({len(gene_list)} genes) ===")
    all_results = []
    for db in databases:
        try:
            enr = gp.enrichr(gene_list=gene_list, gene_sets=db, organism='mouse', no_plot=True)
            res = enr.results.copy()
            res['Database'] = db
            all_results.append(res)
            n_sig = (res['Adjusted P-value'] < 0.05).sum()
            print(f"  {db}: {len(res)} terms, {n_sig} significant")
        except Exception as e:
            print(f"  {db}: ERROR - {e}")
    
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined['Combined_Score'] = combined['Combined Score'] if 'Combined Score' in combined.columns else 0
        return combined
    return pd.DataFrame()

# LFC >= 0.5
cons05 = pd.read_csv(os.path.join(OUT, 'final_consensus_LFC05.csv'))
genes05 = cons05['symbol'].tolist()
enr05 = run_enrichment(genes05, "LFC≥0.5")
enr05.to_csv(os.path.join(OUT, 'final_enrichment_LFC05.csv'), index=False)
print(f"\nSaved final_enrichment_LFC05.csv ({len(enr05)} rows)")

# LFC >= 0 (use top 3000 by n_groups to keep it manageable for enrichr)
cons0 = pd.read_csv(os.path.join(OUT, 'final_consensus_LFC0.csv'))
genes0 = cons0.head(3000)['symbol'].tolist()  # Top 3000 by group frequency
enr0 = run_enrichment(genes0, "LFC≥0 (top 3000)")
enr0.to_csv(os.path.join(OUT, 'final_enrichment_LFC0.csv'), index=False)
print(f"Saved final_enrichment_LFC0.csv ({len(enr0)} rows)")
