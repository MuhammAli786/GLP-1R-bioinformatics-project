# Enrichr over-representation analysis of the cross-accession consensus gene lists.
# final_consensus_LFC05.csv / final_consensus_LFC0.csv -> final_enrichment_*.csv

import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import pandas as pd
import gseapy as gp
import os

OUT = BASE + "/mnt/outputs"
databases = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']

# Query Enrichr for a mouse gene list across all five libraries and concatenate the results.
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

# |log2FC| >= 0.5 consensus set
cons05 = pd.read_csv(os.path.join(OUT, 'final_consensus_LFC05.csv'))
genes05 = cons05['symbol'].tolist()
enr05 = run_enrichment(genes05, "LFC≥0.5")
enr05.to_csv(os.path.join(OUT, 'final_enrichment_LFC05.csv'), index=False)
print(f"\nSaved final_enrichment_LFC05.csv ({len(enr05)} rows)")

# |log2FC| >= 0 consensus set, capped at the top 3000 genes by group frequency for Enrichr
cons0 = pd.read_csv(os.path.join(OUT, 'final_consensus_LFC0.csv'))
genes0 = cons0.head(3000)['symbol'].tolist()
enr0 = run_enrichment(genes0, "LFC≥0 (top 3000)")
enr0.to_csv(os.path.join(OUT, 'final_enrichment_LFC0.csv'), index=False)
print(f"Saved final_enrichment_LFC0.csv ({len(enr0)} rows)")
