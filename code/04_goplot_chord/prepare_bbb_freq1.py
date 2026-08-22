#!/usr/bin/env python3
"""Build GOplot input tables for the BBB chord at a relaxed gene frequency.

BBB genes significant in >=1 group are used instead of the consensus >=2 rule.
Inputs: master_deg_LFC02.csv -> BBB_freq1_{terms,kegg_terms,reactome_terms,genes}.csv,
which the R driver turns into the GOChord_BBB_*_Freq1 plots.
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, re, sys
import pandas as pd
import gseapy as gp

BASE = BASE + "/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
OUT = BASE + "/mnt/Bulk RNA sequencing/GOPLOT analyis/data"
sys.path.insert(0, os.path.join(BASE, "scripts"))
from cnet_gene_lists import BBB_GENES

ALL_LIBS = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
            'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']
CAT = {'GO_Biological_Process_2023': 'BP', 'GO_Molecular_Function_2023': 'MF',
       'GO_Cellular_Component_2023': 'CC', 'KEGG_2021_Human': 'KEGG', 'Reactome_2022': 'Reactome'}


def parse_id(term):
    m = re.search(r'\(GO:(\d+)\)', term)
    if m:
        return f"GO:{m.group(1)}", re.sub(r'\s*\(GO:\d+\)', '', term).strip()
    m = re.search(r'\bR-HSA-(\d+)\b', term)
    if m:
        return f"R-HSA-{m.group(1)}", re.sub(r'\s*R-HSA-\d+', '', term).strip()
    return "", term.strip()


def terms_from(df):
    rows = []
    for _, r in df.iterrows():
        if r['Gene_set'] not in CAT:
            continue
        tid, name = parse_id(str(r['Term']))
        rows.append({'Category': CAT[r['Gene_set']], 'ID': tid, 'Term': name,
                     'Genes': str(r['Genes']).replace(';', ', '), 'adj_pval': float(r['Adjusted P-value'])})
    return pd.DataFrame(rows)


def main():
    m = pd.read_csv(os.path.join(DATA, "master_deg_LFC02.csv"))   # significant rows (>=1 group)
    bbb = {g.upper() for g in BBB_GENES}
    sub = m[m['symbol_key'].isin(bbb)]
    # one row per gene: number of groups it is significant in, plus mean logFC
    agg = sub.groupby('symbol').agg(n_groups=('group', 'nunique'), logFC=('log2FC', 'mean')).reset_index()
    sel = agg['symbol'].tolist()
    print(f"BBB genes present >=1 group: {len(sel)} (vs consensus >=2)")
    res = gp.enrichr(gene_list=sel, gene_sets=ALL_LIBS, organism='mouse', no_plot=True).results
    terms = terms_from(res)
    # per category keep terms at adj p < 0.05, or the top 12 by p-value if fewer than 3
    keep = []
    for cat in ['BP', 'MF', 'CC', 'KEGG', 'Reactome']:
        c = terms[terms['Category'] == cat]
        s = c[c['adj_pval'] < 0.05]
        keep.append(s if len(s) >= 3 else c.nsmallest(12, 'adj_pval'))
    terms = pd.concat(keep, ignore_index=True)
    go = terms[terms['Category'].isin(['BP', 'MF', 'CC'])]
    go.to_csv(os.path.join(OUT, "BBB_freq1_terms.csv"), index=False)
    terms[terms['Category'] == 'KEGG'].to_csv(os.path.join(OUT, "BBB_freq1_kegg_terms.csv"), index=False)
    terms[terms['Category'] == 'Reactome'].to_csv(os.path.join(OUT, "BBB_freq1_reactome_terms.csv"), index=False)
    agg[['symbol', 'logFC']].rename(columns={'symbol': 'ID'}).to_csv(os.path.join(OUT, "BBB_freq1_genes.csv"), index=False)
    print(f"wrote BBB_freq1 tables: {len(go)} GO, "
          f"{(terms['Category']=='KEGG').sum()} KEGG, {(terms['Category']=='Reactome').sum()} Reactome")


if __name__ == "__main__":
    main()
