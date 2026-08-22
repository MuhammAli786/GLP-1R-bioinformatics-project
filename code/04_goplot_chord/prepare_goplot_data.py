#!/usr/bin/env python3
"""Build the input tables GOplot's circle_dat() needs from the consensus meta-analysis.

Inputs: consensus_LFC02.csv and enrichment_LFC02.csv -> <name>_terms.csv
(Category, ID, Term, Genes, adj_pval) and <name>_genes.csv (ID, logFC).
Two sets are produced: consensus, from the existing consensus GO enrichment; and
restricted, where for each curated pathway (BBB, Inflammatory, Survival,
IonChannel) the consensus genes are intersected with that pathway's gene list and
re-enriched against GO BP/MF/CC, KEGG and Reactome.
Requires a python with pandas and gseapy.
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
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(BASE, "scripts"))
from cnet_gene_lists import BBB_GENES, JAK_STAT3_GENES, PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES

GO_LIBS = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023']
ALL_LIBS = GO_LIBS + ['KEGG_2021_Human', 'Reactome_2022']
CAT = {'GO_Biological_Process_2023': 'BP', 'GO_Molecular_Function_2023': 'MF',
       'GO_Cellular_Component_2023': 'CC',
       'KEGG_2021_Human': 'KEGG', 'Reactome_2022': 'Reactome'}
THR = "LFC02"


def parse_id(term):
    """Strip and capture a GO:nnn or R-HSA-nnn id from a term name."""
    m = re.search(r'\(GO:(\d+)\)', term)
    if m:
        return f"GO:{m.group(1)}", re.sub(r'\s*\(GO:\d+\)', '', term).strip()
    m = re.search(r'\bR-HSA-(\d+)\b', term)
    if m:
        return f"R-HSA-{m.group(1)}", re.sub(r'\s*R-HSA-\d+', '', term).strip()
    return "", term.strip()


def terms_from_enrichment(df):
    """Turn gseapy results into a tidy term table with Category in BP/MF/CC/KEGG/Reactome."""
    rows = []
    for _, r in df.iterrows():
        if r['Gene_set'] not in CAT:
            continue
        tid, name = parse_id(str(r['Term']))
        rows.append({
            'Category': CAT[r['Gene_set']],
            'ID': tid,
            'Term': name,
            'Genes': str(r['Genes']).replace(';', ', '),
            'adj_pval': float(r['Adjusted P-value']),
        })
    return pd.DataFrame(rows)


def write_pair(name, terms_df, genes_df):
    """Write the GO term and gene tables, plus separate KEGG and Reactome term tables."""
    go = terms_df[terms_df['Category'].isin(['BP', 'MF', 'CC'])]
    go.to_csv(os.path.join(OUT, f"{name}_terms.csv"), index=False)
    genes_df.to_csv(os.path.join(OUT, f"{name}_genes.csv"), index=False)
    msg = f"  {name}: {len(go)} GO"
    for cat, suffix in [('KEGG', 'kegg'), ('Reactome', 'reactome')]:
        sub = terms_df[terms_df['Category'] == cat]
        sub.to_csv(os.path.join(OUT, f"{name}_{suffix}_terms.csv"), index=False)
        msg += f", {len(sub)} {cat}"
    print(msg + f", {len(genes_df)} genes")


def main():
    cons = pd.read_csv(os.path.join(DATA, f"consensus_{THR}.csv"))     # gene_symbol, n_groups, mean_log2FC
    cons_lfc = dict(zip(cons['gene_symbol'], cons['mean_log2FC']))

    # Consensus set: reuse the precomputed enrichment, keeping adj p < 0.05
    enr = pd.read_csv(os.path.join(DATA, f"enrichment_{THR}.csv"))
    enr = enr[enr['Adjusted P-value'] < 0.05]
    cterms = terms_from_enrichment(enr)
    cgenes = cons[['gene_symbol', 'mean_log2FC']].rename(columns={'gene_symbol': 'ID', 'mean_log2FC': 'logFC'})
    write_pair("consensus", cterms, cgenes)

    # Restricted sets: one re-enrichment per curated pathway
    MECH = {"BBB": BBB_GENES, "Inflammatory": JAK_STAT3_GENES,
            "Survival": PI3K_AKT_GENES, "IonChannel": ION_CHANNEL_BASE_GENES}
    cons_keys = {g.upper() for g in cons['gene_symbol']}
    for name, glist in MECH.items():
        sel = [g for g in glist if g.upper() in cons_keys]   # mechanism genes that are consensus
        if len(sel) < 3:
            print(f"  {name}: only {len(sel)} consensus genes — skipped")
            continue
        try:
            res = gp.enrichr(gene_list=sel, gene_sets=ALL_LIBS, organism='mouse', no_plot=True).results
        except Exception as e:
            print(f"  {name}: enrichr ERR {str(e)[:80]}"); continue
        terms = terms_from_enrichment(res)
        # per category keep terms at adj p < 0.05, or the top 12 by p-value if fewer than 3
        keep = []
        for cat in ['BP', 'MF', 'CC', 'KEGG', 'Reactome']:
            c = terms[terms['Category'] == cat]
            sig = c[c['adj_pval'] < 0.05]
            keep.append(sig if len(sig) >= 3 else c.nsmallest(12, 'adj_pval'))
        terms = pd.concat(keep, ignore_index=True)
        genes = pd.DataFrame({'ID': sel, 'logFC': [cons_lfc.get(g, 0.0) for g in sel]})
        write_pair(name, terms, genes)


if __name__ == "__main__":
    main()
