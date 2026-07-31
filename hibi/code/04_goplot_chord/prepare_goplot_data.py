#!/usr/bin/env python3
"""
hibi_prepare_goplot_data.py
-------------------------------------------------------------
Port of the GLP-1R project's 04_goplot_chord/prepare_goplot_data.py.

Builds the input tables GOplot's circle_dat() needs, from the HIBI
consensus meta-analysis outputs:

  <name>_terms.csv : Category, ID, Term, Genes, adj_pval   (GO BP/MF/CC)
  <name>_kegg_terms.csv / <name>_reactome_terms.csv
  <name>_genes.csv : ID, logFC

Two kinds of inputs:
  * Consensus  -- the HIBI consensus genes, enriched against GO+KEGG+Reactome.
  * Restricted -- for each curated pathway (BBB, Inflammatory, Survival,
    IonChannel) the consensus genes are intersected with the pathway gene
    list and RE-ENRICHED so the term/gene tables are mechanism-specific.

Note: unlike GLP-1R (which reused a precomputed enrichment_LFC02.csv), the
consensus enrichment is computed here. Enrichr caps submitted lists, so the
consensus set is submitted as its top CONS_MAX genes by recurrence
(n_groups, tie-broken by |mean_log2FC|) -- the same ranking the consensus
Cnet and heatmap use.
"""
import os, re, sys, time
import pandas as pd
import gseapy as gp

DATA = "/sessions/lucid-pensive-ride/mnt/outputs/hibi_data/meta_analysis"
OUT = "/sessions/lucid-pensive-ride/mnt/outputs/hibi_data/goplot"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnet_gene_lists import BBB_GENES, JAK_STAT3_GENES, PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES

GO_LIBS = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023']
ALL_LIBS = GO_LIBS + ['KEGG_2021_Human', 'Reactome_2022']
CAT = {'GO_Biological_Process_2023': 'BP', 'GO_Molecular_Function_2023': 'MF',
       'GO_Cellular_Component_2023': 'CC',
       'KEGG_2021_Human': 'KEGG', 'Reactome_2022': 'Reactome'}
THR = "LFC02"
CONS_MAX = 1500


def parse_id(term):
    m = re.search(r'\(GO:(\d+)\)', term)
    if m:
        return f"GO:{m.group(1)}", re.sub(r'\s*\(GO:\d+\)', '', term).strip()
    m = re.search(r'\bR-HSA-(\d+)\b', term)
    if m:
        return f"R-HSA-{m.group(1)}", re.sub(r'\s*R-HSA-\d+', '', term).strip()
    return "", term.strip()


def terms_from_enrichment(df):
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


def enrich(gene_list, libs=ALL_LIBS):
    """Enrichr with per-library retry (Enrichr 429s on rapid repeat calls)."""
    frames = []
    for db in libs:
        for attempt in range(5):
            try:
                r = gp.enrichr(gene_list=gene_list, gene_sets=db, organism='mouse', no_plot=True).results
                r['Gene_set'] = db
                frames.append(r)
                break
            except Exception as e:
                if attempt == 4:
                    print(f"    {db}: FAILED {str(e)[:70]}")
                time.sleep(3 * (attempt + 1))
        time.sleep(1)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def write_pair(name, terms_df, genes_df):
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
    cons = pd.read_csv(os.path.join(DATA, f"consensus_{THR}.csv"))
    cons_lfc = dict(zip(cons['gene_symbol'], cons['mean_log2FC']))

    # ---- Consensus ----
    ranked = cons.assign(_abs=cons['mean_log2FC'].abs()).sort_values(
        ['n_groups', '_abs'], ascending=[False, False])
    sub = ranked.head(CONS_MAX)
    print(f"Consensus: enriching top {len(sub)} of {len(cons)} consensus genes "
          f"(n_groups {sub['n_groups'].min()}-{sub['n_groups'].max()})")
    res = enrich(sub['gene_symbol'].tolist())
    cterms = terms_from_enrichment(res)
    cterms = cterms[cterms['adj_pval'] < 0.05]
    cgenes = sub[['gene_symbol', 'mean_log2FC']].rename(
        columns={'gene_symbol': 'ID', 'mean_log2FC': 'logFC'})
    write_pair("consensus", cterms, cgenes)

    # ---- Restricted pathways ----
    MECH = {"BBB": BBB_GENES, "Inflammatory": JAK_STAT3_GENES,
            "Survival": PI3K_AKT_GENES, "IonChannel": ION_CHANNEL_BASE_GENES}
    cons_keys = {g.upper() for g in cons['gene_symbol']}
    for name, glist in MECH.items():
        sel = [g for g in glist if g.upper() in cons_keys]
        if len(sel) < 3:
            print(f"  {name}: only {len(sel)} consensus genes - skipped")
            continue
        res = enrich(sel)
        if res.empty:
            print(f"  {name}: enrichment failed"); continue
        terms = terms_from_enrichment(res)
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
