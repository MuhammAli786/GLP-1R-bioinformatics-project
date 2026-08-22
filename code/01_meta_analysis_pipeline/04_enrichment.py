#!/usr/bin/env python3
"""Functional enrichment of the consensus gene lists against GO, KEGG and Reactome via Enrichr (gseapy, organism=mouse).

Data/consensus_<thr>.csv -> Data/enrichment_<thr>.csv, per threshold.
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, time, sys
import pandas as pd
import gseapy as gp

OUTDATA = BASE + "/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/Data"
DATABASES = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']


def run(genes, label):
    """Run Enrichr over every database, retrying up to four times each."""
    allr = []
    for db in DATABASES:
        for attempt in range(4):
            try:
                enr = gp.enrichr(gene_list=genes, gene_sets=db, organism='mouse', no_plot=True)
                res = enr.results.copy()
                res['Database'] = db
                if 'Gene_set' not in res.columns:
                    res['Gene_set'] = db
                allr.append(res)
                nsig = (res['Adjusted P-value'] < 0.05).sum()
                print(f"  [{label}] {db}: {len(res)} terms, {nsig} sig")
                break
            except Exception as e:
                print(f"  [{label}] {db} attempt {attempt+1} ERR {str(e)[:80]}")
                time.sleep(3*(attempt+1))
    return pd.concat(allr, ignore_index=True) if allr else pd.DataFrame()


THRS = [sys.argv[1]] if len(sys.argv) > 1 else ["LFC02", "LFC05", "LFC1"]
for thr in THRS:
    cons = pd.read_csv(os.path.join(OUTDATA, f"consensus_{thr}.csv"))
    genes = cons['gene_symbol'].dropna().astype(str).tolist()
    print(f"=== {thr}: {len(genes)} consensus genes ===")
    enr = run(genes, thr)
    enr.to_csv(os.path.join(OUTDATA, f"enrichment_{thr}.csv"), index=False)
    print(f"Saved enrichment_{thr}.csv ({len(enr)} rows)\n")
