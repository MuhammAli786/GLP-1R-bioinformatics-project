#!/usr/bin/env python3
"""
02_run_gsea.py  — GSEA prerank per group (Hallmark + KEGG), resumable.
Output: GSEA/data/gsea/<group>__<lib>.csv  (Term, NES, NOM p-val, FDR q-val, Lead_genes)
Run repeatedly until it prints ALL_DONE (time-budgeted per call).
"""
import os, glob, time, sys
import pandas as pd
import gseapy as gp

ROOT = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/GSEA"
RNK = os.path.join(ROOT, "data", "rnk")
OUT = os.path.join(ROOT, "data", "gsea")
os.makedirs(OUT, exist_ok=True)
LIBS = {"hallmark": "MSigDB_Hallmark_2020", "kegg": "KEGG_2021_Human"}

deadline = time.time() + 36
done = 0
jobs = [(os.path.splitext(os.path.basename(f))[0], f, k, v)
        for f in sorted(glob.glob(os.path.join(RNK, "*.rnk")))
        for k, v in LIBS.items()]
remaining = [(g, f, k, v) for g, f, k, v in jobs
             if not os.path.exists(os.path.join(OUT, f"{g}__{k}.csv"))]
print(f"{len(jobs)-len(remaining)}/{len(jobs)} already done; {len(remaining)} remaining")

for g, f, key, lib in remaining:
    if time.time() > deadline:
        break
    try:
        pre = gp.prerank(rnk=f, gene_sets=lib, min_size=10, max_size=500,
                         permutation_num=100, seed=42, threads=4, no_plot=True, outdir=None, verbose=False)
        r = pre.res2d[["Term", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]].copy()
        r.to_csv(os.path.join(OUT, f"{g}__{key}.csv"), index=False)
        done += 1
        print(f"  {g} [{key}]: {len(r)} terms")
    except Exception as e:
        print(f"  {g} [{key}] ERR {str(e)[:90]}")
        time.sleep(2)

left = sum(1 for g, f, k, v in jobs if not os.path.exists(os.path.join(OUT, f"{g}__{k}.csv")))
print(f"did {done} this run; {left} still remaining")
if left == 0:
    print("ALL_DONE")
