#!/usr/bin/env python3
"""GSEA prerank per group against MSigDB Hallmark and KEGG 2021 Human; resumable.

Inputs: GSEA/data/rnk/*.rnk -> GSEA/data/gsea/<group>__<lib>.csv (Term, NES,
NOM p-val, FDR q-val, Lead_genes). Each call is time-budgeted; re-run until it
prints ALL_DONE.
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, glob, time, sys
import pandas as pd
import gseapy as gp

ROOT = BASE + "/mnt/Bulk RNA sequencing/GSEA"
RNK = os.path.join(ROOT, "data", "rnk")
OUT = os.path.join(ROOT, "data", "gsea")
os.makedirs(OUT, exist_ok=True)
LIBS = {"hallmark": "MSigDB_Hallmark_2020", "kegg": "KEGG_2021_Human"}

deadline = time.time() + 36  # per-call wall-clock budget in seconds
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
