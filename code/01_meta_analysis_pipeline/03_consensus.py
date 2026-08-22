#!/usr/bin/env python3
"""Identify consensus genes: a gene is consensus if it is a significant DEG in >= 2 distinct groups, computed separately per threshold.

Data/master_deg_<thr>.csv -> consensus_<thr>.csv (gene_symbol, n_groups, predominant_direction, mean_log2FC, groups), gene_lfc_<thr>.csv (long per-group lookup for consensus genes) and gene_lfc_comprehensive_<thr>.csv (mean log2FC for all significant genes).
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os
import pandas as pd
import numpy as np

OUTDATA = BASE + "/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/Data"
MIN_GROUPS = 2

for thr in ["LFC02", "LFC05", "LFC1"]:
    m = pd.read_csv(os.path.join(OUTDATA, f"master_deg_{thr}.csv"))
    if m.empty:
        print(thr, "empty master"); continue
    g = m.groupby("symbol_key")
    rows = []
    for key, sub in g:
        groups = sub["group"].unique()
        n = len(groups)
        if n < MIN_GROUPS:
            continue
        disp = sub["symbol"].iloc[0]
        up = (sub["direction"] == "UP").sum()
        down = (sub["direction"] == "DOWN").sum()
        rows.append({
            "gene_symbol": disp,
            "n_groups": n,
            "predominant_direction": "UP" if up >= down else "DOWN",
            "n_up": int(up), "n_down": int(down),
            "mean_log2FC": sub["log2FC"].mean(),
            "groups": ";".join(sorted(groups)),
        })
    cons = pd.DataFrame(rows).sort_values(["n_groups", "mean_log2FC"],
                                          ascending=[False, False])
    cons.to_csv(os.path.join(OUTDATA, f"consensus_{thr}.csv"), index=False)

    cons_keys = set(cons["gene_symbol"].str.upper())
    lfc = m[m["symbol_key"].isin(cons_keys)][["symbol", "group", "log2FC", "padj", "direction"]]
    lfc = lfc.rename(columns={"symbol": "gene_symbol"})
    lfc.to_csv(os.path.join(OUTDATA, f"gene_lfc_{thr}.csv"), index=False)

    # Mean log2FC over all significant genes, used to colour the cnet plots.
    comp = m.groupby("symbol").agg(mean_log2FC=("log2FC", "mean")).reset_index()
    comp = comp.rename(columns={"symbol": "gene_symbol"})
    comp.to_csv(os.path.join(OUTDATA, f"gene_lfc_comprehensive_{thr}.csv"), index=False)

    print(f"{thr}: {len(cons)} consensus genes (>= {MIN_GROUPS} groups); "
          f"max n_groups={cons['n_groups'].max() if len(cons) else 0}")
    print("  top 10:", ", ".join(f"{r.gene_symbol}({r.n_groups})"
                                  for r in cons.head(10).itertuples()))
