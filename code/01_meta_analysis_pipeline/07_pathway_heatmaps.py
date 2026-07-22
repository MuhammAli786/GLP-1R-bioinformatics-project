#!/usr/bin/env python3
"""
07_pathway_heatmaps.py
-------------------------------------------------------------
Task 4: mechanism-specific heatmaps built from the curated
cnet_gene_lists.py sets (BBB, JAK-STAT3 Inflammatory, PI3K/Akt
Survival).  Shows the TOP-10 most relevant genes in each mechanism
(ranked by group recurrence, then mean |log2FC|).

Special rule (BBB): force-include ALL MMP genes (and MMP-9) that are
significant in the data, in addition to the top-10.

Same styling as the consensus heatmap: Y = genes, X = groups labelled
A, B, C ... with a separate group_legend.txt; transparent PNG + PDF.
Saved to Plots/Heatmaps/<mechanism>/<LFCfolder>/<PDF|PNG>/
"""
import os, string, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnet_gene_lists import BBB_GENES, JAK_STAT3_GENES, PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES

BASE = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
PLOTS = os.path.join(BASE, "Plots", "Heatmaps")
THR_FOLDER = {"LFC02": "LFC0.2", "LFC05": "LFC0.5", "LFC1": "LFC1"}

MECHANISMS = {
    "BBB": BBB_GENES,
    "Inflammatory": JAK_STAT3_GENES,
    "Survival": PI3K_AKT_GENES,
    "IonChannel": ION_CHANNEL_BASE_GENES,
}


def excel_letters(n):
    out = []
    for i in range(n):
        s = ""; x = i
        while True:
            s = string.ascii_uppercase[x % 26] + s
            x = x // 26 - 1
            if x < 0:
                break
        out.append(s)
    return out


def make(mech, genes, thr, master, cons):
    gene_set = {g.upper() for g in genes}
    # candidate mechanism genes that are significant in the data
    sub_all = master[master["symbol_key"].isin(gene_set)]
    if sub_all.empty:
        print(f"  {mech}/{thr}: no significant mechanism genes"); return
    # rank by recurrence then mean|lfc|
    rank = (sub_all.groupby("symbol")
            .agg(n_groups=("group", "nunique"), m=("log2FC", lambda x: np.mean(np.abs(x))))
            .reset_index().sort_values(["n_groups", "m"], ascending=[False, False]))
    top = rank.head(20)["symbol"].tolist()
    sel = list(top)
    forced = []
    if mech == "BBB":
        # force-include ALL curated MMP proteins (incl. Mmp9) as explicit rows,
        # even when not significant in any group (shown as blank cells)
        forced = [g for g in genes if g.upper().startswith("MMP")]
        # plus any other significant MMP detected
        forced += [s for s in sub_all["symbol"].unique() if s.upper().startswith("MMP")]
        for m in forced:
            if m not in sel:
                sel.append(m)
    sub = sub_all[sub_all["symbol"].isin(sel)]
    pivot = sub.pivot_table(index="symbol", columns="group", values="log2FC", aggfunc="mean")
    # drop groups (columns) that have no data for the displayed genes
    pivot = pivot.dropna(axis=1, how="all")
    # order rows: top-ranked first, then forced MMP rows (blank rows kept)
    order, seen = [], set()
    for g in top + forced:
        if g not in seen:
            order.append(g); seen.add(g)
    cols = sorted(pivot.columns)
    pivot = pivot.reindex(index=order, columns=cols)
    if pivot.empty or len(cols) == 0:
        print(f"  {mech}/{thr}: empty after pivot"); return
    letters = excel_letters(len(cols))
    outdir = os.path.join(PLOTS, mech, THR_FOLDER[thr])
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "group_legend.txt"), "w") as f:
        f.write(f"Group legend — {mech} mechanism heatmap ({THR_FOLDER[thr]})\nLetter\tGroup\n")
        for l, c in zip(letters, cols):
            f.write(f"{l}\t{c}\n")
    fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.6), max(5, len(order) * 0.5)))
    _v = np.abs(pivot.values); _v = _v[np.isfinite(_v)]
    vmax = max(1.5, float(np.nanpercentile(_v, 90))) if _v.size else 1.5
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(letters, fontsize=12)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(pivot.index, fontsize=11)
    ax.set_xlabel("Group (see group_legend.txt)", fontsize=12)
    ax.set_title(f"{mech}: top relevant genes — mean log2FC by group ({THR_FOLDER[thr]})", fontsize=12)
    plt.colorbar(im, label="mean log2FC", shrink=0.7)
    for fmt in ("PNG", "PDF"):
        d = os.path.join(outdir, fmt); os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, f"heatmap_{mech}.{fmt.lower()}"),
                    dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  {mech}/{thr}: {len(order)} genes x {len(cols)} groups")


def main():
    for thr in ["LFC02", "LFC05", "LFC1"]:
        master = pd.read_csv(os.path.join(DATA, f"master_deg_{thr}.csv"))
        cons = pd.read_csv(os.path.join(DATA, f"consensus_{thr}.csv"))
        print(f"=== {thr} ===")
        for mech, genes in MECHANISMS.items():
            make(mech, genes, thr, master, cons)


if __name__ == "__main__":
    main()
