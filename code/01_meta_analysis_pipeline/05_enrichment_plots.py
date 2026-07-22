#!/usr/bin/env python3
"""
05_enrichment_plots.py
-------------------------------------------------------------
Task 2 deliverables:
  * Enrichment DOT plots  (per database, top terms)
  * Enrichment BAR plots  (per database) + top consensus-gene bar plot
  * HEATMAP of the top-10 most enriched consensus genes
        (ranked by group frequency); Y = genes, X = groups labelled
        A, B, C, D ...  with a separate legend .txt mapping letter->group.

All plots saved as transparent PNG and PDF into the nested folders:
  Plots/<type>/<LFCfolder>/<PDF|PNG>/
"""
import os, string
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow"
DATA = os.path.join(BASE, "Data")
PLOTS = os.path.join(BASE, "Plots")
THR_FOLDER = {"LFC02": "LFC0.2", "LFC05": "LFC0.5", "LFC1": "LFC1"}
DATABASES = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
             'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']
DB_SHORT = {'GO_Biological_Process_2023': 'GO_BP', 'GO_Molecular_Function_2023': 'GO_MF',
            'GO_Cellular_Component_2023': 'GO_CC', 'KEGG_2021_Human': 'KEGG',
            'Reactome_2022': 'Reactome'}


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


def save(fig, ptype, thr, fname):
    for fmt in ("PNG", "PDF"):
        d = os.path.join(PLOTS, ptype, THR_FOLDER[thr], fmt)
        os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, f"{fname}.{fmt.lower()}"),
                    dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)


def clean_term(t):
    import re
    t = re.sub(r"\s*\(GO:\d+\)", "", t)
    t = re.sub(r"\s*R-HSA-\d+", "", t)
    return t[:60] + "..." if len(t) > 63 else t


def dotplot(enr, db, thr, top_n=20):
    sub = enr[(enr["Gene_set"] == db) & (enr["Adjusted P-value"] < 0.05)].copy()
    if sub.empty:
        return
    sub["n_genes"] = sub["Overlap"].apply(lambda x: int(str(x).split("/")[0]))
    sub = sub.nlargest(top_n, "Combined Score").sort_values("Combined Score")
    fig, ax = plt.subplots(figsize=(12, max(5, len(sub) * 0.42)))
    sc = ax.scatter(sub["Combined Score"], range(len(sub)),
                    s=sub["n_genes"] * 22,
                    c=-np.log10(sub["Adjusted P-value"]),
                    cmap="plasma", alpha=0.9, edgecolors="black", linewidths=0.6)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([clean_term(t) for t in sub["Term"]], fontsize=9)
    ax.set_xlabel("Combined Score", fontsize=12)
    ax.set_title(f"{DB_SHORT[db]} enrichment — {THR_FOLDER[thr]}\n"
                 f"Top {len(sub)} terms (padj<0.05)", fontsize=13)
    cb = plt.colorbar(sc, label="-log10(adj p)", shrink=0.7)
    save(fig, "Dot plots", thr, f"dotplot_{DB_SHORT[db]}")


def barplot_db(enr, db, thr, top_n=20):
    sub = enr[(enr["Gene_set"] == db) & (enr["Adjusted P-value"] < 0.05)].copy()
    if sub.empty:
        return
    sub["score"] = -np.log10(sub["Adjusted P-value"])
    sub = sub.nlargest(top_n, "score").sort_values("score")
    fig, ax = plt.subplots(figsize=(11, max(4, len(sub) * 0.42)))
    ax.barh(range(len(sub)), sub["score"], color="#4575b4", edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([clean_term(t) for t in sub["Term"]], fontsize=9)
    ax.set_xlabel("-log10(adj p-value)", fontsize=12)
    ax.set_title(f"{DB_SHORT[db]} — Top {len(sub)} enriched terms ({THR_FOLDER[thr]})", fontsize=12)
    save(fig, "Bar plots", thr, f"barplot_{DB_SHORT[db]}")


def top_gene_barplot(cons, thr, top_n=20):
    top = cons.head(top_n)
    fig, ax = plt.subplots(figsize=(11, 8))
    colors = ["#d73027" if d == "UP" else "#4575b4" for d in top["predominant_direction"]]
    ax.barh(range(len(top)), top["n_groups"], color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene_symbol"], fontsize=11)
    ax.set_xlabel("Number of groups (recurrence)", fontsize=12)
    ax.set_title(f"Top {top_n} consensus genes by group frequency ({THR_FOLDER[thr]})", fontsize=13)
    ax.invert_yaxis()
    ax.legend(handles=[Patch(facecolor="#d73027", label="Predominantly UP"),
                       Patch(facecolor="#4575b4", label="Predominantly DOWN")],
              loc="lower right", fontsize=10)
    save(fig, "Bar plots", thr, "top20_consensus_genes_barplot")


def consensus_heatmap(cons, master, thr, top_n=20):
    genes = cons.head(top_n)["gene_symbol"].tolist()
    sub = master[master["symbol"].isin(genes)]
    pivot = sub.pivot_table(index="symbol", columns="group", values="log2FC", aggfunc="mean")
    # drop groups (columns) that have no data for the displayed genes
    pivot = pivot.reindex(index=genes).dropna(axis=1, how="all")
    cols = sorted(pivot.columns)
    pivot = pivot[cols]
    letters = excel_letters(len(cols))
    # legend txt
    leg_dir = os.path.join(PLOTS, "Heatmaps", "Consensus", THR_FOLDER[thr])
    os.makedirs(leg_dir, exist_ok=True)
    with open(os.path.join(leg_dir, "group_legend.txt"), "w") as f:
        f.write(f"Group legend — Top-{top_n} consensus heatmap ({THR_FOLDER[thr]})\n")
        f.write("Letter\tGroup\n")
        for l, c in zip(letters, cols):
            f.write(f"{l}\t{c}\n")
    fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.6), max(5, len(genes) * 0.55)))
    _v = np.abs(pivot.values); _v = _v[np.isfinite(_v)]
    vmax = max(1.5, float(np.nanpercentile(_v, 90))) if _v.size else 1.5
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(letters, fontsize=12)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(pivot.index, fontsize=11)
    ax.set_xlabel("Group (see group_legend.txt)", fontsize=12)
    ax.set_title(f"Top {top_n} consensus genes — mean log2FC by group ({THR_FOLDER[thr]})", fontsize=12)
    plt.colorbar(im, label="mean log2FC", shrink=0.7)
    save(fig, os.path.join("Heatmaps", "Consensus"), thr, "heatmap_top20_consensus")


def main():
    for thr in ["LFC02", "LFC05", "LFC1"]:
        enr = pd.read_csv(os.path.join(DATA, f"enrichment_{thr}.csv"))
        cons = pd.read_csv(os.path.join(DATA, f"consensus_{thr}.csv"))
        master = pd.read_csv(os.path.join(DATA, f"master_deg_{thr}.csv"))
        print(f"=== {thr} ===")
        for db in DATABASES:
            dotplot(enr, db, thr); barplot_db(enr, db, thr)
        top_gene_barplot(cons, thr)
        consensus_heatmap(cons, master, thr)
        print(f"  plots written for {thr}")


if __name__ == "__main__":
    main()
