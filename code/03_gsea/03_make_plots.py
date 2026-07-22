#!/usr/bin/env python3
"""
03_make_plots.py — all GSEA figures (saved to GSEA/ root):
  01 venn gene overlaps (by condition)      02 upset across accessions
  03 pathway heatmap Hallmark (NES)         04 pathway heatmap KEGG (NES)
  05 ORA dotplot Hallmark                   05b ORA dotplot KEGG
  06 cross-condition Hallmark (mean NES)    07 cross-condition KEGG
  08 GSEA bubble Hallmark (summary)
"""
import os, glob, json, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import upsetplot

ROOT = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/GSEA"
GSEA = os.path.join(ROOT, "data", "gsea")
ORA = os.path.join(ROOT, "data", "ora")
cat = pd.read_csv("/sessions/amazing-zen-bardeen/work/gsea_groups.csv")
COND = dict(zip(cat["group"], cat["condition"]))
ACC = dict(zip(cat["group"], cat["accession"]))
COND_ORDER = ["Obesity", "Aging", "TBI", "Neurodegen"]
sig = json.load(open(os.path.join(ROOT, "data", "sig_genes.json")))


def save(fig, name):
    fig.savefig(os.path.join(ROOT, name), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  saved", name)


def nes_matrix(lib):
    """pathways x groups NES matrix."""
    data = {}
    for f in glob.glob(os.path.join(GSEA, f"*__{lib}.csv")):
        g = os.path.basename(f).replace(f"__{lib}.csv", "")
        d = pd.read_csv(f)
        data[g] = dict(zip(d["Term"], pd.to_numeric(d["NES"], errors="coerce")))
    mat = pd.DataFrame(data)
    # order columns by condition then name
    cols = sorted(mat.columns, key=lambda g: (COND_ORDER.index(COND.get(g, "Neurodegen")), g))
    return mat[cols]


def clean(t):
    return re.sub(r"^HALLMARK_", "", str(t)).replace("_", " ").title()[:42]


def heatmap_nes(lib, fname, title, top=30):
    mat = nes_matrix(lib).fillna(0)
    var = mat.var(axis=1).sort_values(ascending=False)
    mat = mat.loc[var.head(top).index]
    mat.index = [clean(t) for t in mat.index]
    fig, ax = plt.subplots(figsize=(max(12, mat.shape[1] * 0.42), max(8, len(mat) * 0.32)))
    vmax = np.nanpercentile(np.abs(mat.values), 98) or 2
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index, fontsize=8)
    # condition colour bar along the top
    cond_col = {"Obesity": "#1f77b4", "Aging": "#ff7f0e", "TBI": "#2ca02c", "Neurodegen": "#d62728"}
    for j, g in enumerate(mat.columns):
        ax.add_patch(plt.Rectangle((j - 0.5, -1.2), 1, 0.8, color=cond_col[COND.get(g, "Neurodegen")], clip_on=False))
    ax.set_title(title, fontsize=13)
    plt.colorbar(im, label="NES", shrink=0.6)
    handles = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=10, label=k)
               for k, c in cond_col.items()]
    ax.legend(handles=handles, title="Condition", bbox_to_anchor=(1.12, 1), loc="upper left", fontsize=8)
    save(fig, fname)


def cross_condition(lib, fname, title, top=25):
    mat = nes_matrix(lib).fillna(0)
    cond_means = {}
    for cond in COND_ORDER:
        cols = [g for g in mat.columns if COND.get(g) == cond]
        if cols:
            cond_means[cond] = mat[cols].mean(axis=1)
    cm = pd.DataFrame(cond_means)
    var = cm.var(axis=1).sort_values(ascending=False)
    cm = cm.loc[var.head(top).index]; cm.index = [clean(t) for t in cm.index]
    fig, ax = plt.subplots(figsize=(7, max(8, len(cm) * 0.34)))
    vmax = np.nanpercentile(np.abs(cm.values), 98) or 2
    im = ax.imshow(cm.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(cm.shape[1])); ax.set_xticklabels(cm.columns, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(cm))); ax.set_yticklabels(cm.index, fontsize=8)
    for i in range(len(cm)):
        for j in range(cm.shape[1]):
            ax.text(j, i, f"{cm.values[i,j]:.1f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(cm.values[i, j]) > vmax * 0.6 else "black")
    ax.set_title(title, fontsize=13)
    plt.colorbar(im, label="mean NES", shrink=0.6)
    save(fig, fname)


def ora_dotplot(key, fname, title, top=18):
    d = pd.read_csv(os.path.join(ORA, f"consensus_ora_{key}.csv"))
    d = d[d["Adjusted P-value"] < 0.05].copy()
    if d.empty:
        d = pd.read_csv(os.path.join(ORA, f"consensus_ora_{key}.csv")).nsmallest(top, "Adjusted P-value")
    d["n"] = d["Overlap"].map(lambda x: int(str(x).split("/")[0]))
    d["nlp"] = -np.log10(d["Adjusted P-value"])
    d = d.nlargest(top, "Combined Score").sort_values("Combined Score")
    fig, ax = plt.subplots(figsize=(10, max(5, len(d) * 0.45)))
    sc = ax.scatter(d["Combined Score"], range(len(d)), s=d["n"] * 18, c=d["nlp"],
                    cmap="plasma", edgecolors="black", linewidths=0.5)
    ax.set_yticks(range(len(d))); ax.set_yticklabels([clean(t) for t in d["Term"]], fontsize=9)
    ax.set_xlabel("Combined Score"); ax.set_title(title, fontsize=12)
    plt.colorbar(sc, label="-log10(adj p)", shrink=0.7)
    h = [plt.scatter([], [], s=s*18, c="grey", label=str(s)) for s in [5, 15, 30]]
    ax.legend(handles=h, title="genes", loc="lower right", fontsize=8)
    save(fig, fname)


def gsea_bubble(lib, fname, title):
    """Hallmark bubble: mean NES per pathway across groups; size = # groups |NES|>1 & FDR<0.25."""
    files = glob.glob(os.path.join(GSEA, f"*__{lib}.csv"))
    nes, hits = {}, {}
    for f in files:
        d = pd.read_csv(f)
        d["NES"] = pd.to_numeric(d["NES"], errors="coerce")
        d["FDR q-val"] = pd.to_numeric(d["FDR q-val"], errors="coerce")
        for _, r in d.iterrows():
            nes.setdefault(r["Term"], []).append(r["NES"])
            if abs(r["NES"]) > 1 and r["FDR q-val"] < 0.25:
                hits[r["Term"]] = hits.get(r["Term"], 0) + 1
    rows = [(t, np.nanmean(v), hits.get(t, 0)) for t, v in nes.items()]
    df = pd.DataFrame(rows, columns=["Term", "meanNES", "nhit"]).sort_values("meanNES")
    df = pd.concat([df.head(12), df.tail(12)]).drop_duplicates("Term")
    fig, ax = plt.subplots(figsize=(9, max(7, len(df) * 0.36)))
    vmax = max(1, np.nanpercentile(np.abs(df["meanNES"]), 98))
    sc = ax.scatter(df["meanNES"], range(len(df)), s=df["nhit"] * 22 + 25,
                    c=df["meanNES"], cmap="RdBu_r", vmin=-vmax, vmax=vmax, edgecolors="black", linewidths=0.5)
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    ax.set_yticks(range(len(df))); ax.set_yticklabels([clean(t) for t in df["Term"]], fontsize=9)
    ax.set_xlabel("mean NES across groups"); ax.set_title(title, fontsize=12)
    plt.colorbar(sc, label="mean NES", shrink=0.7)
    h = [plt.scatter([], [], s=n*22+25, c="grey", label=str(n)) for n in [1, 5, 10]]
    ax.legend(handles=h, title="# groups sig", loc="lower right", fontsize=8)
    save(fig, fname)


def venn_overlaps():
    sets = {c: set() for c in ["Obesity", "Aging", "TBI"]}
    for g, genes in sig.items():
        c = COND.get(g)
        if c in sets:
            sets[c].update(genes)
    fig, ax = plt.subplots(figsize=(8, 7))
    venn3([sets["Obesity"], sets["Aging"], sets["TBI"]],
          set_labels=("Obesity", "Aging", "TBI"), ax=ax)
    ax.set_title("DEG overlap across conditions (LFC0.2)", fontsize=13)
    save(fig, "01_venn_gene_overlaps.png")


def upset_accessions():
    contents = {}
    for g, genes in sig.items():
        a = f"Acc{ACC[g]}"
        contents.setdefault(a, set()).update(genes)
    data = upsetplot.from_contents(contents)
    fig = plt.figure(figsize=(13, 7))
    upsetplot.UpSet(data, subset_size="count", show_counts=True, sort_by="cardinality",
                    min_subset_size=8).plot(fig=fig)
    fig.suptitle("DEG set intersections across accessions", fontsize=13)
    save(fig, "02_upset_all_accessions.png")


if __name__ == "__main__":
    print("plots:")
    venn_overlaps()
    upset_accessions()
    heatmap_nes("hallmark", "03_pathway_heatmap_hallmark.png", "GSEA NES — Hallmark (top-variable pathways x groups)")
    heatmap_nes("kegg", "04_pathway_heatmap_kegg.png", "GSEA NES — KEGG (top-variable pathways x groups)")
    ora_dotplot("hallmark", "05_ora_dotplot_hallmark.png", "ORA — Hallmark (consensus DEGs)")
    ora_dotplot("kegg", "05b_ora_dotplot_kegg.png", "ORA — KEGG (consensus DEGs)")
    cross_condition("hallmark", "06_cross_condition_hallmark.png", "Cross-condition mean NES — Hallmark")
    cross_condition("kegg", "07_cross_condition_kegg.png", "Cross-condition mean NES — KEGG")
    gsea_bubble("hallmark", "08_gsea_bubble_hallmark.png", "GSEA summary — Hallmark (mean NES across groups)")
    print("DONE")
