#!/usr/bin/env python3
"""The four integrated GLP-1R vs HIBI panels, each rendered as a standalone figure.

Outputs: Fig12A_Consensus_Overlap_Venn (consensus gene overlap),
Fig12B_Shared_Gene_Direction (shared-gene log2FC, HIBI vs GLP-1R),
Fig12C_Shared_Pathways (pathways enriched in both arms) and
Fig12D_Mechanism_Programs (curated programs at gene level).
Panel B differs from the composite version deliberately: opposing quadrants are
not shaded, because shading them while labelling only discordant genes
manufactured an antagonism the statistics do not support; gene labels are drawn
from the most recurrent shared genes regardless of direction; the measured
relationship is drawn on the panel as binned medians plus the Spearman
correlation and its p-value; and both axes use the same range.
Measured relationship: Spearman rho = +0.055 (p = 0.02) across 1,715 shared
genes, i.e. no systematic antagonism. Genes strongly induced by HI have mean
GLP-1R log2FC of +0.31, not negative. 840/1715 discordant = 49%, which is chance
(binomial p = 0.41).
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, re, sys, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import spearmanr

REPO = BASE + "/mnt/GLP-1R-bioinformatics-project"
GLP_DATA = os.path.join(REPO, "data", "meta_analysis")
HIBI_DATA = os.path.join(REPO, "hibi", "data")
OUTFIG = os.path.join(REPO, "hibi", "figures")
sys.path.insert(0, os.path.join(REPO, "hibi", "code", "01_meta_analysis_pipeline"))
from cnet_gene_lists import (BBB_GENES, JAK_STAT3_GENES,
                             PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES)

VMAX = 2.0
C_HIBI, C_GLP = "#C0392B", "#2471A3"
C_REV, C_CONC = "#7D3C98", "#7F8C8D"
os.makedirs(OUTFIG, exist_ok=True)


def save(fig, name):
    for fmt in ("png", "pdf"):
        fig.savefig(os.path.join(OUTFIG, f"{name}.{fmt}"),
                    dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {name} (png + pdf)")


def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", str(t))
    return re.sub(r"\s*R-HSA-\d+", "", t).strip()


# Data
glp = pd.read_csv(os.path.join(GLP_DATA, "consensus_LFC02.csv"))
hib = pd.read_csv(os.path.join(HIBI_DATA, "meta_analysis", "consensus_LFC02.csv"))
glp["key"] = glp["gene_symbol"].str.upper()
hib["key"] = hib["gene_symbol"].str.upper()
m = hib.merge(glp, on="key", suffixes=("_hibi", "_glp"))
m["concordant"] = m["predominant_direction_hibi"] == m["predominant_direction_glp"]
m["recur"] = m["n_groups_hibi"] + m["n_groups_glp"]
n_shared = len(m)
n_conc = int(m["concordant"].sum())
n_disc = n_shared - n_conc

LIM = 4.0
m["px"] = m["mean_log2FC_hibi"].clip(-LIM, LIM)
m["py"] = m["mean_log2FC_glp"].clip(-LIM, LIM)

rho, prho = spearmanr(m["mean_log2FC_hibi"], m["mean_log2FC_glp"])
print(f"shared {n_shared}  concordant {n_conc}  discordant {n_disc}  "
      f"({n_disc/n_shared*100:.1f}%)   Spearman rho={rho:+.3f} p={prho:.3g}")

# Fig A
fig, ax = plt.subplots(figsize=(10, 9))
v = venn2(subsets=(len(hib) - n_shared, len(glp) - n_shared, n_shared),
          set_labels=("HIBI\n(injury)", "GLP-1R\n(agonism)"), ax=ax)
for pid, col in [("10", C_HIBI), ("01", C_GLP), ("11", C_REV)]:
    p = v.get_patch_by_id(pid)
    if p:
        p.set_color(col); p.set_alpha(0.68)
        p.set_edgecolor("white"); p.set_linewidth(2.5)
for pid in ("10", "01", "11"):
    lb = v.get_label_by_id(pid)
    if lb:
        lb.set_fontsize(22); lb.set_fontweight("bold"); lb.set_color("white")
for t, col in zip(v.set_labels, (C_HIBI, C_GLP)):
    if t:
        t.set_fontsize(19); t.set_fontweight("bold"); t.set_color(col)
ax.set_title("Consensus gene overlap:\nGLP-1R agonism vs neonatal HI brain injury",
             fontsize=19, fontweight="bold", pad=18)
ax.text(0.5, -0.14,
        f"{n_shared:,} shared consensus genes\n"
        f"{n_conc:,} same direction  ({n_conc/n_shared*100:.0f}%)   •   "
        f"{n_disc:,} opposite  ({n_disc/n_shared*100:.0f}%)",
        transform=ax.transAxes, ha="center", va="top", fontsize=15,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F4ECF7",
                  edgecolor=C_REV, linewidth=1.8))
ax.text(0.5, -0.30,
        "A ~50/50 direction split is what independent gene sets produce\n"
        "(binomial vs 50:50, $p$ = 0.41): the split alone is not evidence of reversal.",
        transform=ax.transAxes, ha="center", va="top", fontsize=12.5,
        style="italic", color="#555")
fig.subplots_adjust(bottom=0.22)
save(fig, "Fig12A_Consensus_Overlap_Venn")

# Fig B
fig, ax = plt.subplots(figsize=(11.5, 10.5))
opp = ~m["concordant"]
ax.axhline(0, color="#555", lw=1.2, zorder=1)
ax.axvline(0, color="#555", lw=1.2, zorder=1)
ax.scatter(m["px"][~opp], m["py"][~opp], s=26, c=C_CONC, alpha=0.45,
           edgecolors="none", zorder=2, label=f"Same direction (n={n_conc:,})")
ax.scatter(m["px"][opp], m["py"][opp], s=30, c=C_REV, alpha=0.50,
           edgecolors="none", zorder=2, label=f"Opposite direction (n={n_disc:,})")

# Measured relationship: median GLP-1R log2FC in bins of HIBI log2FC
bins = np.array([-4, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 4])
cen, med = [], []
for lo, hi in zip(bins[:-1], bins[1:]):
    s = m[(m["mean_log2FC_hibi"] > lo) & (m["mean_log2FC_hibi"] <= hi)]
    if len(s) >= 15:
        cen.append((lo + hi) / 2); med.append(s["mean_log2FC_glp"].median())
ax.plot(cen, med, "-o", color="#E67E22", lw=3.2, ms=9, zorder=6,
        markeredgecolor="white", markeredgewidth=1.5,
        label="Median GLP-1R log2FC per bin")

# Label the most recurrent shared genes, direction-agnostic, to stay unbiased
PAD, MIN_SEP = 0.62, 0.062
placed = []
for r in m.nlargest(220, "recur").itertuples():
    xi, yi = r.px, r.py
    if not (-LIM + PAD < xi < LIM - PAD and -LIM + PAD < yi < LIM - PAD):
        continue
    fx, fy = (xi + LIM) / (2 * LIM), (yi + LIM) / (2 * LIM)
    if any((fx - a) ** 2 + (fy - b) ** 2 < MIN_SEP ** 2 for a, b in placed):
        continue
    placed.append((fx, fy))
    ax.annotate(r.gene_symbol_hibi, (xi, yi), fontsize=12, fontweight="bold",
                xytext=(6, 5), textcoords="offset points", zorder=7,
                color="#4A235A" if not r.concordant else "#2C3E50", clip_on=True,
                bbox=dict(boxstyle="round,pad=0.13", facecolor="white",
                          alpha=0.75, edgecolor="none"))
    if len(placed) >= 24:
        break

ax.set_xlabel("HIBI mean log2FC   (injury)", fontsize=17, color=C_HIBI, fontweight="bold")
ax.set_ylabel("GLP-1R mean log2FC   (agonism)", fontsize=17, color=C_GLP, fontweight="bold")
ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM)
ax.set_aspect("equal")
ax.tick_params(labelsize=14)
ax.legend(fontsize=13, loc="upper left", framealpha=0.93)
ax.set_title("Shared consensus genes: direction in each arm",
             fontsize=19, fontweight="bold", pad=16)
ax.text(0.975, 0.03,
        f"Spearman $\\rho$ = {rho:+.3f}  ($p$ = {prho:.2g}, n = {n_shared:,})\n"
        "No systematic antagonism: a reversal effect\nwould give a clearly NEGATIVE $\\rho$.",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=12.5,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#FEF9E7",
                  edgecolor="#B7950B", linewidth=1.5))
ax.text(0.02, -0.115,
        "Labels are the most recurrent shared genes regardless of direction "
        "(purple = opposite, dark = same), so the labelling does not favour either story.",
        transform=ax.transAxes, ha="left", va="top", fontsize=11.5, style="italic", color="#555")
for s in ax.spines.values():
    s.set_linewidth(1.3)
save(fig, "Fig12B_Shared_Gene_Direction")

# Fig C
ge = pd.read_csv(os.path.join(GLP_DATA, "enrichment_LFC02.csv"))
ge = ge[ge["Adjusted P-value"] < 0.05].copy()
ge["term_key"] = ge["Term"].map(clean_term).str.lower()
ge = ge.sort_values("Adjusted P-value").drop_duplicates("term_key")

hf = ["consensus_terms_FULL.csv", "consensus_kegg_terms.csv", "consensus_reactome_terms.csv"]
he = pd.concat([pd.read_csv(os.path.join(HIBI_DATA, "goplot", f))
                for f in hf if os.path.exists(os.path.join(HIBI_DATA, "goplot", f))],
               ignore_index=True)
he = he[he["adj_pval"] < 0.05].copy()
he["term_key"] = he["Term"].map(clean_term).str.lower()
he = he.sort_values("adj_pval").drop_duplicates("term_key")

sp = he.merge(ge[["term_key", "Adjusted P-value"]], on="term_key")
sp = sp.rename(columns={"adj_pval": "padj_hibi", "Adjusted P-value": "padj_glp"})
sp["nlp_hibi"] = -np.log10(sp["padj_hibi"].clip(lower=1e-300))
sp["nlp_glp"] = -np.log10(sp["padj_glp"].clip(lower=1e-300))
sp["combined"] = sp["nlp_hibi"] + sp["nlp_glp"]
sp = sp.sort_values("combined", ascending=False)

fig, ax = plt.subplots(figsize=(13, 11))
top = sp.head(18).iloc[::-1]
yy = np.arange(len(top)); h = 0.40
ax.barh(yy + h / 2, top["nlp_hibi"], height=h, color=C_HIBI, alpha=0.88, label="HIBI (injury)")
ax.barh(yy - h / 2, top["nlp_glp"], height=h, color=C_GLP, alpha=0.88, label="GLP-1R (agonism)")
ax.set_yticks(yy)
ax.set_yticklabels(["\n".join(textwrap.wrap(t, 38)) for t in top["Term"]], fontsize=13)
ax.set_xlabel("$-\\log_{10}$ adjusted $p$", fontsize=17, fontweight="bold")
ax.axvline(-np.log10(0.05), color="#333", ls="--", lw=1.6)
ax.text(-np.log10(0.05), len(top) - 0.4, " padj = 0.05", fontsize=12,
        color="#333", rotation=90, va="top", ha="left")
ax.legend(fontsize=14, loc="lower right", framealpha=0.93)
ax.set_title(f"Pathways enriched in BOTH arms  ({len(sp)} total)",
             fontsize=19, fontweight="bold", pad=16)
ax.tick_params(axis="x", labelsize=14)
ax.text(0.0, -0.085,
        "Both arms are brain tissue, so synaptic and neurodevelopmental terms are expected to\n"
        "enrich in both. Shared enrichment indicates the same programs are engaged, not a shared direction.",
        transform=ax.transAxes, ha="left", va="top", fontsize=11.5, style="italic", color="#555")
for s in ax.spines.values():
    s.set_linewidth(1.3)
save(fig, "Fig12C_Shared_Pathways")

# Fig D
PROGRAMS = [("BBB / MMP", BBB_GENES), ("JAK-STAT3\ninflammatory", JAK_STAT3_GENES),
            ("PI3K-Akt\npro-survival", PI3K_AKT_GENES)]
lfc_h = dict(zip(m["key"], m["mean_log2FC_hibi"]))
lfc_g = dict(zip(m["key"], m["mean_log2FC_glp"]))
recur = dict(zip(m["key"], m["recur"]))
disc = dict(zip(m["key"], ~m["concordant"]))

rows, ylabels, bounds, centres = [], [], [], []
cursor, used = 0, set()
for pname, glist in PROGRAMS:
    sel = [g for g in dict.fromkeys(glist) if g.upper() in lfc_h and g.upper() not in used]
    sel = sorted(sel, key=lambda g: (disc[g.upper()], recur[g.upper()]), reverse=True)[:11]
    sel = sorted(sel, key=lambda g: lfc_h[g.upper()], reverse=True)
    used.update(g.upper() for g in sel)
    for g in sel:
        rows.append([lfc_h[g.upper()], lfc_g[g.upper()]]); ylabels.append(g)
    centres.append((cursor + len(sel) / 2 - 0.5, pname))
    cursor += len(sel); bounds.append(cursor - 0.5)

fig, ax = plt.subplots(figsize=(9.5, 13))
M = np.clip(np.array(rows), -VMAX, VMAX)
im = ax.imshow(M, cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, aspect="auto")
ax.set_xticks([0, 1])
ax.set_xticklabels(["HIBI\n(injury)", "GLP-1R\n(agonism)"], fontsize=15, fontweight="bold")
for t, c in zip(ax.get_xticklabels(), (C_HIBI, C_GLP)):
    t.set_color(c)
ax.set_yticks(range(len(ylabels)))
ax.set_yticklabels(ylabels, fontsize=11.5)
for b in bounds[:-1]:
    ax.axhline(b, color="black", lw=2.4)
ax.axvline(0.5, color="white", lw=2.4)
for centre, pname in centres:
    ax.text(1.78, centre, pname, fontsize=13, fontweight="bold", va="center", ha="center")
cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.30)
cb.set_label("mean log2FC", fontsize=14)
cb.ax.tick_params(labelsize=12)
ax.set_title("Shared mechanism programs,\ngene level", fontsize=19, fontweight="bold", pad=16)
ax.set_xlim(-0.5, 1.5)
ax.text(0.0, -0.085,
        "Values clipped to ±2. Genes appearing in more than one curated list are shown once.",
        transform=ax.transAxes, ha="left", va="top", fontsize=10.5, style="italic", color="#555")
fig.subplots_adjust(bottom=0.10)
save(fig, "Fig12D_Mechanism_Programs")

print("\nall four standalone figures written to hibi/figures/")
