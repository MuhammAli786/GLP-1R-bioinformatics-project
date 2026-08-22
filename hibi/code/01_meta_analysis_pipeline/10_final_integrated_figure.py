#!/usr/bin/env python3
"""Integrated summary figure of shared genes, pathways and signalling between the GLP-1R agonism arm and the HIBI (neonatal hypoxic-ischemic brain injury) arm.

Panel A: consensus gene overlap (Venn), shared set split concordant/discordant.
Panel B: shared-gene log2FC scatter, HIBI (x) vs GLP-1R (y); off-diagonal
quadrants hold genes the injury and the drug move in opposite directions, and the
top reversals are labelled.
Panel C: pathways enriched in both arms, -log10 adjusted p side by side.
Panel D: per-gene log2FC in each arm for the four curated mechanism programs also
used by the GLP-1R Cnet/GOChord figures.
Inputs: GLP-1R data/meta_analysis/consensus_LFC02.csv and enrichment_LFC02.csv;
HIBI hibi/data/meta_analysis/consensus_LFC02.csv and
hibi/data/goplot/consensus_{terms_FULL,kegg_terms,reactome_terms}.csv.
Outputs: hibi/figures/Fig12_Integrated_GLP1R_vs_HIBI.{png,pdf} and
hibi/data/meta_analysis/shared_pathways_GLP1R_HIBI.csv
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
from matplotlib.patches import Patch
from matplotlib_venn import venn2

REPO = BASE + "/mnt/GLP-1R-bioinformatics-project"
GLP_DATA = os.path.join(REPO, "data", "meta_analysis")
HIBI_DATA = os.path.join(REPO, "hibi", "data")
OUTFIG = os.path.join(REPO, "hibi", "figures")
sys.path.insert(0, os.path.join(REPO, "hibi", "code", "01_meta_analysis_pipeline"))
from cnet_gene_lists import (BBB_GENES, JAK_STAT3_GENES,
                             PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES)

VMAX = 2.0
C_HIBI = "#C0392B"     # injury  (red)
C_GLP = "#2471A3"      # drug    (blue)
C_REV = "#7D3C98"      # reversal / discordant (purple)
C_CONC = "#7F8C8D"     # concordant (grey)


def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", str(t))
    t = re.sub(r"\s*R-HSA-\d+", "", t)
    return t.strip()


# Load consensus sets
glp = pd.read_csv(os.path.join(GLP_DATA, "consensus_LFC02.csv"))
hib = pd.read_csv(os.path.join(HIBI_DATA, "meta_analysis", "consensus_LFC02.csv"))
glp["key"] = glp["gene_symbol"].str.upper()
hib["key"] = hib["gene_symbol"].str.upper()

m = hib.merge(glp, on="key", suffixes=("_hibi", "_glp"))
m["concordant"] = m["predominant_direction_hibi"] == m["predominant_direction_glp"]
n_shared, n_conc, n_disc = len(m), int(m["concordant"].sum()), int((~m["concordant"]).sum())
print(f"HIBI {len(hib)}  GLP-1R {len(glp)}  shared {n_shared} "
      f"(concordant {n_conc} / discordant {n_disc})")

# Direction comes from predominant_direction (the per-group UP/DOWN vote), not
# from sign(mean_log2FC); the two disagree for ~350 genes per arm because a single
# extreme group can flip the mean, and the vote is the robust call used for the
# concordant/discordant counts above.
# Reversals are ranked by recurrence (groups supporting the call in each arm), not
# by |log2FC| gap: a few genes carry |mean log2FC| > 10 and are low-expression or
# single-group artefacts that would otherwise dominate.
m["recur"] = m["n_groups_hibi"] + m["n_groups_glp"]
m["reversal"] = m["recur"]

# Clip plotted values to the panel range so extreme genes sit at the axis edge
# instead of compressing everything else.
PLIM_X, PLIM_Y = 5.0, 2.2
m["plot_x"] = m["mean_log2FC_hibi"].clip(-PLIM_X, PLIM_X)
m["plot_y"] = m["mean_log2FC_glp"].clip(-PLIM_Y, PLIM_Y)
n_clip = int(((m["mean_log2FC_hibi"].abs() > PLIM_X) |
              (m["mean_log2FC_glp"].abs() > PLIM_Y)).sum())
print(f"{n_clip} of {n_shared} shared genes clipped to panel-B axis limits")

# Shared pathways, significant in both arms
ge = pd.read_csv(os.path.join(GLP_DATA, "enrichment_LFC02.csv"))
ge = ge[ge["Adjusted P-value"] < 0.05].copy()
ge["term_key"] = ge["Term"].map(clean_term).str.lower()
ge = ge.sort_values("Adjusted P-value").drop_duplicates("term_key")

hfiles = ["consensus_terms_FULL.csv", "consensus_kegg_terms.csv", "consensus_reactome_terms.csv"]
he = pd.concat([pd.read_csv(os.path.join(HIBI_DATA, "goplot", f))
                for f in hfiles if os.path.exists(os.path.join(HIBI_DATA, "goplot", f))],
               ignore_index=True)
he = he[he["adj_pval"] < 0.05].copy()
he["term_key"] = he["Term"].map(clean_term).str.lower()
he = he.sort_values("adj_pval").drop_duplicates("term_key")

sp = he.merge(ge[["term_key", "Term", "Adjusted P-value", "Gene_set"]],
              on="term_key", suffixes=("", "_glp"))
sp = sp.rename(columns={"adj_pval": "padj_hibi", "Adjusted P-value": "padj_glp"})
sp["nlp_hibi"] = -np.log10(sp["padj_hibi"].clip(lower=1e-300))
sp["nlp_glp"] = -np.log10(sp["padj_glp"].clip(lower=1e-300))
sp["combined"] = sp["nlp_hibi"] + sp["nlp_glp"]
sp = sp.sort_values("combined", ascending=False)
print(f"pathways significant in BOTH arms: {len(sp)}")

sp_out = sp[["Category", "Term", "padj_hibi", "padj_glp", "nlp_hibi", "nlp_glp", "combined"]]
sp_out.to_csv(os.path.join(HIBI_DATA, "meta_analysis", "shared_pathways_GLP1R_HIBI.csv"), index=False)

# Figure
fig = plt.figure(figsize=(22, 19))
gs = fig.add_gridspec(2, 2, hspace=0.26, wspace=0.30,
                      left=0.075, right=0.955, top=0.90, bottom=0.05)

fig.suptitle("Convergence of GLP-1R agonism and neonatal hypoxic-ischemic brain injury\n"
             "consensus transcriptomic programs  (padj < 0.05, |log2FC| $\\geq$ 0.2, gene significant in $\\geq$2 groups)",
             fontsize=25, fontweight="bold", y=0.965)

# Panel A: Venn
axA = fig.add_subplot(gs[0, 0])
v = venn2(subsets=(len(hib) - n_shared, len(glp) - n_shared, n_shared),
          set_labels=("HIBI\n(injury)", "GLP-1R\n(agonism)"), ax=axA)
for pid, col in [("10", C_HIBI), ("01", C_GLP), ("11", C_REV)]:
    p = v.get_patch_by_id(pid)
    if p:
        p.set_color(col); p.set_alpha(0.68); p.set_edgecolor("white"); p.set_linewidth(2.5)
for pid in ("10", "01", "11"):
    lb = v.get_label_by_id(pid)
    if lb:
        lb.set_fontsize(21); lb.set_fontweight("bold"); lb.set_color("white")
for t, col in zip(v.set_labels, (C_HIBI, C_GLP)):
    if t:
        t.set_fontsize(20); t.set_fontweight("bold"); t.set_color(col)
axA.text(0.5, -0.10,
         f"{n_shared:,} shared consensus genes\n"
         f"{n_conc:,} same direction   •   {n_disc:,} OPPOSITE direction",
         transform=axA.transAxes, ha="center", va="top", fontsize=18,
         bbox=dict(boxstyle="round,pad=0.55", facecolor="#F4ECF7",
                   edgecolor=C_REV, linewidth=2))
axA.set_title("A   Consensus gene overlap", fontsize=23, fontweight="bold", loc="left", pad=16)

# Panel B: log2FC scatter
axB = fig.add_subplot(gs[0, 1])
x, y = m["plot_x"], m["plot_y"]
opp = ~m["concordant"]          # same definition as Panel A
axB.axhline(0, color="#555", lw=1.2, zorder=1)
axB.axvline(0, color="#555", lw=1.2, zorder=1)
lim = PLIM_X
axB.add_patch(plt.Rectangle((0, -PLIM_Y), lim, PLIM_Y, facecolor=C_REV, alpha=0.07, zorder=0))
axB.add_patch(plt.Rectangle((-lim, 0), lim, PLIM_Y, facecolor=C_REV, alpha=0.07, zorder=0))
axB.scatter(x[~opp], y[~opp], s=26, c=C_CONC, alpha=0.45,
            edgecolors="none", zorder=2, label=f"Same direction (n={n_conc:,})")
axB.scatter(x[opp], y[opp], s=34, c=C_REV, alpha=0.60,
            edgecolors="none", zorder=3, label=f"Opposite direction (n={n_disc:,})")

# Label the most recurrent reversals, those supported by the most groups in both
# arms. Only points well inside the axes are labelled, and spacing is measured in
# axis-fraction units so it behaves the same on both axes.
PAD_X, PAD_Y = 0.35, 0.16
MIN_SEP = 0.052          # fraction of the panel
lab = m[opp].nlargest(200, "recur").copy()
placed = []
for r in lab.itertuples():
    xi, yi = r.plot_x, r.plot_y
    if not (-PLIM_X + PAD_X < xi < PLIM_X - PAD_X and -PLIM_Y + PAD_Y < yi < PLIM_Y - PAD_Y):
        continue
    fx, fy = (xi + PLIM_X) / (2 * PLIM_X), (yi + PLIM_Y) / (2 * PLIM_Y)
    if any((fx - px) ** 2 + (fy - py) ** 2 < MIN_SEP ** 2 for px, py in placed):
        continue
    placed.append((fx, fy))
    dx = 7 if xi < PLIM_X - 1.4 else -7
    axB.annotate(r.gene_symbol_hibi, (xi, yi), fontsize=13.5, fontweight="bold",
                 xytext=(dx, 5), textcoords="offset points", zorder=5,
                 color="#4A235A", clip_on=True,
                 ha="left" if dx > 0 else "right",
                 bbox=dict(boxstyle="round,pad=0.14", facecolor="white",
                           alpha=0.72, edgecolor="none"))
    if len(placed) >= 22:
        break
axB.set_xlabel("HIBI mean log2FC   (injury)", fontsize=19, color=C_HIBI, fontweight="bold")
axB.set_ylabel("GLP-1R mean log2FC   (agonism)", fontsize=19, color=C_GLP, fontweight="bold")
axB.set_xlim(-lim, lim); axB.set_ylim(-PLIM_Y, PLIM_Y)
axB.tick_params(labelsize=15)
axB.legend(fontsize=15, loc="upper left", framealpha=0.92)
axB.text(0.985, 0.03, "shaded quadrants =\ninjury and drug oppose",
         transform=axB.transAxes, ha="right", va="bottom", fontsize=14,
         style="italic", color=C_REV)
axB.set_title("B   Shared genes: direction in each arm", fontsize=23, fontweight="bold", loc="left", pad=16)
for s in axB.spines.values():
    s.set_linewidth(1.3)

# Panel C: shared pathways
axC = fig.add_subplot(gs[1, 0])
top = sp.head(16).iloc[::-1]
yy = np.arange(len(top))
h = 0.40
axC.barh(yy + h / 2, top["nlp_hibi"], height=h, color=C_HIBI, alpha=0.88, label="HIBI (injury)")
axC.barh(yy - h / 2, top["nlp_glp"], height=h, color=C_GLP, alpha=0.88, label="GLP-1R (agonism)")
axC.set_yticks(yy)
axC.set_yticklabels(["\n".join(textwrap.wrap(t, 34)) for t in top["Term"]], fontsize=13.5)
axC.set_xlabel("$-\\log_{10}$ adjusted $p$", fontsize=19, fontweight="bold")
axC.axvline(-np.log10(0.05), color="#333", ls="--", lw=1.6)
axC.text(-np.log10(0.05), len(top) - 0.3, " padj = 0.05", fontsize=13, color="#333", rotation=90,
         va="top", ha="left")
axC.legend(fontsize=15, loc="lower right", framealpha=0.92)
axC.set_title(f"C   Pathways enriched in BOTH arms  ({len(sp)} total)",
              fontsize=23, fontweight="bold", loc="left", pad=16)
axC.tick_params(axis="x", labelsize=15)
for s in axC.spines.values():
    s.set_linewidth(1.3)

# Panel D: curated mechanism programs
axD = fig.add_subplot(gs[1, 1])
PROGRAMS = [("BBB / MMP", BBB_GENES), ("JAK-STAT3\ninflammatory", JAK_STAT3_GENES),
            ("PI3K-Akt\npro-survival", PI3K_AKT_GENES)]
lfc_h = dict(zip(m["key"], m["mean_log2FC_hibi"]))
lfc_g = dict(zip(m["key"], m["mean_log2FC_glp"]))
recur = dict(zip(m["key"], m["recur"]))
disc = dict(zip(m["key"], ~m["concordant"]))

rows, ylabels, boundaries, prog_centres = [], [], [], []
cursor = 0
used = set()   # a gene can sit in several curated lists (e.g. Vegfa); show it once
for pname, glist in PROGRAMS:
    sel = [g for g in dict.fromkeys(glist)
           if g.upper() in lfc_h and g.upper() not in used]
    # most recurrent genes first, discordant (reversed) ones prioritised
    sel = sorted(sel, key=lambda g: (disc[g.upper()], recur[g.upper()]), reverse=True)[:11]
    sel = sorted(sel, key=lambda g: lfc_h[g.upper()], reverse=True)
    used.update(g.upper() for g in sel)
    for g in sel:
        rows.append([lfc_h[g.upper()], lfc_g[g.upper()]])
        ylabels.append(g)
    prog_centres.append((cursor + len(sel) / 2 - 0.5, pname, len(sel)))
    cursor += len(sel)
    boundaries.append(cursor - 0.5)

M = np.clip(np.array(rows), -VMAX, VMAX)
im = axD.imshow(M, cmap="RdBu_r", vmin=-VMAX, vmax=VMAX, aspect="auto")
axD.set_xticks([0, 1])
axD.set_xticklabels(["HIBI\n(injury)", "GLP-1R\n(agonism)"], fontsize=17, fontweight="bold")
for t, c in zip(axD.get_xticklabels(), (C_HIBI, C_GLP)):
    t.set_color(c)
axD.set_yticks(range(len(ylabels)))
axD.set_yticklabels(ylabels, fontsize=12.5)
for b in boundaries[:-1]:
    axD.axhline(b, color="black", lw=2.4)
axD.axvline(0.5, color="white", lw=2.4)
for centre, pname, nsel in prog_centres:
    axD.text(1.72, centre, pname, fontsize=15, fontweight="bold", va="center", ha="center")
cb = fig.colorbar(im, ax=axD, shrink=0.62, pad=0.20)
cb.set_label("mean log2FC", fontsize=17)
cb.ax.tick_params(labelsize=13)
axD.set_title("D   Shared mechanism programs, gene-level", fontsize=23, fontweight="bold", loc="left", pad=16)
axD.set_xlim(-0.5, 1.5)

os.makedirs(OUTFIG, exist_ok=True)
for fmt in ("png", "pdf"):
    fig.savefig(os.path.join(OUTFIG, f"Fig12_Integrated_GLP1R_vs_HIBI.{fmt}"),
                dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("saved Fig12_Integrated_GLP1R_vs_HIBI (png + pdf)")

# Text summary for the caption
print("\nTop 12 shared pathways (both arms):")
for r in sp.head(12).itertuples():
    print(f"  {r.Category:9s} {r.Term[:58]:58s} HIBI padj={r.padj_hibi:.2e}  GLP padj={r.padj_glp:.2e}")
print("\nTop 15 reversals (opposite direction, most recurrent across groups):")
for r in m[opp].nlargest(15, "recur").itertuples():
    print(f"  {r.gene_symbol_hibi:10s} HIBI {r.predominant_direction_hibi:4s} "
          f"({r.n_groups_hibi:2d} grp, {r.mean_log2FC_hibi:+.2f})   "
          f"GLP-1R {r.predominant_direction_glp:4s} ({r.n_groups_glp:2d} grp, {r.mean_log2FC_glp:+.2f})")

m_out = m[["gene_symbol_hibi", "n_groups_hibi", "predominant_direction_hibi", "mean_log2FC_hibi",
           "n_groups_glp", "predominant_direction_glp", "mean_log2FC_glp", "concordant", "recur"]]
m_out.rename(columns={"gene_symbol_hibi": "gene_symbol"}).sort_values(
    ["concordant", "recur"], ascending=[True, False]).to_csv(
    os.path.join(HIBI_DATA, "meta_analysis", "shared_genes_GLP1R_HIBI.csv"), index=False)
print("\nwrote shared_pathways_GLP1R_HIBI.csv and shared_genes_GLP1R_HIBI.csv")
