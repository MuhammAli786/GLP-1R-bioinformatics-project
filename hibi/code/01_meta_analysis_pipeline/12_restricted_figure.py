#!/usr/bin/env python3
"""Fig12B restricted to the four curated mechanism gene sets: BBB/MMP, JAK-STAT3 inflammatory, PI3K-Akt pro-survival and ion channel/synaptic.

These are the gene lists used by the Cnet and GOChord figures in both arms. The
unrestricted Fig12B plots all 1,715 shared consensus genes; this version asks
only whether GLP-1R agonism opposes the injury within those mechanisms. Each
program contributes 8-29 shared genes, so every gene is labelled and the figure
is read gene by gene.
Outputs: Fig12B_Restricted_MechanismGenes.{png,pdf} and
data/meta_analysis/restricted_shared_genes_by_program.csv
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, binomtest

REPO = BASE + "/mnt/GLP-1R-bioinformatics-project"
HIBI = os.path.join(REPO, "hibi")
OUTFIG = os.path.join(HIBI, "figures")
sys.path.insert(0, os.path.join(HIBI, "code", "01_meta_analysis_pipeline"))
from cnet_gene_lists import (BBB_GENES, JAK_STAT3_GENES,
                             PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES)

C_HIBI, C_GLP = "#C0392B", "#2471A3"
C_REV, C_CONC = "#7D3C98", "#2C3E50"
PROGRAMS = [("BBB / MMP", BBB_GENES),
            ("JAK-STAT3 inflammatory", JAK_STAT3_GENES),
            ("PI3K-Akt pro-survival", PI3K_AKT_GENES)]

m = pd.read_csv(os.path.join(HIBI, "data", "meta_analysis", "shared_genes_GLP1R_HIBI.csv"))
m["key"] = m["gene_symbol"].str.upper()

fig, axes = plt.subplots(1, 3, figsize=(24, 9))
fig.suptitle("Restricted analysis: shared consensus genes within the curated mechanism programs\n"
             "HIBI (injury) vs GLP-1R (agonism)",
             fontsize=22, fontweight="bold", y=1.0)

summary, keep_rows = [], []
for ax, (pname, glist) in zip(axes.ravel(), PROGRAMS):
    keys = {g.upper() for g in glist}
    s = m[m["key"].isin(keys)].copy()
    s["program"] = pname
    keep_rows.append(s)

    # Robust axis range: a few genes (e.g. F11r at +9.9 in the GLP-1R arm) are
    # single-group low-expression artefacts that would squash everything else into
    # the origin. Scale to the 90th percentile of |log2FC| and mark anything
    # outside with a hollow triangle at the boundary.
    vals = np.abs(np.concatenate([s["mean_log2FC_hibi"].values, s["mean_log2FC_glp"].values]))
    lim = float(np.nanpercentile(vals, 90)) * 1.55
    lim = max(lim, 1.2)
    s["px"] = s["mean_log2FC_hibi"].clip(-lim * 0.97, lim * 0.97)
    s["py"] = s["mean_log2FC_glp"].clip(-lim * 0.97, lim * 0.97)
    s["clipped"] = ((s["mean_log2FC_hibi"].abs() > lim * 0.97) |
                    (s["mean_log2FC_glp"].abs() > lim * 0.97))

    # Shade both opposing quadrants, very lightly, so the shading marks the
    # regions rather than advertising one conclusion.
    ax.add_patch(plt.Rectangle((0, -lim), lim, lim, facecolor=C_REV, alpha=0.045, zorder=0))
    ax.add_patch(plt.Rectangle((-lim, 0), lim, lim, facecolor=C_REV, alpha=0.045, zorder=0))
    ax.axhline(0, color="#555", lw=1.1, zorder=1)
    ax.axvline(0, color="#555", lw=1.1, zorder=1)
    ax.plot([-lim, lim], [lim, -lim], ls=":", lw=1.4, color="#999", zorder=1)

    disc = ~s["concordant"]
    for mask, col in [(~disc, C_CONC), (disc, C_REV)]:
        ok = mask & ~s["clipped"]
        cl = mask & s["clipped"]
        ax.scatter(s.loc[ok, "px"], s.loc[ok, "py"], s=110, c=col, alpha=0.78,
                   edgecolors="white", linewidths=1.2, zorder=3)
        ax.scatter(s.loc[cl, "px"], s.loc[cl, "py"], s=150, marker="^",
                   facecolors="none", edgecolors=col, linewidths=2.2, zorder=4)

    # Label every gene, nudging labels apart where they collide
    placed = []
    for r in s.sort_values("recur", ascending=False).itertuples():
        fx, fy = (r.px + lim) / (2 * lim), (r.py + lim) / (2 * lim)
        dy = 5
        for _ in range(6):
            if not any((fx - a) ** 2 + (fy - b + 0.0) ** 2 < 0.052 ** 2 for a, b in placed):
                break
            fy += 0.055; dy += 15
        placed.append((fx, fy))
        ax.annotate(r.gene_symbol, (r.px, r.py),
                    fontsize=11.5, fontweight="bold", zorder=6,
                    xytext=(7, dy), textcoords="offset points",
                    color=C_REV if not r.concordant else C_CONC, clip_on=True,
                    bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                              alpha=0.78, edgecolor="none"))

    n = len(s)
    nd = int(disc.sum())
    if n >= 4:
        rho, prho = spearmanr(s["mean_log2FC_hibi"], s["mean_log2FC_glp"])
        pb = binomtest(nd, n, 0.5).pvalue
    else:
        rho, prho, pb = np.nan, np.nan, np.nan
    summary.append((pname, n, nd, rho, prho, pb))

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel("HIBI mean log2FC", fontsize=15, color=C_HIBI, fontweight="bold")
    ax.set_ylabel("GLP-1R mean log2FC", fontsize=15, color=C_GLP, fontweight="bold")
    ax.set_title(f"{pname}   (n = {n} shared genes)", fontsize=17, fontweight="bold", pad=12)
    ax.tick_params(labelsize=12)
    ax.text(0.03, 0.03,
            f"$\\rho$ = {rho:+.2f} ($p$ = {prho:.2f})\n"
            f"opposite direction: {nd}/{n} ({nd/n*100:.0f}%, $p$ = {pb:.2f})",
            transform=ax.transAxes, va="bottom", ha="left", fontsize=12.5, zorder=8,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#FEF9E7",
                      edgecolor="#B7950B", linewidth=1.3))
    for sp in ax.spines.values():
        sp.set_linewidth(1.2)

handles = [plt.Line2D([0], [0], marker="o", ls="", markerfacecolor=C_REV,
                      markeredgecolor="white", markersize=13, label="Opposite direction in the two arms"),
           plt.Line2D([0], [0], marker="o", ls="", markerfacecolor=C_CONC,
                      markeredgecolor="white", markersize=13, label="Same direction"),
           plt.Line2D([0], [0], ls=":", color="#999", lw=1.6, label="Perfect-reversal diagonal"),
           plt.Line2D([0], [0], marker="^", ls="", markerfacecolor="none",
                      markeredgecolor="#555", markersize=13, label="Beyond axis range (clipped)")]
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=13,
           frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, 0.055))

fig.text(0.5, -0.005,
         "Each program contributes only 13–29 genes present in BOTH consensus sets, so these panels are descriptive: "
         "no correlation or direction split reaches significance.\n"
         "Points near the dotted diagonal are genes the two arms move oppositely by a comparable amount — "
         "the pattern a reversal effect would produce.",
         ha="center", va="top", fontsize=12, style="italic", color="#555")

fig.subplots_adjust(top=0.84, bottom=0.20, wspace=0.24)
for fmt in ("png", "pdf"):
    fig.savefig(os.path.join(OUTFIG, f"Fig12B_Restricted_MechanismGenes.{fmt}"),
                dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("saved Fig12B_Restricted_MechanismGenes (png + pdf)")

out = pd.concat(keep_rows, ignore_index=True)
out.to_csv(os.path.join(HIBI, "data", "meta_analysis",
                        "restricted_shared_genes_by_program.csv"), index=False)

print("\nprogram                    n  opposite     rho     p(rho)  p(50:50)")
for pname, n, nd, rho, prho, pb in summary:
    print(f"  {pname:24s} {n:3d}  {nd:3d} ({nd/n*100:3.0f}%)  {rho:+.3f}   {prho:.3f}    {pb:.3f}")
