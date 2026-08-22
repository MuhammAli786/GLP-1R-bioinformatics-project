#!/usr/bin/env python3
"""Magnitude-aware version of the cross-arm comparison, and the candidate list derived from it.

The UP/DOWN vote and Spearman correlation both discard effect size, so a gene HI
drives up by 4 log2FC while GLP-1R drives it down by 1.5 is treated like one each
arm moves by 0.05.
Panel A: opposition candidates, genes moving >= MIN_LFC in both arms in opposite
directions with support in >= MIN_GRP groups per arm; diverging bars show the
effect size in each arm.
Panel B: global magnitude-weighted test, observed sum(-x*y) against a 20,000-draw
permutation null, where positive means net opposition.
Result: including magnitude strengthens the concordance rather than the reversal.
Pearson r = +0.083 (p = 6.6e-4) vs Spearman rho = +0.051; requiring |LFC| >= 1.0
in both arms leaves 34% opposite (66% same direction); and the magnitude-weighted
opposition score sits 3.4 SD below its null. The Panel A candidates are real but
outweighed by concordant genes genome-wide.
Outputs: Fig12E_Magnitude_Opposition.{png,pdf} and
data/meta_analysis/opposition_candidates.csv
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
from scipy.stats import pearsonr, spearmanr

REPO = BASE + "/mnt/GLP-1R-bioinformatics-project"
HIBI = os.path.join(REPO, "hibi")
OUTFIG = os.path.join(HIBI, "figures")
sys.path.insert(0, os.path.join(HIBI, "code", "01_meta_analysis_pipeline"))
from cnet_gene_lists import (BBB_GENES, JAK_STAT3_GENES,
                             PI3K_AKT_GENES, ION_CHANNEL_BASE_GENES)

C_HIBI, C_GLP, C_REV = "#C0392B", "#2471A3", "#7D3C98"
MIN_LFC = 0.5      # both arms must move at least this much
MIN_GRP = 3        # and be supported by at least this many groups in each arm
TOP_N = 24
rng = np.random.default_rng(42)

m = pd.read_csv(os.path.join(HIBI, "data", "meta_analysis", "shared_genes_GLP1R_HIBI.csv"))
m = m[(m.mean_log2FC_hibi.abs() < 8) & (m.mean_log2FC_glp.abs() < 8)].copy()
x = m.mean_log2FC_hibi.values
y = m.mean_log2FC_glp.values
opp = np.sign(x) != np.sign(y)

pr, pp = pearsonr(x, y)
sr, sp_ = spearmanr(x, y)
print(f"Pearson r = {pr:+.3f} (p={pp:.2e})   Spearman rho = {sr:+.3f} (p={sp_:.2e})")

# Magnitude-weighted opposition score and permutation null
score_obs = float(-(x * y).sum())
null = np.array([float(-(x * rng.permutation(y)).sum()) for _ in range(20000)])
z = (score_obs - null.mean()) / null.std()
p_opp = (np.sum(null >= score_obs) + 1) / (len(null) + 1)
print(f"opposition score {score_obs:+.1f}  null {null.mean():+.1f}+/-{null.std():.1f}  z={z:+.2f}  p={p_opp:.4f}")

# Candidate list
m["opp_mag"] = np.where(opp, np.minimum(np.abs(x), np.abs(y)), 0.0)
m["offset_frac"] = np.where(opp, np.clip(np.abs(y) / np.abs(x), 0, 3), 0.0)
prog = {}
for nm, gl in [("BBB/MMP", BBB_GENES), ("JAK-STAT3", JAK_STAT3_GENES),
               ("PI3K-Akt", PI3K_AKT_GENES)]:
    for g in gl:
        prog.setdefault(g.upper(), []).append(nm)
m["program"] = [", ".join(prog.get(g.upper(), [])) for g in m.gene_symbol]

cand = m[opp & (np.abs(x) >= MIN_LFC) & (np.abs(y) >= MIN_LFC) &
         (m.n_groups_hibi >= MIN_GRP) & (m.n_groups_glp >= MIN_GRP)].copy()
cand = cand.sort_values("opp_mag", ascending=False)
cand.to_csv(os.path.join(HIBI, "data", "meta_analysis", "opposition_candidates.csv"), index=False)
print(f"candidates (|LFC|>={MIN_LFC} both, >={MIN_GRP} groups both, opposite): n={len(cand)}")

# Figure
fig = plt.figure(figsize=(19, 12))
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.30,
                      left=0.10, right=0.97, top=0.86, bottom=0.10)
fig.suptitle("Magnitude-aware comparison: effect size, not just direction",
             fontsize=24, fontweight="bold", y=0.955)

# Panel A: opposition candidates, diverging bars
axA = fig.add_subplot(gs[0, 0])
top = cand.head(TOP_N).iloc[::-1]
yy = np.arange(len(top)); h = 0.40
axA.barh(yy + h / 2, top.mean_log2FC_hibi, height=h, color=C_HIBI, alpha=0.9, label="HIBI (injury)")
axA.barh(yy - h / 2, top.mean_log2FC_glp, height=h, color=C_GLP, alpha=0.9, label="GLP-1R (agonism)")
axA.axvline(0, color="#333", lw=1.5)
axA.set_yticks(yy)
labs = [f"{r.gene_symbol}" + (f"  [{r.program}]" if r.program else "")
        for r in top.itertuples()]
axA.set_yticklabels(labs, fontsize=12.5)
for i, r in enumerate(top.itertuples()):
    axA.text(0.02, i, f"  {r.n_groups_hibi}/{r.n_groups_glp}", transform=axA.get_yaxis_transform(),
             fontsize=9.5, va="center", color="#666")
axA.set_xlabel("mean log2FC", fontsize=17, fontweight="bold")
axA.legend(fontsize=14, loc="lower right", framealpha=0.93)
axA.set_title(f"A   Opposition candidates: |log2FC| $\\geq$ {MIN_LFC} in BOTH arms,\n"
              f"opposite directions, $\\geq${MIN_GRP} groups each  (n = {len(cand)}; top {len(top)} shown)",
              fontsize=17, fontweight="bold", loc="left", pad=14)
axA.tick_params(axis="x", labelsize=13)
axA.text(1.0, -0.075, "grey numbers = groups supporting the call (HIBI/GLP-1R)",
         transform=axA.transAxes, ha="right", va="top", fontsize=11.5, style="italic", color="#666")
for s in axA.spines.values():
    s.set_linewidth(1.2)

# Panel B: permutation null for the magnitude-weighted score
axB = fig.add_subplot(gs[0, 1])
axB.hist(null, bins=60, color="#BDC3C7", edgecolor="white", linewidth=0.5)
axB.axvline(score_obs, color=C_REV, lw=3.5, zorder=5)
axB.annotate(f"observed\n{score_obs:+.0f}", xy=(score_obs, axB.get_ylim()[1] * 0.72),
             xytext=(-95, 0), textcoords="offset points", fontsize=14, fontweight="bold",
             color=C_REV, ha="center",
             arrowprops=dict(arrowstyle="->", color=C_REV, lw=2.2))
axB.axvline(null.mean(), color="#555", ls="--", lw=2)
axB.text(null.mean(), axB.get_ylim()[1] * 0.96, " chance\n expectation", fontsize=12.5,
         color="#555", ha="left", va="top")
axB.set_xlabel("magnitude-weighted opposition score,  $\\Sigma\\,(-x\\cdot y)$",
               fontsize=15, fontweight="bold")
axB.set_ylabel("permutations", fontsize=15, fontweight="bold")
axB.set_title("B   Global test, weighted by effect size\n(20,000 permutations)",
              fontsize=17, fontweight="bold", loc="left", pad=14)
axB.tick_params(labelsize=12)
axB.text(0.5, -0.155,
         f"Higher score = more net opposition. Observed sits $z$ = {z:+.2f} BELOW chance:\n"
         "weighting by magnitude makes the two arms look MORE concordant, not less.\n"
         f"Pearson $r$ = {pr:+.3f} ($p$ = {pp:.1e}) vs Spearman $\\rho$ = {sr:+.3f}; requiring\n"
         "|log2FC| $\\geq$ 1.0 in both arms leaves 34% opposite, i.e. 66% same direction.",
         transform=axB.transAxes, ha="center", va="top", fontsize=12,
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#FEF9E7",
                   edgecolor="#B7950B", linewidth=1.5))
for s in axB.spines.values():
    s.set_linewidth(1.2)

for fmt in ("png", "pdf"):
    fig.savefig(os.path.join(OUTFIG, f"Fig12E_Magnitude_Opposition.{fmt}"),
                dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("saved Fig12E_Magnitude_Opposition (png + pdf)")
