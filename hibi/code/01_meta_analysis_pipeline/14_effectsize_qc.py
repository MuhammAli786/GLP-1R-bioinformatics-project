#!/usr/bin/env python3
"""Quality-control figure asking how large the compared effect sizes are and whether the cross-arm comparison survives removing near-zero genes.

Panel A: |mean log2FC| distribution per arm; the GLP-1R arm is compressed
(median 0.27 vs 0.87), with 34% of genes below 0.2.
Panel B: GLP-1R mean log2FC collapses toward zero when a gene's groups disagree
in direction, since averaging 31 heterogeneous groups (different agonists,
regions, disease contexts) cancels out.
Panel C: effect-size floor sweep. Requiring a minimum |mean log2FC| in both arms,
the direction split stays near 50% and the correlation stays weakly positive;
nominal significance at low thresholds is carried by near-zero genes and
disappears by 0.4.
Near-zero genes affect both sides and are more common among discordant genes
(19% vs 10% below 0.05), so removing them does not reveal a reversal signal.
Outputs: Fig12F_EffectSize_QC.{png,pdf}
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, binomtest

REPO = BASE + "/mnt/GLP-1R-bioinformatics-project"
HIBI = os.path.join(REPO, "hibi")
OUTFIG = os.path.join(HIBI, "figures")
C_HIBI, C_GLP, C_REV = "#C0392B", "#2471A3", "#7D3C98"

g = pd.read_csv(os.path.join(REPO, "data", "meta_analysis", "consensus_LFC02.csv"))
h = pd.read_csv(os.path.join(HIBI, "data", "meta_analysis", "consensus_LFC02.csv"))
for d in (g, h):
    d["key"] = d.gene_symbol.str.upper()
m = h.merge(g, on="key", suffixes=("_h", "_g"))
m = m[(m.mean_log2FC_h.abs() < 8) & (m.mean_log2FC_g.abs() < 8)].copy()
m["cons_g"] = m[["n_up_g", "n_down_g"]].max(axis=1) / (m.n_up_g + m.n_down_g)
m["amin"] = np.minimum(m.mean_log2FC_h.abs(), m.mean_log2FC_g.abs())
m["conc"] = m.predominant_direction_h == m.predominant_direction_g

fig, axes = plt.subplots(1, 3, figsize=(23, 7.5))
fig.suptitle("Effect-size quality control: how much movement is actually being compared?",
             fontsize=23, fontweight="bold", y=1.02)

# Panel A
ax = axes[0]
bins = np.linspace(0, 3, 55)
ax.hist(m.mean_log2FC_h.abs().clip(0, 3), bins=bins, alpha=0.62, color=C_HIBI,
        label=f"HIBI  (median {m.mean_log2FC_h.abs().median():.2f})")
ax.hist(m.mean_log2FC_g.abs().clip(0, 3), bins=bins, alpha=0.62, color=C_GLP,
        label=f"GLP-1R  (median {m.mean_log2FC_g.abs().median():.2f})")
ax.axvline(0.2, color="#333", ls="--", lw=2)
ax.text(0.21, ax.get_ylim()[1] * 0.93, " |LFC| = 0.2\n (the DEG threshold)", fontsize=12, va="top")
ax.set_xlabel("|mean log2FC|  across the shared genes", fontsize=15, fontweight="bold")
ax.set_ylabel("genes", fontsize=15, fontweight="bold")
ax.legend(fontsize=13)
ax.set_title("A   The GLP-1R arm is compressed\n"
             f"{(m.mean_log2FC_g.abs()<0.2).mean()*100:.0f}% of GLP-1R values are below 0.2, "
             f"vs {(m.mean_log2FC_h.abs()<0.2).mean()*100:.0f}% for HIBI",
             fontsize=15.5, fontweight="bold", loc="left", pad=12)
ax.tick_params(labelsize=12)

# Panel B
ax = axes[1]
groups = [("groups disagree\n(<70% agree)", m[m.cons_g < 0.7]),
          ("mostly agree\n(70-99%)", m[(m.cons_g >= 0.7) & (m.cons_g < 0.99)]),
          ("unanimous\n(100%)", m[m.cons_g >= 0.99])]
data = [s.mean_log2FC_g.abs().values for _, s in groups]
bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False)
for patch, col in zip(bp["boxes"], ["#AEB6BF", "#7FB3D5", C_GLP]):
    patch.set_facecolor(col); patch.set_alpha(0.8)
for med in bp["medians"]:
    med.set_color("black"); med.set_linewidth(2)
ax.set_xticklabels([f"{lab}\nn={len(s)}" for lab, s in groups], fontsize=12.5)
ax.set_ylabel("|mean log2FC|  (GLP-1R arm)", fontsize=15, fontweight="bold")
ax.set_title("B   Why: cancellation across the 31 GLP-1R groups\n"
             "genes whose groups disagree average out to ~0",
             fontsize=15.5, fontweight="bold", loc="left", pad=12)
ax.tick_params(labelsize=12)
for i, (_, s) in enumerate(groups, 1):
    ax.text(i, s.mean_log2FC_g.abs().median(), f"  {s.mean_log2FC_g.abs().median():.2f}",
            fontsize=12, fontweight="bold", va="bottom")

# Panel C
ax = axes[2]
ts = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]
pct, rs, ns, ps = [], [], [], []
for t in ts:
    s = m[m.amin >= t]
    o = int((np.sign(s.mean_log2FC_h) != np.sign(s.mean_log2FC_g)).sum())
    pct.append(o / len(s) * 100); ns.append(len(s))
    rs.append(pearsonr(s.mean_log2FC_h, s.mean_log2FC_g)[0])
    ps.append(binomtest(o, len(s), 0.5).pvalue)
ax.plot(ts, pct, "-o", color=C_REV, lw=3, ms=10, label="% opposite direction")
ax.axhline(50, color="#333", ls="--", lw=2)
ax.text(1.0, 50.6, "chance (50%)", fontsize=12, ha="right", color="#333")
for t, p_, n_, pv in zip(ts, pct, ns, ps):
    ax.annotate(f"n={n_}", (t, p_), fontsize=10, xytext=(0, -18),
                textcoords="offset points", ha="center", color="#555")
ax.set_ylim(25, 62)
ax.set_xlabel("minimum |mean log2FC| required in BOTH arms", fontsize=15, fontweight="bold")
ax.set_ylabel("% of shared genes moving oppositely", fontsize=15, fontweight="bold", color=C_REV)
ax2 = ax.twinx()
ax2.plot(ts, rs, "-s", color="#E67E22", lw=2.5, ms=8, label="Pearson r")
ax2.axhline(0, color="#E67E22", ls=":", lw=1.5)
ax2.set_ylabel("Pearson $r$", fontsize=15, fontweight="bold", color="#E67E22")
ax2.set_ylim(-0.25, 0.25)
ax2.tick_params(labelsize=12, colors="#E67E22")
ax.tick_params(labelsize=12)
ax.set_title("C   Removing near-zero genes does not\nreveal a reversal",
             fontsize=15.5, fontweight="bold", loc="left", pad=12)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=12, loc="lower left")

fig.text(0.5, -0.06,
         "Near-zero genes are in fact MORE common among discordant genes (19% vs 10% below |LFC| 0.05), so they were inflating the opposite-direction count, not the concordant one.\n"
         "Direction calls are also unreliable there: vote-direction and sign(mean) disagree for 32% of genes below |LFC| 0.2, versus 7% above 0.5.",
         ha="center", va="top", fontsize=13, style="italic", color="#444")

fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(os.path.join(OUTFIG, f"Fig12F_EffectSize_QC.{fmt}"),
                dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("saved Fig12F_EffectSize_QC (png + pdf)")
for t, p_, r_, n_, pv in zip(ts, pct, rs, ns, ps):
    print(f"  floor {t:4.2f}: n={n_:5d}  opposite {p_:5.1f}%  binom p={pv:.4f}  Pearson r={r_:+.3f}")
