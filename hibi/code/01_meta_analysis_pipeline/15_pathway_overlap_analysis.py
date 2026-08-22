#!/usr/bin/env python3
"""Pathway-level comparison of HIBI (neonatal hypoxic-ischemic brain injury) against GLP-1R agonism across the CNS.

Where Fig12C showed only which terms were significant in both arms, this asks for
every shared pathway whether the two conditions move it in the same or opposite
direction. It is more robust than the gene-level direction test because averaging
over a pathway's member genes suppresses the per-gene noise that made gene-level
means unreliable (see Fig12F).
It computes: term-level overlap overall and per database, with a hypergeometric
test against the terms actually tested; per shared pathway a direction score per
arm (mean log2FC of the pathway's member genes within that arm's consensus set)
plus the fraction of member genes up-regulated; the correlation of those pathway
scores between arms and a permutation test on the magnitude-weighted opposition
score; gene membership overlap (Jaccard) per shared pathway; and pathways unique
to each arm, ranked by significance.
Outputs: data/meta_analysis/pathway_overlap_full.csv (every shared pathway),
data/meta_analysis/pathway_arm_specific.csv (arm-unique pathways) and
figures/Fig13_Pathway_Overlap.{png,pdf} (4-panel summary)
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
from scipy.stats import hypergeom, pearsonr, spearmanr, binomtest

REPO = BASE + "/mnt/GLP-1R-bioinformatics-project"
HIBI = os.path.join(REPO, "hibi")
GLP_DATA = os.path.join(REPO, "data", "meta_analysis")
HIBI_DATA = os.path.join(HIBI, "data")
OUTFIG = os.path.join(HIBI, "figures")
OUTDATA = os.path.join(HIBI_DATA, "meta_analysis")

C_HIBI, C_GLP, C_REV, C_CONC = "#C0392B", "#2471A3", "#7D3C98", "#7F8C8D"
SIG = 0.05
rng = np.random.default_rng(42)

CATMAP = {"GO_Biological_Process_2023": "GO-BP", "GO_Molecular_Function_2023": "GO-MF",
          "GO_Cellular_Component_2023": "GO-CC", "KEGG_2021_Human": "KEGG",
          "Reactome_2022": "Reactome"}
CAT_SHORT = {"BP": "GO-BP", "MF": "GO-MF", "CC": "GO-CC",
             "KEGG": "KEGG", "Reactome": "Reactome"}


def clean_term(t):
    t = re.sub(r"\s*\(GO:\d+\)", "", str(t))
    return re.sub(r"\s*R-HSA-\d+", "", t).strip()


# Load
def load_glp():
    e = pd.read_csv(os.path.join(GLP_DATA, "enrichment_LFC02.csv"))
    e["term"] = e["Term"].map(clean_term)
    e["term_key"] = e["term"].str.lower()
    e["db"] = e["Gene_set"].map(CATMAP)
    e["padj"] = e["Adjusted P-value"]
    e["genes"] = e["Genes"].fillna("").str.upper().str.split(";")
    return e[["term", "term_key", "db", "padj", "genes"]]


def load_hibi():
    frames = []
    for f, default in [("consensus_terms_FULL.csv", None),
                       ("consensus_kegg_terms.csv", "KEGG"),
                       ("consensus_reactome_terms.csv", "Reactome")]:
        p = os.path.join(HIBI_DATA, "goplot", f)
        if os.path.exists(p):
            frames.append(pd.read_csv(p))
    e = pd.concat(frames, ignore_index=True)
    e["term"] = e["Term"].map(clean_term)
    e["term_key"] = e["term"].str.lower()
    e["db"] = e["Category"].map(CAT_SHORT)
    e["padj"] = e["adj_pval"]
    e["genes"] = (e["Genes"].fillna("").str.upper()
                  .str.replace(", ", ";", regex=False).str.split(";"))
    return e[["term", "term_key", "db", "padj", "genes"]]


glp_all, hib_all = load_glp(), load_hibi()
# Keep the most significant record per term
glp_all = glp_all.sort_values("padj").drop_duplicates("term_key")
hib_all = hib_all.sort_values("padj").drop_duplicates("term_key")
glp = glp_all[glp_all.padj < SIG].copy()
hib = hib_all[hib_all.padj < SIG].copy()

print(f"terms tested   : HIBI {len(hib_all):5d}   GLP-1R {len(glp_all):5d}")
print(f"terms significant (padj<{SIG}): HIBI {len(hib):4d}   GLP-1R {len(glp):4d}")

# 1. Overlap.
# Caveat: the HIBI enrichment tables were written already filtered to padj<0.05,
# so the HIBI arm contributes only its significant terms to the universe, while
# the GLP-1R table retains all 8,384 tested terms and dominates it. The universe
# is therefore an under-estimate and the fold-enrichment below an upper bound;
# re-run prepare_goplot_data.py without the significance filter for an exact
# figure.
universe = len(set(hib_all.term_key) | set(glp_all.term_key))
shared_keys = set(hib.term_key) & set(glp.term_key)
exp = len(hib) * len(glp) / universe
p_over = hypergeom.sf(len(shared_keys) - 1, universe, len(hib), len(glp))
print(f"\nshared significant terms: {len(shared_keys)}  "
      f"(expected {exp:.1f} by chance, {len(shared_keys)/exp:.2f}x, p={p_over:.2e})")

per_db = []
for db in ["GO-BP", "GO-MF", "GO-CC", "KEGG", "Reactome"]:
    h = set(hib[hib.db == db].term_key)
    g = set(glp[glp.db == db].term_key)
    u = len(set(hib_all[hib_all.db == db].term_key) | set(glp_all[glp_all.db == db].term_key))
    sh = h & g
    e_ = (len(h) * len(g) / u) if u else np.nan
    p_ = hypergeom.sf(len(sh) - 1, u, len(h), len(g)) if u and len(h) and len(g) else np.nan
    per_db.append(dict(db=db, hibi=len(h), glp=len(g), shared=len(sh),
                       expected=e_, fold=len(sh) / e_ if e_ else np.nan, p=p_))
    print(f"  {db:9s} HIBI {len(h):4d}  GLP-1R {len(g):4d}  shared {len(sh):4d}  "
          f"exp {e_:6.1f}  {len(sh)/e_ if e_ else float('nan'):.2f}x  p={p_:.2e}")
per_db = pd.DataFrame(per_db)

# 2. Direction
cg = pd.read_csv(os.path.join(GLP_DATA, "consensus_LFC02.csv"))
ch = pd.read_csv(os.path.join(OUTDATA, "consensus_LFC02.csv"))
lfc_g = {r.gene_symbol.upper(): r.mean_log2FC for r in cg.itertuples()}
lfc_h = {r.gene_symbol.upper(): r.mean_log2FC for r in ch.itertuples()}

hmap = hib.set_index("term_key")
gmap = glp.set_index("term_key")
rows = []
for k in shared_keys:
    rh, rg = hmap.loc[k], gmap.loc[k]
    gh = {x for x in rh.genes if x}
    gg = {x for x in rg.genes if x}
    members = gh | gg
    vh = [lfc_h[x] for x in members if x in lfc_h]
    vg = [lfc_g[x] for x in members if x in lfc_g]
    if len(vh) < 3 or len(vg) < 3:
        continue
    jac = len(gh & gg) / len(gh | gg) if (gh | gg) else np.nan
    rows.append(dict(
        term=rh.term, db=rh.db,
        padj_hibi=rh.padj, padj_glp=rg.padj,
        n_genes_hibi=len(vh), n_genes_glp=len(vg),
        jaccard_genes=jac,
        score_hibi=float(np.mean(vh)), score_glp=float(np.mean(vg)),
        frac_up_hibi=float(np.mean([v > 0 for v in vh])),
        frac_up_glp=float(np.mean([v > 0 for v in vg])),
    ))
P = pd.DataFrame(rows)
P["opposite"] = np.sign(P.score_hibi) != np.sign(P.score_glp)
P["combined_sig"] = -np.log10(P.padj_hibi.clip(lower=1e-300)) - np.log10(P.padj_glp.clip(lower=1e-300))
P = P.sort_values("combined_sig", ascending=False)
P.to_csv(os.path.join(OUTDATA, "pathway_overlap_full.csv"), index=False)

nopp = int(P.opposite.sum())
pr, ppr = pearsonr(P.score_hibi, P.score_glp)
sr, psr = spearmanr(P.score_hibi, P.score_glp)
pb = binomtest(nopp, len(P), 0.5).pvalue
print(f"\npathway-level direction ({len(P)} shared pathways with >=3 genes per arm):")
print(f"  opposite direction: {nopp}/{len(P)} ({nopp/len(P)*100:.0f}%)  binomial p={pb:.4f}")
print(f"  Pearson r  = {pr:+.3f} (p={ppr:.3g})")
print(f"  Spearman rho = {sr:+.3f} (p={psr:.3g})")

obs = float(-(P.score_hibi.values * P.score_glp.values).sum())
null = np.array([float(-(P.score_hibi.values * rng.permutation(P.score_glp.values)).sum())
                 for _ in range(20000)])
z = (obs - null.mean()) / null.std()
p_perm = (np.sum(null >= obs) + 1) / (len(null) + 1)
print(f"  magnitude-weighted opposition: obs {obs:+.2f}, null {null.mean():+.2f}+/-{null.std():.2f}, "
      f"z={z:+.2f}, p={p_perm:.4f}")

# 3. Arm-specific pathways
spec = []
for nm, d, other in [("HIBI", hib, set(glp.term_key)), ("GLP-1R", glp, set(hib.term_key))]:
    u = d[~d.term_key.isin(other)].nsmallest(400, "padj")
    for r in u.itertuples():
        spec.append(dict(arm=nm, term=r.term, db=r.db, padj=r.padj))
spec = pd.DataFrame(spec)
spec.to_csv(os.path.join(OUTDATA, "pathway_arm_specific.csv"), index=False)
print(f"\narm-specific pathways written: HIBI {int((spec.arm=='HIBI').sum())}, "
      f"GLP-1R {int((spec.arm=='GLP-1R').sum())} (top 400 each)")

# Figure
fig = plt.figure(figsize=(21, 17))
gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.26,
                      left=0.075, right=0.96, top=0.855, bottom=0.06)
fig.suptitle("Pathway-level comparison of the two conditions\n"
             "HIBI (neonatal hypoxic-ischemic injury) vs GLP-1R agonism",
             fontsize=24, fontweight="bold", y=0.955)

# Panel A: overlap per database
axA = fig.add_subplot(gs[0, 0])
yy = np.arange(len(per_db)); h = 0.26
axA.barh(yy + h, per_db.hibi, height=h, color=C_HIBI, alpha=0.9, label="HIBI significant")
axA.barh(yy, per_db.glp, height=h, color=C_GLP, alpha=0.9, label="GLP-1R significant")
axA.barh(yy - h, per_db.shared, height=h, color=C_REV, alpha=0.95, label="shared")
for i, r in per_db.iterrows():
    if r.shared > 0:
        axA.text(r.shared + max(per_db.hibi) * 0.012, i - h,
                 f"{r.fold:.1f}× exp", va="center", fontsize=11, color=C_REV, fontweight="bold")
axA.set_yticks(yy); axA.set_yticklabels(per_db.db, fontsize=14)
axA.set_xlabel("number of enriched terms (padj < 0.05)", fontsize=16, fontweight="bold")
axA.legend(fontsize=13, loc="lower right", framealpha=0.93)
axA.set_title(f"A   Term overlap by database\n{len(shared_keys)} shared overall "
              f"({len(shared_keys)/exp:.2f}× expected, p = {p_over:.1e})",
              fontsize=18, fontweight="bold", loc="left", pad=12)
axA.tick_params(axis="x", labelsize=12)

# Panel B: pathway direction scatter
axB = fig.add_subplot(gs[0, 1])
lim = float(np.nanpercentile(np.abs(np.r_[P.score_hibi, P.score_glp]), 98)) * 1.35
axB.axhline(0, color="#555", lw=1.2); axB.axvline(0, color="#555", lw=1.2)
axB.add_patch(plt.Rectangle((0, -lim), lim, lim, facecolor=C_REV, alpha=0.05, zorder=0))
axB.add_patch(plt.Rectangle((-lim, 0), lim, lim, facecolor=C_REV, alpha=0.05, zorder=0))
opp = P.opposite
axB.scatter(P.score_hibi[~opp].clip(-lim, lim), P.score_glp[~opp].clip(-lim, lim),
            s=70, c=C_CONC, alpha=0.65, edgecolors="white", linewidths=0.8,
            label=f"same direction (n={len(P)-nopp})")
axB.scatter(P.score_hibi[opp].clip(-lim, lim), P.score_glp[opp].clip(-lim, lim),
            s=80, c=C_REV, alpha=0.75, edgecolors="white", linewidths=0.8,
            label=f"opposite (n={nopp})")
placed = []
PADL = lim * 0.30          # keep labels clear of the panel edges
for r in P.head(45).itertuples():
    x, y = np.clip(r.score_hibi, -lim, lim), np.clip(r.score_glp, -lim, lim)
    if not (-lim + PADL * 0.5 < x < lim - PADL and -lim + PADL * 0.3 < y < lim - PADL * 0.3):
        continue
    fx, fy = (x + lim) / (2 * lim), (y + lim) / (2 * lim)
    if any((fx - a) ** 2 + (fy - b) ** 2 < 0.062 ** 2 for a, b in placed):
        continue
    placed.append((fx, fy))
    axB.annotate(textwrap.shorten(r.term, 24, placeholder="…"), (x, y), fontsize=10.5,
                 fontweight="bold", xytext=(6, 4), textcoords="offset points",
                 color=C_REV if r.opposite else "#2C3E50", clip_on=True,
                 bbox=dict(boxstyle="round,pad=0.1", facecolor="white", alpha=0.72, edgecolor="none"))
    if len(placed) >= 12:
        break
axB.set_xlim(-lim, lim); axB.set_ylim(-lim, lim); axB.set_aspect("equal")
axB.set_xlabel("HIBI pathway score  (mean log2FC of members)", fontsize=15,
               color=C_HIBI, fontweight="bold")
axB.set_ylabel("GLP-1R pathway score", fontsize=15, color=C_GLP, fontweight="bold")
axB.legend(fontsize=12, loc="upper left", framealpha=0.93)
axB.set_title("B   Do shared pathways move the same way?", fontsize=18,
              fontweight="bold", loc="left", pad=12)
axB.text(0.975, 0.03,
         f"opposite: {nopp}/{len(P)} ({nopp/len(P)*100:.0f}%), $p$ = {pb:.2f}\n"
         f"Pearson $r$ = {pr:+.2f} ($p$ = {ppr:.2g})\n"
         f"permutation $z$ = {z:+.2f}",
         transform=axB.transAxes, ha="right", va="bottom", fontsize=12,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#FEF9E7",
                   edgecolor="#B7950B", linewidth=1.4))
axB.tick_params(labelsize=12)

# Panel C: top shared pathways, significance in both arms
axC = fig.add_subplot(gs[1, 0])
top = P.head(15).iloc[::-1]
yy = np.arange(len(top)); h = 0.4
axC.barh(yy + h / 2, -np.log10(top.padj_hibi.clip(lower=1e-300)), height=h,
         color=C_HIBI, alpha=0.9, label="HIBI")
axC.barh(yy - h / 2, -np.log10(top.padj_glp.clip(lower=1e-300)), height=h,
         color=C_GLP, alpha=0.9, label="GLP-1R")
axC.set_yticks(yy)
axC.set_yticklabels(["\n".join(textwrap.wrap(t, 34)) for t in top.term], fontsize=11.5)
axC.axvline(-np.log10(0.05), color="#333", ls="--", lw=1.5)
axC.set_xlabel("$-\\log_{10}$ adjusted $p$", fontsize=16, fontweight="bold")
axC.legend(fontsize=13, loc="lower right", framealpha=0.93)
axC.set_title("C   Most significant shared pathways", fontsize=18,
              fontweight="bold", loc="left", pad=12)
axC.tick_params(axis="x", labelsize=12)

# Panel D: gene-membership agreement
axD = fig.add_subplot(gs[1, 1])
axD.hist(P.jaccard_genes.dropna(), bins=28, color="#5D6D7E", alpha=0.85, edgecolor="white")
axD.axvline(P.jaccard_genes.median(), color=C_REV, lw=3)
axD.text(P.jaccard_genes.median(), axD.get_ylim()[1] * 0.93,
         f"  median {P.jaccard_genes.median():.2f}", fontsize=13,
         color=C_REV, fontweight="bold", va="top")
axD.set_xlabel("Jaccard overlap of the pathway's driving genes\n(HIBI vs GLP-1R)",
               fontsize=15, fontweight="bold")
axD.set_ylabel("shared pathways", fontsize=15, fontweight="bold")
axD.set_title("D   Same pathway, same genes?", fontsize=18, fontweight="bold", loc="left", pad=12)
axD.tick_params(labelsize=12)
axD.text(0.5, -0.20,
         "Low Jaccard means the two conditions hit the same pathway through largely DIFFERENT member genes,\n"
         "which is itself a substantive result: shared enrichment does not imply a shared molecular route.",
         transform=axD.transAxes, ha="center", va="top", fontsize=12, style="italic", color="#555")

for fmt in ("png", "pdf"):
    fig.savefig(os.path.join(OUTFIG, f"Fig13_Pathway_Overlap.{fmt}"),
                dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("\nsaved Fig13_Pathway_Overlap (png + pdf)")
print(f"median gene-membership Jaccard across shared pathways: {P.jaccard_genes.median():.3f}")
print("\nTop 10 shared pathways by combined significance:")
for r in P.head(10).itertuples():
    d = "OPPOSITE" if r.opposite else "same"
    print(f"  {r.db:9s} {r.term[:44]:44s} HIBI {r.score_hibi:+.2f} / GLP {r.score_glp:+.2f}  {d}")
