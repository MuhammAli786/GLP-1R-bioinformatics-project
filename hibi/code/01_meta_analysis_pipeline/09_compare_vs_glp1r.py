#!/usr/bin/env python3
"""
Compare HIBI consensus DEGs (neonatal hypoxic-ischemic brain injury, this
project) vs GLP-1R-agonism consensus DEGs (GLP-1R-bioinformatics-project).
Both at LFC02 threshold (padj<0.05, |log2FC|>=0.2, gene sig in >=2 groups).
"""
import pandas as pd

HIBI = "/sessions/lucid-pensive-ride/mnt/outputs/hibi_data/meta_analysis/consensus_LFC02.csv"
GLP1R = "/sessions/lucid-pensive-ride/mnt/GLP-1R-bioinformatics-project/data/meta_analysis/consensus_LFC02.csv"
OUT = "/sessions/lucid-pensive-ride/mnt/outputs/hibi_vs_glp1r_consensus_comparison.csv"
SUMMARY = "/sessions/lucid-pensive-ride/mnt/outputs/comparison_summary.txt"

hibi = pd.read_csv(HIBI)
glp = pd.read_csv(GLP1R)

hibi["key"] = hibi["gene_symbol"].str.upper()
glp["key"] = glp["gene_symbol"].str.upper()

hibi_set = set(hibi["key"])
glp_set = set(glp["key"])
shared = hibi_set & glp_set
hibi_only = hibi_set - glp_set
glp_only = glp_set - hibi_set

merged = hibi.merge(glp, on="key", suffixes=("_hibi", "_glp1r"), how="inner")
merged["direction_concordant"] = merged["predominant_direction_hibi"] == merged["predominant_direction_glp1r"]
merged_out = merged[[
    "gene_symbol_hibi", "n_groups_hibi", "predominant_direction_hibi", "mean_log2FC_hibi",
    "n_groups_glp1r", "predominant_direction_glp1r", "mean_log2FC_glp1r", "direction_concordant"
]].rename(columns={"gene_symbol_hibi": "gene_symbol"})
merged_out = merged_out.sort_values(["direction_concordant", "n_groups_hibi", "n_groups_glp1r"],
                                     ascending=[False, False, False])
merged_out.to_csv(OUT, index=False)

n_conc = merged_out["direction_concordant"].sum()
n_disc = (~merged_out["direction_concordant"]).sum()

lines = []
lines.append("HIBI (hypoxic-ischemic brain injury) vs GLP-1R (GLP-1R agonism, CNS) consensus DEG comparison")
lines.append("Threshold: padj<0.05, |log2FC|>=0.2, gene significant in >=2 groups within each project")
lines.append("")
lines.append(f"HIBI consensus genes:   {len(hibi_set)}  (from {12} HIBI groups: 6 GEO datasets)")
lines.append(f"GLP-1R consensus genes: {len(glp_set)}  (from GLP-1R-bioinformatics-project, 9 accessions)")
lines.append(f"Shared consensus genes: {len(shared)}")
lines.append(f"  - same predominant direction (concordant):    {n_conc}")
lines.append(f"  - opposite predominant direction (discordant): {n_disc}")
lines.append(f"HIBI-only consensus genes:   {len(hibi_only)}")
lines.append(f"GLP-1R-only consensus genes: {len(glp_only)}")
lines.append("")
lines.append("Top 20 shared, concordant-direction genes (by combined recurrence):")
top_conc = merged_out[merged_out["direction_concordant"]].copy()
top_conc["combined_n"] = top_conc["n_groups_hibi"] + top_conc["n_groups_glp1r"]
top_conc = top_conc.sort_values("combined_n", ascending=False).head(20)
for _, r in top_conc.iterrows():
    lines.append(f"  {r['gene_symbol']}: HIBI n_groups={r['n_groups_hibi']} ({r['predominant_direction_hibi']}), "
                 f"GLP-1R n_groups={r['n_groups_glp1r']} ({r['predominant_direction_glp1r']})")
lines.append("")
lines.append("Discordant-direction shared genes (opposite regulation between injury and GLP-1R agonism):")
disc = merged_out[~merged_out["direction_concordant"]].sort_values(
    ["n_groups_hibi", "n_groups_glp1r"], ascending=False).head(20)
for _, r in disc.iterrows():
    lines.append(f"  {r['gene_symbol']}: HIBI={r['predominant_direction_hibi']} (n={r['n_groups_hibi']}), "
                 f"GLP-1R={r['predominant_direction_glp1r']} (n={r['n_groups_glp1r']})")

summary = "\n".join(lines)
print(summary)
with open(SUMMARY, "w") as f:
    f.write(summary + "\n")
