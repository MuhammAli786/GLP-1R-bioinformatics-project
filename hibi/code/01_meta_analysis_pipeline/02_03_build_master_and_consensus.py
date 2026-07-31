#!/usr/bin/env python3
"""
Convert HIBI vs Vehicle dataset into the GLP-1R-bioinformatics-project layout
(meta_analysis stage: master_deg / group_catalog / consensus), then compare
HIBI consensus vs GLP-1R consensus.

Only the 6 HIBI GEO datasets with existing finished DEG tables are used
(per user decision): GSE144455, GSE144456, GSE23317, GSE23319, GSE236133,
GSE238220. GSE97299/294909/312452/36215 excluded (no saved DEG table, raw
data not present, R/Bioconductor unavailable in this environment).

Threshold: padj < 0.05 & |log2FC| >= 0.2 only (LFC02), per user decision.
Consensus: gene significant in >=2 groups, mirrors GLP-1R's 03_consensus.py.
Dedup rule (per group, per gene symbol): keep row with min padj, then max |lfc|
  -- mirrors GLP-1R's 02_build_master.py.
"""
import os, glob, gzip, csv
import pandas as pd
import numpy as np

HIBI_ROOT = "/sessions/lucid-pensive-ride/mnt/HIBI vs Vehicle/bulk_rna_hibi"
OUT = "/sessions/lucid-pensive-ride/mnt/outputs/hibi_data/meta_analysis"
os.makedirs(OUT, exist_ok=True)

PADJ_CUT = 0.05
LFC_CUT = 0.2
MIN_GROUPS = 2

def mouse_sym(s):
    s = str(s).strip()
    if not s or s.lower() in ("nan", "na", "none", ""):
        return None
    return s[0].upper() + s[1:].lower() if len(s) > 1 else s.upper()

rows = []          # collected sig rows across all groups
catalog = []        # group metadata

def add_group(group, accession, region, treatment, gse, kind, df):
    """df must have columns: raw_symbol, lfc, padj"""
    df = df.copy()
    df["symbol"] = df["raw_symbol"].apply(mouse_sym)
    df = df.dropna(subset=["symbol"])
    df = df[np.isfinite(df["lfc"]) & np.isfinite(df["padj"])]
    df["abslfc"] = df["lfc"].abs()
    df = df.sort_values(["padj", "abslfc"], ascending=[True, False])
    df = df.drop_duplicates(subset="symbol", keep="first")
    sig = df[(df["padj"] < PADJ_CUT) & (df["abslfc"] >= LFC_CUT)]
    for _, r in sig.iterrows():
        rows.append({
            "group": group, "accession": accession, "region": region,
            "treatment": treatment, "gse": gse, "kind": kind,
            "symbol": r["symbol"], "symbol_key": r["symbol"].upper(),
            "log2FC": r["lfc"], "padj": r["padj"],
            "direction": "UP" if r["lfc"] > 0 else "DOWN",
        })
    catalog.append({"group": group, "accession": accession, "region": region,
                     "treatment": treatment, "gse": gse, "kind": kind,
                     "n_sig_LFC02": len(sig)})
    print(f"{group}: {len(sig)} sig genes (of {len(df)} tested)")

# ---------------------------------------------------------------
# 1. GSE144455 -- forebrain, HI vs naive, 3h  (array, MgSO4 study 2-channel)
# ---------------------------------------------------------------
f = os.path.join(HIBI_ROOT, "results/GSE144455_HI_vs_naive_3h_all_annotated.csv")
d = pd.read_csv(f)
d = d[d["GENE_SYMBOL"].notna() & (d["GENE_SYMBOL"].astype(str).str.strip() != "")]
d = d.rename(columns={"GENE_SYMBOL": "raw_symbol", "logFC": "lfc", "adj.P.Val": "padj"})
add_group("HIBI144455_forebrain_HI_vs_naive_3h", "144455", "forebrain",
          "HI_vs_naive_3h", "GSE144455", "array", d[["raw_symbol", "lfc", "padj"]])

# ---------------------------------------------------------------
# 2. GSE144456 -- forebrain, P5 HI vs control, 3h (array, 2-channel)
# ---------------------------------------------------------------
f = os.path.join(HIBI_ROOT, "results/GSE144456_P5_3h_all_annotated.csv")
d = pd.read_csv(f)
d = d[d["GENE_SYMBOL"].notna() & (d["GENE_SYMBOL"].astype(str).str.strip() != "")]
d = d.rename(columns={"GENE_SYMBOL": "raw_symbol", "logFC": "lfc", "adj.P.Val": "padj"})
add_group("HIBI144456_forebrain_P5_HI_vs_control_3h", "144456", "forebrain",
          "P5_HI_vs_control_3h", "GSE144456", "array", d[["raw_symbol", "lfc", "padj"]])

# ---------------------------------------------------------------
# 3+4. GSE23317 (cortex) / GSE23319 (striatum) -- HI vs Sham, 3h
#      Illumina GPL6885 -- probe IDs need annotation mapping
# ---------------------------------------------------------------
annot = {}
with gzip.open("/tmp/gpl6885.annot.gz", "rt") as fh:
    started = False
    reader = csv.reader(fh, delimiter="\t")
    for row in reader:
        if not row:
            continue
        if row[0] == "ID":
            started = True
            continue
        if started and len(row) > 2 and row[0].startswith("ILMN"):
            annot[row[0]] = row[2]  # Gene symbol column

for gse_id, region in [("GSE23317", "cortex"), ("GSE23319", "striatum")]:
    f = os.path.join(HIBI_ROOT, f"results/{gse_id}_3h_all_genes.csv")
    d = pd.read_csv(f)
    d.columns = ["probe_id"] + list(d.columns[1:])
    d["raw_symbol"] = d["probe_id"].map(annot)
    d = d.rename(columns={"logFC": "lfc", "adj.P.Val": "padj"})
    acc = gse_id.replace("GSE", "")
    add_group(f"HIBI{acc}_{region}_HI_vs_Sham_3h", acc, region, "HI_vs_Sham_3h",
              gse_id, "array", d[["raw_symbol", "lfc", "padj"]])

# ---------------------------------------------------------------
# 5. GSE236133 -- hippocampus, ipsilateral vs contralateral, WT/NEIL1KO/NEIL2KO x 3h/6h
# ---------------------------------------------------------------
for genotype in ["WT", "NEIL1KO", "NEIL2KO"]:
    for time in ["3h", "6h"]:
        f = os.path.join(HIBI_ROOT, f"results/GSE236133/DE/{genotype}_{time}_ipsilateral_vs_contralateral_all_genes.csv")
        d = pd.read_csv(f)
        d = d[d["SYMBOL"].notna() & (d["SYMBOL"].astype(str).str.strip() != "")]
        d = d.rename(columns={"SYMBOL": "raw_symbol", "logFC": "lfc", "adj.P.Val": "padj"})
        add_group(f"HIBI236133_hippocampus_{genotype}_ipsi_vs_contra_{time}", "236133",
                  "hippocampus", f"{genotype}_ipsi_vs_contra_{time}", "GSE236133", "rnaseq",
                  d[["raw_symbol", "lfc", "padj"]])

# ---------------------------------------------------------------
# 6. GSE238220 -- hippocampus (sorted microglia), HI vs control, 1d/3d
# ---------------------------------------------------------------
for time in ["1d", "3d"]:
    f = os.path.join(HIBI_ROOT, f"results/GSE238220/DE/Microglia_{time}_HI_vs_Control_all_genes.csv")
    d = pd.read_csv(f)
    d = d[d["SYMBOL"].notna() & (d["SYMBOL"].astype(str).str.strip() != "")]
    d = d.rename(columns={"SYMBOL": "raw_symbol", "log2FoldChange": "lfc", "padj": "padj"})
    add_group(f"HIBI238220_hippocampus_microglia_HI_vs_control_{time}", "238220",
              "hippocampus_microglia", f"HI_vs_control_{time}", "GSE238220", "rnaseq",
              d[["raw_symbol", "lfc", "padj"]])

# ---------------------------------------------------------------
# Write master_deg + group_catalog
# ---------------------------------------------------------------
master = pd.DataFrame(rows, columns=["group", "accession", "region", "treatment", "gse",
                                      "kind", "symbol", "symbol_key", "log2FC", "padj", "direction"])
master.to_csv(os.path.join(OUT, "master_deg_LFC02.csv"), index=False)

cat = pd.DataFrame(catalog)
cat.to_csv(os.path.join(OUT, "group_catalog.csv"), index=False)

print(f"\nTotal groups: {len(cat)}; total sig rows: {len(master)}")

# ---------------------------------------------------------------
# Consensus (mirrors 03_consensus.py)
# ---------------------------------------------------------------
g = master.groupby("symbol_key")
cons_rows = []
for key, sub in g:
    groups = sub["group"].unique()
    n = len(groups)
    if n < MIN_GROUPS:
        continue
    disp = sub["symbol"].iloc[0]
    up = (sub["direction"] == "UP").sum()
    down = (sub["direction"] == "DOWN").sum()
    cons_rows.append({
        "gene_symbol": disp, "n_groups": n,
        "predominant_direction": "UP" if up >= down else "DOWN",
        "n_up": int(up), "n_down": int(down),
        "mean_log2FC": sub["log2FC"].mean(),
        "groups": ";".join(sorted(groups)),
    })
cons = pd.DataFrame(cons_rows).sort_values(["n_groups", "mean_log2FC"], ascending=[False, False])
cons.to_csv(os.path.join(OUT, "consensus_LFC02.csv"), index=False)

print(f"HIBI consensus (>= {MIN_GROUPS} groups): {len(cons)} genes; "
      f"max n_groups={cons['n_groups'].max() if len(cons) else 0}")
if len(cons):
    print("Top 15:", ", ".join(f"{r.gene_symbol}({r.n_groups},{r.predominant_direction})"
                                for r in cons.head(15).itertuples()))
