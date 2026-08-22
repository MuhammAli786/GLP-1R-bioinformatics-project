#!/usr/bin/env python3
"""Standardize every GLP-1R-agonist / multi-agonist CNS DEG file to mouse gene symbols and build the master DEG table.

Per-group DEG CSVs (RNA-seq DESeq2 + limma arrays) -> Data/master_deg_<thr>.csv (long: group x gene) and Data/group_catalog.csv.
Thresholds are padj < 0.05 with |log2FC| >= 0.2 (LFC02), 0.5 (LFC05) and 1.0 (LFC1). Ensembl mouse IDs map to symbols via mygene; symbol columns are title-cased; arrays are already symbol-mapped by 01b. Each treatment-vs-control comparison within a region is its own group, never merged across regions or accessions.
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, glob, json, csv
import pandas as pd
import numpy as np

STAGE   = BASE + "/stage"
ARRAYS  = os.path.join(STAGE, "arrays")
OUTDATA = BASE + "/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/Data"
os.makedirs(OUTDATA, exist_ok=True)

# (g-index, accession, region, treatment, symbol source). Source is 'index' (symbol is the
# row index), 'gene_name' (use the gene_name column) or 'ensembl' (map gene_id ENSMUSG -> symbol).
RNASEQ = [
 # Accession 1 is excluded.
 (6,3,"accumbens","GLP_1","gene_name"), (7,3,"accumbens","GLP_1_MK_801","gene_name"),
 (8,3,"brainstem","GLP_1","gene_name"), (9,3,"brainstem","GLP_1_MK_801","gene_name"),
 (10,3,"hypothalamus","GLP_1","gene_name"), (11,3,"hypothalamus","GLP_1_MK_801","gene_name"),
 (12,3,"hypothalamus","conjugate_2","gene_name"),
 (13,4,"hypothalamus","Lira","gene_name"), (14,4,"hypothalamus","PYY_Lira","gene_name"),
 (15,5,"arc","IP118_PY115","ensembl"), (16,5,"nts","IP118_PY115","ensembl"),
 (17,5,"pvn","IP118_PY115","ensembl"),
 (18,6,"dvc","AC3174","ensembl"), (19,6,"dvc","AC710222","ensembl"), (20,6,"dvc","COMBO","ensembl"),
 (21,6,"mbh","AC3174","ensembl"), (22,6,"mbh","AC710222","ensembl"), (23,6,"mbh","COMBO","ensembl"),
 (24,7,"lumbar_spinal_cord","LIRA_EAE","gene_name"),
 (25,8,"frontalcortex","YA_Ex","gene_name"), (26,8,"hippocampus","YA_Ex","gene_name"),
 (27,8,"hypothalamus","YA_Ex","gene_name"), (28,8,"frontalcortex","LT_A_Ex","gene_name"),
 (29,8,"hippocampus","LT_A_Ex","gene_name"), (30,8,"hypothalamus","LT_A_Ex","gene_name"),
 (31,8,"frontalcortex","ST_A_Ex","gene_name"), (32,8,"frontalcortex","ST_A_Ex_KD","gene_name"),
 (33,8,"hippocampus","ST_A_Ex","gene_name"), (34,8,"hippocampus","ST_A_Ex_KD","gene_name"),
]

THRESHOLDS = {"LFC02": 0.2, "LFC05": 0.5, "LFC1": 1.0}


def mouse_sym(s):
    """Normalize a gene symbol to mouse title-case, or None if it is missing."""
    s = str(s).strip()
    if not s or s.lower() in ("nan", "na", "none"):
        return None
    return s[0].upper() + s[1:].lower() if len(s) > 1 else s.upper()


def build_ensembl_map(ens_ids):
    """Query mygene in batches of 1000 for mouse ENSMUSG -> symbol."""
    import mygene
    mg = mygene.MyGeneInfo()
    ens_ids = sorted(ens_ids)
    m = {}
    B = 1000
    for i in range(0, len(ens_ids), B):
        batch = ens_ids[i:i+B]
        res = mg.querymany(batch, scopes="ensembl.gene", fields="symbol",
                           species="mouse", returnall=True, verbose=False)
        for hit in res["out"]:
            if "symbol" in hit and "query" in hit:
                m[hit["query"]] = hit["symbol"]
    return m


def load_group_df(gidx, acc, source):
    """Load one group's DEG CSV as raw_sym/lfc/padj, or None if unusable."""
    path = os.path.join(STAGE, f"g{gidx:02d}_acc{acc}.csv")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    if source == "index":
        df = pd.read_csv(path, index_col=0)
        df = df.rename(columns={"log2FoldChange": "lfc"})
        df["raw_sym"] = df.index.astype(str)
    else:
        df = pd.read_csv(path)
        df = df.rename(columns={"log2FoldChange": "lfc"})
        if source == "gene_name":
            df["raw_sym"] = df["gene_name"].astype(str)
        elif source == "ensembl":
            df["raw_sym"] = df["gene_id"].astype(str)
    if "padj" not in df.columns or "lfc" not in df.columns:
        return None
    return df[["raw_sym", "lfc", "padj"]].copy()


def main():
    """Build and write the master DEG tables and group catalog."""
    ens_ids = set()
    for gidx, acc, region, treat, source in RNASEQ:
        if source != "ensembl":
            continue
        df = load_group_df(gidx, acc, source)
        if df is not None:
            ens_ids.update(x for x in df["raw_sym"] if str(x).startswith("ENSMUSG"))
    print(f"Ensembl IDs to map: {len(ens_ids)}")
    cache = BASE + "/work/ens_map.json"
    if os.path.exists(cache):
        ens_map = json.load(open(cache))
        print(f"Ensembl map loaded from cache: {len(ens_map)}")
    else:
        ens_map = build_ensembl_map(ens_ids) if ens_ids else {}
        print(f"Ensembl mapped: {len(ens_map)}")

    group_rows = {t: [] for t in THRESHOLDS}
    catalog = []
    skipped = []

    def add_group(acc, region, treat, df, gse, kind):
        gname = f"Acc{acc}_{region}_{treat}"
        syms = []
        for raw in df["raw_sym"]:
            if kind == "ensembl":
                mapped = ens_map.get(str(raw))
                syms.append(mouse_sym(mapped) if mapped else None)
            else:
                syms.append(mouse_sym(raw))
        df = df.assign(symbol=syms).dropna(subset=["symbol"])
        df = df[np.isfinite(df["lfc"]) & np.isfinite(df["padj"])]
        # dedupe per symbol: keep most significant (min padj, then max |lfc|)
        df["abslfc"] = df["lfc"].abs()
        df = df.sort_values(["padj", "abslfc"], ascending=[True, False])
        df = df.drop_duplicates(subset="symbol", keep="first")
        n_any = 0
        for thr, cut in THRESHOLDS.items():
            sig = df[(df["padj"] < 0.05) & (df["abslfc"] >= cut)]
            for _, r in sig.iterrows():
                group_rows[thr].append({
                    "group": gname, "accession": acc, "region": region,
                    "treatment": treat, "gse": gse, "kind": kind,
                    "symbol": r["symbol"], "symbol_key": r["symbol"].upper(),
                    "log2FC": r["lfc"], "padj": r["padj"],
                    "direction": "UP" if r["lfc"] > 0 else "DOWN",
                })
            if thr == "LFC02":
                n_any = len(sig)
        catalog.append({"group": gname, "accession": acc, "region": region,
                        "treatment": treat, "gse": gse, "kind": kind,
                        "n_sig_LFC02": n_any})

    gse_map = {1:"GSE314029",3:"GSE244607",4:"GSE190218",5:"GSE160802",
               6:"GSE135862",7:"GSE186318",8:"GSE280382"}
    for gidx, acc, region, treat, source in RNASEQ:
        df = load_group_df(gidx, acc, source)
        if df is None:
            skipped.append(f"Acc{acc}_{region}_{treat} (g{gidx:02d}) [LOCKED/empty]")
            continue
        add_group(acc, region, treat, df, gse_map.get(acc, ""), source)

    amani = json.load(open(os.path.join(ARRAYS, "array_manifest.json")))
    for a in amani:
        df = pd.read_csv(a["file"]).rename(columns={"log2FoldChange": "lfc"})
        df["raw_sym"] = df["gene_name"].astype(str)
        add_group(a["acc"], a["region"], a["treatment"], df[["raw_sym","lfc","padj"]],
                  a["gse"], "array")

    for thr in THRESHOLDS:
        m = pd.DataFrame(group_rows[thr])
        m.to_csv(os.path.join(OUTDATA, f"master_deg_{thr}.csv"), index=False)
        ng = m["group"].nunique() if len(m) else 0
        print(f"{thr}: {len(m)} rows, {ng} groups, {m['symbol_key'].nunique() if len(m) else 0} unique genes")
    pd.DataFrame(catalog).to_csv(os.path.join(OUTDATA, "group_catalog.csv"), index=False)
    print(f"\nGroups built: {len(catalog)}")
    if skipped:
        print("SKIPPED (locked/empty inputs):")
        for s in skipped: print("  -", s)
    json.dump({"skipped": skipped}, open(os.path.join(OUTDATA, "_skipped.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
