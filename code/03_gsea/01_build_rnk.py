#!/usr/bin/env python3
"""
01_build_rnk.py
-------------------------------------------------------------
Build a FULL ranked gene list (.rnk) per group for GSEA prerank.
Unlike the consensus pipeline (which keeps only significant DEGs), GSEA
needs EVERY gene ranked. Genes are standardized to mouse symbols and
then UPPER-cased so they match the human MSigDB/KEGG gene-set libraries
(ortholog-by-name, the standard mouse->human approximation).

Groups with no DEGs are skipped (Acc7 spinal cord, Acc8 hippocampus LT_A_Ex).
Output: GSEA/data/rnk/<group>.rnk   (gene<TAB>log2FC, sorted desc)
        GSEA/data/sig_genes.json    (per-group significant DEGs for ORA/Venn)
"""
import os, json
import numpy as np
import pandas as pd

STAGE = "/sessions/amazing-zen-bardeen/stage"
ARR = os.path.join(STAGE, "arrays")
ENSMAP = "/sessions/amazing-zen-bardeen/work/ens_map.json"
OUT = "/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/GSEA/data"
RNK = os.path.join(OUT, "rnk")
os.makedirs(RNK, exist_ok=True)

# gidx, acc, region, treatment, source ; exclude the 2 no-DEG groups (24, 29)
RNASEQ = [
 (6,3,"accumbens","GLP_1","gene_name"),(7,3,"accumbens","GLP_1_MK_801","gene_name"),
 (8,3,"brainstem","GLP_1","gene_name"),(9,3,"brainstem","GLP_1_MK_801","gene_name"),
 (10,3,"hypothalamus","GLP_1","gene_name"),(11,3,"hypothalamus","GLP_1_MK_801","gene_name"),
 (12,3,"hypothalamus","conjugate_2","gene_name"),
 (13,4,"hypothalamus","Lira","gene_name_rat"),(14,4,"hypothalamus","PYY_Lira","gene_name_rat"),
 (15,5,"arc","IP118_PY115","ensembl"),(16,5,"nts","IP118_PY115","ensembl"),(17,5,"pvn","IP118_PY115","ensembl"),
 (18,6,"dvc","AC3174","ensembl"),(20,6,"dvc","COMBO","ensembl"),
 (21,6,"mbh","AC3174","ensembl"),(23,6,"mbh","COMBO","ensembl"),
 # AC710222 (Acc6 dvc/mbh) EXCLUDED: it is the CCK1R agonist, not a GLP-1R agonist.
 (25,8,"frontalcortex","YA_Ex","gene_name"),(26,8,"hippocampus","YA_Ex","gene_name"),
 (27,8,"hypothalamus","YA_Ex","gene_name"),(28,8,"frontalcortex","LT_A_Ex","gene_name"),
 (30,8,"hypothalamus","LT_A_Ex","gene_name"),(31,8,"frontalcortex","ST_A_Ex","gene_name"),
 (32,8,"frontalcortex","ST_A_Ex_KD","gene_name"),(33,8,"hippocampus","ST_A_Ex","gene_name"),
 (34,8,"hippocampus","ST_A_Ex_KD","gene_name"),
]
ens_map = json.load(open(ENSMAP)) if os.path.exists(ENSMAP) else {}


def mouse_upper(s):
    s = str(s).strip()
    if not s or s.lower() in ("nan", "na", "none"):
        return None
    return s.upper()


def build_one(df, source):
    df = df.rename(columns={"log2FoldChange": "lfc"})
    if "lfc" not in df.columns:
        return None
    if source == "ensembl":
        df["sym"] = df["gene_id"].astype(str).map(lambda x: ens_map.get(x))
    else:  # gene_name / gene_name_rat / array all carry a symbol in gene_name
        df["sym"] = df["gene_name"].astype(str)
    df = df.dropna(subset=["sym"])
    df["sym"] = df["sym"].map(mouse_upper)
    df = df.dropna(subset=["sym"])
    df = df[np.isfinite(df["lfc"])]
    df["abslfc"] = df["lfc"].abs()
    df = df.sort_values("abslfc", ascending=False).drop_duplicates("sym", keep="first")
    return df[["sym", "lfc"]].sort_values("lfc", ascending=False)


def main():
    cat = pd.read_csv("/sessions/amazing-zen-bardeen/work/gsea_groups.csv")
    keep_groups = set(cat["group"])
    n = 0
    for gidx, acc, region, treat, source in RNASEQ:
        gname = f"Acc{acc}_{region}_{treat}"
        if gname not in keep_groups:
            continue
        p = os.path.join(STAGE, f"g{gidx:02d}_acc{acc}.csv")
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            print("missing", gname); continue
        rnk = build_one(pd.read_csv(p), source)
        if rnk is None or rnk.empty:
            print("empty", gname); continue
        rnk.to_csv(os.path.join(RNK, f"{gname}.rnk"), sep="\t", header=False, index=False)
        n += 1
    # arrays
    amani = json.load(open(os.path.join(ARR, "array_manifest.json")))
    for a in amani:
        gname = f"Acc{a['acc']}_{a['region']}_{a['treatment']}"
        if gname not in keep_groups:
            continue
        rnk = build_one(pd.read_csv(a["file"]), "array")
        if rnk is not None and not rnk.empty:
            rnk.to_csv(os.path.join(RNK, f"{gname}.rnk"), sep="\t", header=False, index=False)
            n += 1
    print(f"wrote {n} rnk files")

    # per-group significant DEGs (for ORA / Venn / UpSet), UPPER-cased
    m = pd.read_csv("/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/Data/master_deg_LFC02.csv")
    m = m[m["group"].isin(keep_groups)]
    sig = {g: sorted(set(s["symbol_key"])) for g, s in m.groupby("group")}
    json.dump(sig, open(os.path.join(OUT, "sig_genes.json"), "w"))
    print(f"sig_genes.json: {len(sig)} groups, "
          f"{sum(len(v) for v in sig.values())} total DEG entries")


if __name__ == "__main__":
    main()
