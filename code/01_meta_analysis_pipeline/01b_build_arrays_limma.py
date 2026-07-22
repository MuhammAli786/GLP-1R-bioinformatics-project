#!/usr/bin/env python3
"""
01b_build_arrays_limma.py
-------------------------------------------------------------
GEO2R / limma re-analysis of the three Illumina GPL6885 array
datasets (Acc10 GSE71850, Acc11 GSE41345, Acc12 GSE39586).

For each dataset we:
  1. parse the GEO series matrix -> expression matrix (probe x GSM)
  2. assign samples to groups (positional, from sample order)
  3. run limma (limma_array_engine.R): lmFit(~0+group), eBayes(0.01),
     topTable(adjust="fdr") for each GLP-1R-agonist-vs-control contrast
  4. map ILMN probes -> mouse Gene symbol via the GPL6885 AnnotGPL file
  5. write a DEG CSV per contrast in the same schema as the RNA-seq
     DESeq2 outputs (gene_name, log2FoldChange, padj, direction)

Outputs go to <STAGE>/arrays/ and an array manifest is written.
"""
import os, gzip, csv, io, json, subprocess, sys

WORK  = "/sessions/amazing-zen-bardeen/work"
STAGE = "/sessions/amazing-zen-bardeen/stage/arrays"
RSCRIPT = "/sessions/amazing-zen-bardeen/mamba/envs/bio/bin/Rscript"
ENGINE  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "limma_array_engine.R")
ANNOT   = os.path.join(WORK, "GPL6885.annot.gz")
os.makedirs(STAGE, exist_ok=True)

# ---- dataset configs: positional sample groups + GLP contrasts ----
# group lists are 0-based sample indices in series-matrix column order
DATASETS = {
  "GSE71850": {   # Acc10, hippocampus (blast cohort)
    "acc": 10, "region": "hippocampus",
    "groups": {
        "Control":     list(range(16, 21)),
        "Exendin4":    list(range(21, 26)),
        "Ex4_PreTBI":  list(range(26, 31)),
        "TBI":         list(range(31, 36)),
        "Ex4_PostTBI": list(range(36, 40)),
    },
    "contrasts": [   # (contrast string, output treatment label)
        ("Exendin4-Control",    "Exendin4"),
        ("Ex4_PreTBI-Control",  "Ex4_Pre_TBI"),
        ("Ex4_PostTBI-Control", "Ex4_Po_t_TBI"),
    ],
  },
  "GSE41345": {   # Acc11, hippocampus
    "acc": 11, "region": "hippocampus",
    "groups": {
        "Ex4":      list(range(0, 5)),
        "mTBI":     list(range(5, 9)),
        "Ex4_mTBI": list(range(9, 13)),
        "sham":     list(range(13, 18)),
    },
    "contrasts": [
        ("Ex4-sham",      "Ex4"),
        ("Ex4_mTBI-sham", "Ex4_mTBI"),
    ],
  },
  "GSE39586": {   # Acc12, hypothalamus -- user-supplied gsms "00000111XXX..."
    "acc": 12, "region": "hypothalamus",
    "groups": {
        "nORM": list(range(0, 5)),
        "glp":  list(range(5, 8)),
    },
    "contrasts": [
        ("glp-nORM", "GLP"),   # treatment (glp) vs control (nORM)
    ],
  },
}


def parse_series_matrix(path):
    """Return (gsm_ids list, dict probe-> [values str]) from a GEO series matrix."""
    gsms = None
    rows = {}
    with gzip.open(path, "rt", errors="replace") as f:
        in_tbl = False
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                gsms = next(csv.reader(io.StringIO(line.strip()), delimiter="\t"))[1:]
            elif line.startswith("!series_matrix_table_begin"):
                in_tbl = True
                header = next(f)
                hdr = next(csv.reader(io.StringIO(header.strip()), delimiter="\t"))
                col_gsms = hdr[1:]
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_tbl:
                parts = next(csv.reader(io.StringIO(line.rstrip("\n")), delimiter="\t"))
                if not parts or not parts[0]:
                    continue
                rows[parts[0].strip('"')] = parts[1:]
    if gsms is None:
        gsms = col_gsms
    return col_gsms, rows


def load_probe_symbol(annot_path):
    m = {}
    with gzip.open(annot_path, "rt", errors="replace") as f:
        header = None
        idx_id = idx_sym = None
        for line in f:
            if line.startswith(("#", "^", "!")):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                if "ID" in parts and "Gene symbol" in parts:
                    header = parts
                    idx_id = header.index("ID")
                    idx_sym = header.index("Gene symbol")
                continue
            if len(parts) > idx_sym:
                pid = parts[idx_id].strip()
                sym = parts[idx_sym].strip()
                if pid and sym:
                    # AnnotGPL sometimes lists multiple symbols sep by ///
                    m[pid] = sym.split("///")[0].strip()
    return m


def main():
    probe2sym = load_probe_symbol(ANNOT)
    print(f"GPL6885 annot: {len(probe2sym)} probe->symbol")
    manifest = []
    for gse, cfg in DATASETS.items():
        sm_path = os.path.join(WORK, f"{gse}_series_matrix.txt.gz")
        col_gsms, rows = parse_series_matrix(sm_path)
        n = len(col_gsms)
        print(f"\n=== {gse} (Acc{cfg['acc']}): {n} samples, {len(rows)} probes ===")
        # sample -> group
        sample_group = {}
        for grp, idxs in cfg["groups"].items():
            for i in idxs:
                if i < n:
                    sample_group[col_gsms[i]] = grp
        used = [g for g in col_gsms if g in sample_group]
        tag = gse
        # write exprs.tsv (probe x used samples)
        with open(os.path.join(WORK, f"{tag}_exprs.tsv"), "w") as out:
            out.write("ID_REF\t" + "\t".join(used) + "\n")
            col_idx = [col_gsms.index(g) for g in used]
            for pid, vals in rows.items():
                row = [vals[i] if i < len(vals) else "" for i in col_idx]
                if any(v == "" or v.lower() == "null" for v in row):
                    continue
                out.write(pid + "\t" + "\t".join(row) + "\n")
        with open(os.path.join(WORK, f"{tag}_groups.tsv"), "w") as out:
            out.write("sample\tgroup\n")
            for g in used:
                out.write(f"{g}\t{sample_group[g]}\n")
        with open(os.path.join(WORK, f"{tag}_contrasts.tsv"), "w") as out:
            out.write("\n".join(c for c, _ in cfg["contrasts"]) + "\n")
        # run limma engine
        r = subprocess.run([RSCRIPT, ENGINE, WORK, tag], capture_output=True, text=True)
        print(r.stdout.strip())
        if r.returncode != 0:
            print("R ERROR:\n", r.stderr[-1500:]); continue
        # map + write DEG csv per contrast
        for contrast, treat in cfg["contrasts"]:
            safe = "".join(ch if ch.isalnum() else "_" for ch in contrast)
            tt_path = os.path.join(WORK, f"{tag}__{safe}__topTable.tsv")
            if not os.path.exists(tt_path):
                print("  MISSING", tt_path); continue
            out_csv = os.path.join(STAGE, f"acc{cfg['acc']}_{cfg['region']}_{treat}.csv")
            n_map = 0
            with open(tt_path) as f, open(out_csv, "w", newline="") as o:
                rd = csv.DictReader(f, delimiter="\t")
                w = csv.writer(o)
                w.writerow(["gene_id", "gene_name", "log2FoldChange", "pvalue", "padj", "direction"])
                for row in rd:
                    pid = row["ID"]
                    sym = probe2sym.get(pid)
                    if not sym:
                        continue
                    lfc = float(row["logFC"]); padj = float(row["adj.P.Val"])
                    direction = "UP" if lfc > 0 else "DOWN"
                    w.writerow([pid, sym, lfc, row["P.Value"], padj, direction])
                    n_map += 1
            print(f"  {treat}: {n_map} probes mapped -> {os.path.basename(out_csv)}")
            manifest.append({"acc": cfg["acc"], "region": cfg["region"], "treatment": treat,
                             "gse": gse, "file": out_csv, "source": "limma_array"})
    json.dump(manifest, open(os.path.join(STAGE, "array_manifest.json"), "w"), indent=1)
    print(f"\nArray manifest: {len(manifest)} groups")


if __name__ == "__main__":
    main()
