# GSEA Analysis — GLP-1R agonist CNS meta-analysis

GSEA / enrichment across the **GLP-1R-agonist (incl. multi-agonist) groups
that have DEGs** — **31 of 35**. Excluded: the two empty groups (Acc7 spinal
cord, Acc8 hippocampus LT_A_Ex) **and AC710222** (Acc6 dvc + mbh), which is
the **CCK1R agonist** in the GLP-1/CCK study — not a GLP-1R agonist (AC3174 =
the GLP-1R agonist and COMBO = GLP-1R+CCK are kept). Gene-set libraries:
**MSigDB Hallmark 2020** and **KEGG 2021**.

## Method

1. **Ranked lists** (`01_build_rnk.py`): one full ranked `.rnk` per group —
   *every* gene (not just DEGs), standardized to mouse symbols then
   UPPER-cased to match the human Hallmark/KEGG libraries (ortholog-by-name),
   ranked by log2FC. → `data/rnk/<group>.rnk` (33 files).
2. **GSEA prerank** (`02_run_gsea.py`): `gseapy.prerank` per group × library
   (NES, nominal p, FDR). → `data/gsea/<group>__{hallmark,kegg}.csv` (66 files).
3. **ORA** on the consensus DEGs (genes significant in ≥2 groups) via Enrichr.
   → `data/ora/consensus_ora_{hallmark,kegg}.csv`.
4. **Figures** (`03_make_plots.py`).

## Condition groups (by accession)

| Condition | Accessions | # groups |
|-----------|-----------|----------|
| Obesity | 3, 4, 5, 6 | 16 |
| Aging | 8 | 9 |
| TBI | 10, 11 | 5 |
| Neurodegeneration | 12 | 1 |

(AC710222 / CCK1R-agonist groups excluded → Obesity 16, total 31.)

## Figures

| File | What it shows |
|------|---------------|
| `01_venn_gene_overlaps.png` | DEG overlap across Obesity / Aging / TBI conditions |
| `02_upset_all_accessions.png` | DEG set intersections across the 8 accessions |
| `03_pathway_heatmap_hallmark.png` | GSEA NES, top-variable Hallmark pathways × groups (condition colour bar) |
| `04_pathway_heatmap_kegg.png` | GSEA NES, top-variable KEGG pathways × groups |
| `05_ora_dotplot_hallmark.png` | ORA of consensus DEGs vs Hallmark (dot size = genes, colour = -log10 padj) |
| `05b_ora_dotplot_kegg.png` | ORA of consensus DEGs vs KEGG |
| `06_cross_condition_hallmark.png` | Mean NES per condition × Hallmark pathway |
| `07_cross_condition_kegg.png` | Mean NES per condition × KEGG pathway |
| `08_gsea_bubble_hallmark.png` | Summary: mean Hallmark NES across groups (size = # groups significant) |

## Headline result (08 / bubble)

Consistently **up** across GLP-1R-agonist CNS groups: Oxidative
Phosphorylation, Fatty-Acid Metabolism, Adipogenesis, ROS pathway, DNA
repair, Peroxisome. Consistently **down**: TNF-α/NF-κB signalling,
Angiogenesis, PI3K/Akt/mTOR, Wnt/β-catenin, TGF-β, Hedgehog, Glycolysis —
i.e. a pro-oxidative-metabolism, anti-inflammatory shift, matching the
consensus DEG / c-net findings.

## Reproduce

```bash
python scripts/01_build_rnk.py        # ranked lists
python scripts/02_run_gsea.py         # prerank (run until ALL_DONE)
python scripts/03_make_plots.py       # figures
```

## Notes

* GSEA prerank used 100 permutations (fast) — NES/ranking are robust; treat
  FDR as approximate. Increase `permutation_num` in `02_run_gsea.py` for
  publication-grade q-values.
* Hallmark/KEGG are human libraries; mouse genes are matched by upper-casing
  the symbol (standard ortholog-by-name approximation).
