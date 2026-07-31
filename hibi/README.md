# HIBI arm — neonatal hypoxic-ischemic brain injury

> **The HIBI transcriptomic work in this folder was carried out by Melody and
> Kate.**

The **injury-side** transcriptomic meta-analysis of this thesis: which genes and
pathways are consistently dysregulated after neonatal hypoxic-ischemic brain
injury (HIBI), analysed with the same pipeline, thresholds and figure styling as
the GLP-1R arm in the repository root — so the two are directly comparable.

Where the GLP-1R arm asks *what does GLP-1R agonism do to the CNS*, this arm asks
*what does the injury do*, and the comparison at the end asks *which of the
injury's changes does GLP-1R agonism reverse*.

## Datasets

Ten HIBI GEO accessions were collected (`source_pipeline/`); **six** carry a
completed per-comparison DEG table and feed the meta-analysis:

| GSE | Platform | Region / model | Comparisons |
|-----|----------|----------------|-------------|
| GSE144455 | Agilent 2-colour array | Forebrain | HI vs naive, 3 h |
| GSE144456 | Agilent 2-colour array | Forebrain | P5 HI vs control, 3 h |
| GSE23317  | Illumina GPL6885 array | Cortex | HI vs Sham, 3 h |
| GSE23319  | Illumina GPL6885 array | Striatum | HI vs Sham, 3 h |
| GSE236133 | RNA-seq | Hippocampus | ipsi vs contra; WT / NEIL1KO / NEIL2KO × 3 h / 6 h |
| GSE238220 | RNA-seq (sorted microglia) | Hippocampus | HI vs control, 1 d / 3 d |

That yields **12 groups** (accession × region × comparison), matching the
"one row per distinct comparison" structure the GLP-1R arm uses:

| Group | n sig genes (LFC0.2) |
|---|---|
| GSE144455 forebrain HI vs naive 3h | 264 |
| GSE144456 forebrain P5 HI vs control 3h | 7 |
| GSE23317 cortex HI vs Sham 3h | 536 |
| GSE23319 striatum HI vs Sham 3h | 36 |
| GSE236133 hippocampus WT ipsi/contra 3h / 6h | 1732 / 2364 |
| GSE236133 hippocampus NEIL1KO ipsi/contra 3h / 6h | 1277 / 1747 |
| GSE236133 hippocampus NEIL2KO ipsi/contra 3h / 6h | 2632 / 3709 |
| GSE238220 hippocampus microglia HI vs control 1d / 3d | 4002 / 4444 |

**GSE97299, GSE294909, GSE312452 and GSE36215 are not in the meta-analysis.**
Their `source_pipeline/scripts/` and figures are kept for provenance, but no DEG
table was saved for them and their raw inputs are not in the repository, so they
could not be folded in without re-running the analysis from GEO. Two of the four
(GSE97299, GSE294909) are rat rather than mouse.

## Method

Identical to the GLP-1R arm, so the two consensus sets are comparable:

- DEG threshold **padj < 0.05 and |log2FC| ≥ 0.2** (LFC0.2).
- Gene symbols standardised to mouse title-case (`Fos`, not `FOS`).
- Per group, duplicate probes/rows for a gene collapse to the most significant
  (min padj, then max |log2FC|).
- **Consensus = significant in ≥ 2 groups**, as in `03_consensus.py`.
- GSE23317 / GSE23319 arrived as Illumina probe IDs with no local symbol
  annotation; mapped through the GPL6885 platform annotation table from NCBI GEO.

## Results

**5,346 consensus genes**, topping out at 11 of 12 groups (`Maff`, `Cxcl1`,
`Slc2a1`). The signature is an immediate-early and cytokine response — the whole
AP-1 family (`Fos`, `Fosb`, `Jun`, `Junb`, `Jund`), `Egr2`, `Nr4a1`, `Socs3`,
`Nfil3`, `Bhlhe40`, `Tnfrsf12a` — almost all strongly upregulated.

Against the four curated mechanism gene sets shared with the GLP-1R arm:

| Pathway | In HIBI consensus | Dominant direction |
|---|---|---|
| BBB / MMP | 35 / 44 | Mixed — `Spp1`, `Fn1`, `Cdh5`, `Nfkb2` up; `Rock2`, `F11r`, `Cav1` down |
| JAK-STAT3 inflammatory | 22 / 33 | Mostly up — `Socs3`, `Il1r1`, `Nfkb2`, `Stat1`, `Cdkn1a` |
| PI3K-Akt pro-survival | 41 / 68 | Mostly down — `Fgfr1/2`, `Akt3`, `Hras`, `Foxo1` |
| Ion channel | 18 / 30 | Uniformly down — `Gria1/2`, `Grin1/2a`, `Gabra1`, `Gabrb3`, `Grm1/5` |

The cleanest single result: every glutamate/GABA receptor and voltage-gated
channel in the set is suppressed after injury, against upregulated inflammatory
signalling and downregulated PI3K-Akt survival signalling.

### Comparison with the GLP-1R arm

| | genes |
|---|---|
| HIBI consensus | 5,346 |
| GLP-1R consensus | 4,769 |
| **Shared** | **1,715** |
| — same direction (concordant) | 875 |
| — opposite direction (discordant) | 840 |

The discordant set is the one that matters for the neuroprotection hypothesis:
genes the injury moves one way and GLP-1R agonism moves back. The immediate-early
and stress genes `Fos`, `Fosb`, `Egr2`, `Nfil3`, `Bhlhe40` are all HI-induced but
GLP-1R-suppressed; `Fkbp5` is HI-suppressed but is the single strongest GLP-1R
consensus gene (up in 12 groups).

Full gene-by-gene table in `hibi_vs_glp1r_consensus_comparison.csv`; ranked
top-20 lists in `comparison_summary.txt`.

## Layout

```
hibi/
├── README.md                     # this file
├── FIGURE_INDEX.md               # every figure → generating script
├── figures/                      # 5 Cnets (PNG + PDF)
├── figures_goplot/               # 31 GOplot figures (PNG + PDF)
│   └── Consensus/ BBB/ Inflammatory/ Survival/ IonChannel/
├── code/
│   ├── 01_meta_analysis_pipeline/
│   │   ├── 02_03_build_master_and_consensus.py
│   │   ├── 06f_consensus_frequency_cnet.py     # consensus signature Cnet
│   │   ├── 06d_focused_pathway_cnet.py         # 4 restricted pathway Cnets
│   │   ├── 09_compare_vs_glp1r.py              # consensus overlap vs GLP-1R arm
│   │   ├── cnet_style.py                       # verbatim from the GLP-1R arm
│   │   └── cnet_gene_lists.py                  # verbatim from the GLP-1R arm
│   └── 04_goplot_chord/
│       ├── prepare_goplot_data.py
│       ├── goplot_analysis.R
│       └── README.txt
├── data/
│   ├── meta_analysis/            # master_deg, group_catalog, consensus
│   └── goplot/                   # circle_dat input tables
├── source_pipeline/              # the original per-dataset HIBI R pipeline
│   ├── scripts/<GSE>/            # 01_download → 08_cnetplot per accession
│   ├── metadata/                 # GEO sample metadata, annotations
│   ├── results/                  # per-dataset DE + enrichment tables
│   └── figures/                  # per-dataset volcano / cnet / enrichment plots
├── hibi_vs_glp1r_consensus_comparison.csv
└── comparison_summary.txt
```

`source_pipeline/` is the upstream, per-dataset work that produced the DEG tables
`code/01_meta_analysis_pipeline/` consumes. Serialized R objects (`.rds` —
ExpressionSets, DESeq2 `dds`, enrichment result objects, ~313 MB) are excluded by
`.gitignore`; they are regenerable by re-running the stage scripts, and their
derived CSV/TXT tables are tracked.

## Reproducing

- **Python** ≥ 3.10: `pandas`, `numpy`, `matplotlib`, `networkx`, `gseapy`
- **R** ≥ 4.3 with `GOplot` (chords); the `source_pipeline` scripts additionally
  need `limma`, `DESeq2`, `clusterProfiler`, `org.Mm.eg.db`, `ReactomePA`

```bash
# meta-analysis + consensus
python code/01_meta_analysis_pipeline/02_03_build_master_and_consensus.py
# Cnets
python code/01_meta_analysis_pipeline/06f_consensus_frequency_cnet.py
python code/01_meta_analysis_pipeline/06d_focused_pathway_cnet.py
# GOplot chords
python code/04_goplot_chord/prepare_goplot_data.py
Rscript code/04_goplot_chord/goplot_analysis.R
# comparison against the GLP-1R arm
python code/01_meta_analysis_pipeline/09_compare_vs_glp1r.py
```

Scripts take their paths from constants near the top of each file — point them at
this folder before running.

## Not yet built

Stages 02 (cross-accession Venn / UpSet / ORA overview) and 03 (ranked-list GSEA,
NES heatmaps and bubble plots) exist for the GLP-1R arm but have no HIBI
equivalent yet.

`data/goplot/consensus_terms.csv` is capped to the top 60 terms per GO category:
GOplot's `reduce_overlap()` is O(n²) and does not complete on the full 959-term
table. No figure is affected — the plots draw only the top 20/10/7 terms — and
the uncapped table is kept as `consensus_terms_FULL.csv`.
