# GLP-1R Bioinformatics Project

Transcriptomic **meta-analysis of GLP-1 receptor (GLP-1R) agonism across the central nervous system**, using nine public bulk-RNA / microarray GEO datasets. This repository is the cleaned, reproducible subset of a larger working directory: it contains the analysis code, the processed (non-raw) data tables the code consumes, the figures used in the thesis / presentation, and the written reports.

This meta-analysis is the transcriptomic arm (Aim 3 / "bioinformatic meta-analysis") of an MSc thesis on **semaglutide neuroprotection in neonatal hypoxic-ischemic brain injury**. It asks which pathways and consensus genes are consistently modulated by GLP-1R agonists across brain regions, disease contexts, and drug analogues — with a focused read-out on **blood-brain-barrier / MMP, JAK-STAT3 inflammatory, PI3K-Akt pro-survival, and ion-channel** programs.

## Datasets

Nine GEO accessions, spanning obesity, aging, neuroinflammation, neurodegeneration, and traumatic brain injury contexts. Raw count matrices are **not** included in this repo (they are re-downloadable from GEO by accession); the processed per-group DEG tables that feed the meta-analysis are under `data/`.

| Accession | GSE | Platform | Contexts / regions |
|-----------|-----|----------|--------------------|
| Acc3  | GSE244607 | RNA-seq (gene symbol) | Obesity — accumbens, brainstem, hypothalamus (GLP-1 ± MK-801, conjugates) |
| Acc4  | GSE190218 | RNA-seq (gene symbol) | Obesity — hypothalamus (Liraglutide, PYY+Lira) |
| Acc5  | GSE160802 | RNA-seq (Ensembl) | Obesity — ARC, NTS, PVN (IP118/PY115) |
| Acc6  | GSE135862 | RNA-seq (Ensembl) | Obesity — DVC, MBH (AC3174, AC710222, combo) |
| Acc7  | GSE186318 | RNA-seq (gene symbol) | Neuroinflammation — lumbar spinal cord (Liraglutide, EAE) |
| Acc8  | GSE280382 | RNA-seq (gene symbol) | Aging CNS — frontal cortex, hippocampus, hypothalamus (Exendin, short/long-term) |
| Acc10 | GSE71850  | Microarray | TBI — hippocampus (Exendin-4, pre/post-TBI) |
| Acc11 | GSE41345  | Microarray | mTBI — hippocampus (Exendin-4) |
| Acc12 | GSE39586  | Microarray | Neurodegeneration — hypothalamus (GLP) |

> Note: the presentation deck lists "GSE13596" — that is a transcription typo for **GSE135862** (Acc6), the accession used throughout the pipeline and data tables here.

DEGs were called per group (accession × region × treatment) and thresholded at **padj < 0.05** with log2 fold-change cutoffs of **|LFC| ≥ 0.2 / 0.5 / 1.0** (referred to as LFC0.2 / LFC0.5 / LFC1 in filenames). The LFC0.2 set is the primary one used for the figures. RNA-seq groups were analysed with a DESeq2-style workflow; microarray groups with limma/GEO2R (see `data/meta_analysis/master_group_legend.txt` for the platform of each group).

## Repository layout

```
GLP-1R bioinformatics project/
├── README.md                     # this file
├── FIGURE_INDEX.md               # every figure → source pipeline → generating script
├── figures/                      # figures used in the thesis / deck (PNG + PDF), Fig01–Fig11
├── code/
│   ├── 01_meta_analysis_pipeline/   # core pipeline: build → consensus → enrichment → cnets → heatmaps
│   ├── 02_cross_accession_and_cnets/# cross-accession comparison + pathway-restricted Cnet plots
│   ├── 03_gsea/                     # ranked-list GSEA (KEGG/Hallmark) + NES plots
│   └── 04_goplot_chord/             # GOplot chord / bubble diagrams (R)
├── data/                         # processed tables only — no raw matrices
│   ├── meta_analysis/               # master DEG, consensus, enrichment, gene-LFC lookups, group legend
│   ├── cnets/                       # consensus + curated gene lists used by the Cnet scripts
│   ├── cross_accession/             # cross-accession ORA/GSEA result tables + NES matrices
│   ├── gsea/                        # per-group .rnk inputs and GSEA/ORA result CSVs
│   └── goplot/                      # curated pathway gene/term lists for the chord plots
└── reports/                      # written analysis reports & methods write-ups (.docx)
```

## Pipeline overview

The four code stages are independent folders (each was a distinct working stage) rather than a single fused script. Broad order of execution:

1. **`01_meta_analysis_pipeline/`** — the backbone. `01b_build_arrays_limma.py` builds per-group arrays (microarray via `limma_array_engine.R`), `02_build_master.py` assembles the master DEG table, `03_consensus.py` derives consensus DEGs (genes significant in > 2 datasets), `04`/`05` run enrichment (KEGG, Reactome, GO-BP/CC/MF) and dot/bar plots, `06*` build the concept-network (Cnet) plots, `07` the pathway heatmaps, `08*` the per-group Cnet montages.
2. **`02_cross_accession_and_cnets/`** — `cross_accession_complete_pipeline.py` produces the cross-accession Venn / UpSet / ORA overview; the `cnet_*` and `step*` scripts build the pathway-restricted Cnet plots (BBB/MMP, JAK-STAT3 inflammatory, PI3K-Akt survival, ion channel).
3. **`03_gsea/`** — `01_build_rnk.py` → `02_run_gsea.py` → `03_make_plots.py`: ranked-list GSEA against KEGG and Hallmark gene sets, producing the NES heatmaps and bubble plots.
4. **`04_goplot_chord/`** — `prepare_goplot_data.py` shapes the enrichment output; `goplot_analysis.R` renders the GOChord / GOBubble diagrams.

Each stage folder keeps its own original `README` where one existed. Scripts take input/output paths as arguments or as paths near the top of the file — check the stage README and the script header before running, and point `data-dir` at the matching `data/` subfolder.

## Reproducing / re-running

- **Python** ≥ 3.10: `pandas`, `numpy`, `matplotlib`, `networkx`, `gseapy`, `mygene`.
- **R** ≥ 4.5: `limma`, `GEO2R`, and (for chord plots) `GOplot`.
- Raw GEO matrices must be downloaded separately from GEO by accession if you want to regenerate the per-group DEG tables in `data/`. Everything downstream of the DEG tables runs from the data included here.

## Provenance

Figures were selected to match those used in the thesis draft and the "Thesis figures" presentation deck; the mapping from each figure back to its generating stage and script is in [`FIGURE_INDEX.md`](FIGURE_INDEX.md). Exploratory / superseded analysis folders and the large raw inputs from the original working directory were intentionally left out to keep this repository focused and reviewable.
