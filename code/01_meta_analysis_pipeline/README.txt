================================================================
 FINALIZED BIOINFORMATICS WORKFLOW
 GLP-1R agonist / multi-agonist CNS bulk-RNA meta-analysis
================================================================

WHAT THIS IS
------------
A consensus differential-expression + functional-enrichment +
network pipeline across the GLP-1R-agonist accessions, restricted to
(a) GLP-1R-agonist and multi-agonist (dual / conjugate / combo) arms
and (b) central-nervous-system tissue only. Peripheral tissues
(liver, kidney, lung, heart, spleen, colon, adipose, muscle, WBCs in
Acc 8) and non-GLP arms (MK-801 alone, PYY alone, rapamycin, insulin,
vehicle-only) are excluded. Accession 3 is GLP-1R agonist and is kept.

GROUP DEFINITION
----------------
Each distinct treatment-vs-control comparison within a region is its
own GROUP (groups are NOT merged across regions or accessions). The
analysis is built from 35 groups (see Data/group_catalog.csv and
Data/master_group_legend.txt) across accessions 3, 4, 5, 6, 7, 8, 10,
11, 12. Accession 1 was excluded per user request. Accession 7 (lumbar
spinal cord, liraglutide vs EAE) is now INCLUDED, but that comparison
has zero significant DEGs (min padj ~0.999), so it contributes no genes
to the consensus and appears as a blank column on the heatmaps.

GENE STANDARDIZATION (all to MOUSE symbols)
-------------------------------------------
  Ensembl mouse (ENSMUSG, Acc 5,6)  -> symbol via mygene.info
  rat/mouse symbol columns (Acc 3,4,7,8) -> mouse title-case
  Illumina arrays (Acc 10,11,12)    -> see "ARRAYS" below
Matching for consensus/overlap is case-insensitive (UPPER key); the
displayed name is mouse title-case.

ARRAYS — GEO2R / limma re-analysis (Acc 10, 11, 12; all GPL6885)
----------------------------------------------------------------
Per request, the three Illumina GPL6885 microarray datasets were
re-derived with the GEO2R limma workflow (R 4.5.3 + Bioconductor
limma 3.66.0):
  log2 auto-transform -> lmFit(~0+group) -> makeContrasts ->
  eBayes(0.01) -> topTable(adjust="fdr").
Probe -> gene symbol uses the GPL6885 AnnotGPL annotation (the same
"Gene symbol" field GEO2R uses with AnnotGPL=TRUE). Contrasts:
  GSE71850 (Acc10): Exendin4, Ex4_Pre_TBI, Ex4_Po_t_TBI each vs Control (blast cohort)
  GSE41345 (Acc11): Ex4, Ex4_mTBI each vs sham
  GSE39586 (Acc12): glp vs nORM, using the user-supplied sample
                    grouping gsms = "00000111XXXXXXXXXXXXXXX"
Scripts: limma_array_engine.R, 01b_build_arrays_limma.py.

CONSENSUS GENES
---------------
Consensus = genes that are significant DEGs in >= 2 distinct groups
(overlap across datasets), computed separately per threshold:
  LFC0.2 : padj < 0.05 AND |log2FC| >= 0.2   -> 4,770 consensus genes
  LFC1   : padj < 0.05 AND |log2FC| >= 1.0   ->   297 consensus genes

ENRICHMENT (Enrichr via gseapy, organism = mouse)
-------------------------------------------------
Databases: GO_Biological_Process_2023, GO_Molecular_Function_2023,
GO_Cellular_Component_2023 (= GO:All), KEGG_2021_Human, Reactome_2022.
Significant terms at padj<0.05: 298 at LFC0.2; few at LFC1 (3 — the
strict |LFC|>=1 + recurrence filter leaves a small, mostly KEGG set,
so the LFC1 dot/bar/cnet plots are intentionally sparse).

PLOTS  (all are transparent PNG + PDF; folder = type/threshold/format)
--------------------------------------------------------------------
  Plots/Dot plots/<LFC>/<PDF|PNG>/        enrichment dot plots per DB
  Plots/Bar plots/<LFC>/<PDF|PNG>/        enrichment bar plots per DB
                                          + top-20 consensus-gene bar plot
  Plots/Heatmaps/Consensus/<LFC>/         top-10 consensus genes (by group
                                          frequency) x groups; group_legend.txt
  Plots/Heatmaps/{BBB,Inflammatory,Survival}/<LFC>/   Task-4 mechanism
                                          heatmaps (top-10 relevant genes;
                                          BBB force-includes all MMPs + Mmp9)
  Plots/Cnet plots/Consensus/<LFC>/       cnet of all consensus genes
  Plots/Cnet plots/{BBB,Inflammatory,Survival}/<LFC>/ pathway-restricted
                                          cnets (gene sets from cnet_gene_lists.py)
<LFC> is "LFC0.2" or "LFC1".

HEATMAP / CNET STYLING
----------------------
Heatmaps: Y = genes, X = groups labelled A, B, C ...; each heatmap
folder has a group_legend.txt mapping letter -> group. Colour scale is
RdBu_r centred at 0 with a robust (90th-percentile) cap for readability;
exact log2FC is annotated in every cell.
Cnets: cnet_style.py reference style with ENLARGED gene dots and
THICKER edges (per request); gene nodes coloured by mean log2FC,
term squares coloured by database, transparent background.

FOLDER MAP
----------
  scripts/   all Python + R scripts (01b -> 07) + cnet_style.py,
             cnet_gene_lists.py, limma_array_engine.R
  Data/      master_deg_*, consensus_*, enrichment_*, gene_lfc_*,
             group_catalog.csv, master_group_legend.txt
  Plots/     as above

RUN ORDER
---------
  01b_build_arrays_limma.py   (GEO2R/limma for Acc10/11/12)
  02_build_master.py          (standardize + per-group DEGs, 2 thresholds)
  03_consensus.py             (>=2-group consensus)
  04_enrichment.py [LFC02|LFC1]
  05_enrichment_plots.py      (dot, bar, consensus heatmap)
  06_make_cnets.py            (consensus + 3 pathway cnets)
  07_pathway_heatmaps.py      (Task-4 mechanism heatmaps)

PENDING / EXCLUDED
------------------
EXCLUDED (per user request): all Accession 1 groups.
PENDING (still locked at last run — re-run 02..07 once released):
  Acc7 lumbar_spinal_cord LIRA_EAE.
The Acc3 accumbens GLP_1_MK_801 and Acc4 hypothalamus Lira groups,
previously locked, are now INCLUDED.
