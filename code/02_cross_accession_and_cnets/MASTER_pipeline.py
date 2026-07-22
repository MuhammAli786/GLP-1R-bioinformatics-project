"""
=============================================================================
GLP-1R AGONIST META-ANALYSIS: COMPLETE PIPELINE
=============================================================================
Accessions: 3-8, 10-12 (excluding 1, 9)
Species: Mouse (Acc 4 = rat, converted to mouse orthologs)
Gene ID types: Ensembl (Acc 3,5,6), Gene Symbols (Acc 4,7,8), Illumina probes (Acc 10,11,12)
Platform for probes: GPL6885 (Illumina MouseRef-8 v2)

Pipeline steps:
  1. Load DEG files + standardize gene IDs
  2. Build group-overlap consensus (treatment×region, freq≥2)
  3. Run functional enrichment (GO BP/MF/CC, KEGG, Reactome)
  4. Generate dot plots, barplots, heatmaps
  5. Build cnet plots (greedy set-cover): All, JAK-STAT3, BBB, Akt/PI3K, Combined

Two thresholds:
  - LFC≥0.5: Uses *_significant.csv files (pre-filtered |LFC|≥0.5, padj<0.05)
  - LFC≥0: Uses full DEG files (padj<0.05 only)

Cnet plot style:
  - Font size 20 for all labels
  - Edge width 2.5, alpha 0.35
  - Spring layout k=4.0, iterations=150, seed=42
  - Figure size 30×26, transparent background
  - RdBu_r colormap for log2FC
  - Term labels: white on colored boxes (Set3 palette)
  - Gene labels: white on dark boxes
  - Greedy set-cover term selection (max 15 terms, maximize gene coverage)
=============================================================================
"""

# This is the documentation header. Individual step scripts:
# - final_pipeline_step1.py: Gene ID standardization + data loading
# - step2b_fast.py: LFC≥0 consensus building (full DEG files)
# - step3_enrichment.py: Enrichr analysis (5 databases)
# - step4_plots.py: Dot plots, barplots, heatmaps
# - step5_cnets.py: Cnet plot generation

# PATHWAY GENE LISTS (from GLP1R_Pathway_Gene_Tables.docx):
JAK_STAT3_GENES = [
    'Gfap','Il1r1','Il6ra','Il6st','Lifr','Osmr','Mapk1','Fyn','Akt1','Jak1',
    'Nfkbia','Nfkb2','Ikbkg','Pias2','Pias1','Ptpn11','Ptprd','Ptprt','Egfr',
    'Erbb2','Fgfr1','Socs5','Socs3','Socs4','Socs6','Stat1','Cdkn1a','Bcl2l1',
    'Vegfa','Hif1a','Mmp2','Mcl1','Ccnd1'
]  # 33 genes

BBB_MMP_GENES = [
    'Adam10','Adam9','Adam15','Cdh5','Ctnnb1','Ctnnd1','Vegfa','Hif1a','Aqp4',
    'Gfap','Slc2a1','Lrp1','Slc7a5','Slc16a1','Tfrc','Abcb1b','Nid1','Lamc1',
    'Fn1','Spp1','Kdr','Cav1','Vwf','Pecam1','Eng','Nos3','Mmp2','Mmp16',
    'Mmp14','Mmp25','Nfkbia','Nfkb2','Mapk1','Rock2','Fyn','Rock1','Akt1',
    'Cldn5','Tjp1','Jam2','F11r','Cldn12','Ocln','Tjp2'
]  # 44 genes

AKT_PI3K_GENES = [
    'Irs1','Irs2','Sos1','Gab1','Grb2','Pdpk1','Akt3','Akt1','Bcl2l1','Bad',
    'Mcl1','Bax','Xiap','Birc2','Ulk2','Map1lc3a','Becn1','Ulk1','Creb1',
    'Creb5','Creb3','Cdkn1a','Ccnd1','Foxo1','Foxo3','Foxo4','Foxo6','Gsk3a',
    'Gsk3b','Vegfa','Kdr','Insr','Igf1r','Egfr','Erbb2','Erbb3','Fgfr1',
    'Fgfr2','Glp1r','Nfkbia','Nfkb2','Phlpp1','Phlpp2','Pik3ca','Pik3r1',
    'Pik3r3','Pik3c2a','Pik3cb','Pik3cg','Pik3c2b','Pten','Hras','Kras',
    'Mapk1','Map2k2','Nras','Raf1','Fyn','Mtor','Tsc2','Rptor','Rictor',
    'Rheb','Tsc1','Eif4g1','Eif4ebp1','Rps6kb1','Rps6kb2'
]  # 68 genes

COMBINED_GENES = list(set(JAK_STAT3_GENES + AKT_PI3K_GENES))  # 88 unique

# TREATMENT CLASSIFICATIONS:
# GLP1_Mono: Liraglutide, Exendin-4, GLP-1, Ex4 (±TBI context)
# Dual_Agonist: PYY+Lira, IP118/PY115, AC710222
# Conjugate: GLP-1/MK-801 conjugates
# Combination: Multi-drug combos (COMBO)
# KD_Experiment: Knockdown experiments

# ACCESSION INFO:
# Acc 3 (GSE162614): Mouse, hypothalamus/accumbens/brainstem, GLP-1 conjugates
# Acc 4 (GSE190218): Rat→Mouse, hypothalamus, Liraglutide/PYY
# Acc 5 (GSE155178): Mouse, NTS/PVN/ARC, dual agonist IP118/PY115
# Acc 6 (GSE184435): Mouse, DVC/MBH, AC3174/AC710222/COMBO
# Acc 7 (GSE106543): Mouse, spinal cord, Lira+EAE (empty - no sig DEGs)
# Acc 8 (GSE113071): Mouse, multi-tissue, Exendin-4 aging study
# Acc 10 (GSE71850): Mouse, hippocampus, Ex4 ± blast TBI
# Acc 11 (GSE41345): Mouse, hippocampus, Ex4 ± mTBI
# Acc 12 (GSE39586): Mouse, hypothalamus, Ex4/GLP-1 WT/HD
