#!/usr/bin/env python3
"""
Cross-accession pathway analysis pipeline over accessions 3-8 and 10-12.
Per-accession DEG CSVs -> Cross_Accession_Pathway_Analysis/Complete_Analysis/

Figures generated:
  01_venn_gene_overlaps         — 3-way Venn (Obesity/Neuro/Aging) + DEG bar chart
  02_upset_all_accessions       — UpSet plot of gene overlaps across all accessions
  03_pathway_heatmap_hallmarks  — GSEA Hallmark heatmap × 9 accessions (* FDR<0.25, ** FDR<0.05)
  04_pathway_heatmap_kegg       — GSEA KEGG heatmap × 9 accessions (top 40)
  05_ora_dotplot_hallmark       — ORA Hallmark dot plots per accession
  05b_ora_dotplot_kegg          — ORA KEGG dot plots per accession
  06_cross_condition_hallmark   — Shared Hallmark pathways across conditions
  07_cross_condition_kegg       — Shared KEGG pathways across conditions
  08_gsea_bubble_hallmark       — Bubble plot: mean NES × condition (top 25 Hallmarks)

Prerequisites:
  pip install gseapy mygene matplotlib-venn upsetplot seaborn --break-system-packages

Gene ID mapping:
  - Mouse Ensembl IDs (ENSMUSG*) → gene symbols via Ensembl REST API
  - Illumina probe IDs (ILMN_*)  → gene symbols via GPL6885 annotation (downloaded from GEO)
  - Mouse gene symbols used as-is for Acc3, Acc4, Acc7, Acc8
  - Acc7 (Neuroinflammation/EAE) had 0 significant DEGs at padj<0.05, |LFC|>0.2

Gene sets:
  - Hallmark: MSigDB mouse Hallmarks (mh.all.v2023.2.Mm.symbols.gmt)
  - KEGG: KEGG_2019_Mouse via gseapy/Enrichr (uppercase gene symbols)

GSEA: gseapy prerank, 100 permutations, seed=42
ORA:  Local hypergeometric test with Benjamini-Hochberg FDR correction
"""

import os, sys, glob, pickle, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import seaborn as sns
import gseapy as gp
from scipy import stats
from scipy.stats import hypergeom
from collections import defaultdict
import urllib.request, json, time, io, gzip

warnings.filterwarnings('ignore')
plt.rcParams.update({'font.family': 'DejaVu Sans'})

MNT = '/path/to/Bulk RNA sequencing/'   # <-- UPDATE THIS PATH
OUT = '/path/to/Cross_Accession_Pathway_Analysis/Complete_Analysis/'  # <-- UPDATE
os.makedirs(OUT, exist_ok=True)

CNS_KW = ['hippo','frontal','hypothal','brainstem','accumbens',
           'arc','nts','pvn','dvc','mbh','spinal','lumbar']

ACC_CFG = {
    'Acc3' : {'dir': MNT+'Accsession 3/DEG_analysis/csv/',        'condition': 'Obesity/GLP-1R',    'cns_filter': False},
    'Acc4' : {'dir': MNT+'Accsession 4 /DEG_analysis/csv/',       'condition': 'Obesity/GLP-1R',    'cns_filter': False},
    'Acc5' : {'dir': MNT+'Accsession 5/DEG_analysis/csv/',        'condition': 'Obesity/GLP-1R',    'cns_filter': False},
    'Acc6' : {'dir': MNT+'Accsession 6/DEG_analysis/csv/',        'condition': 'Obesity/GLP-1R',    'cns_filter': False},
    'Acc7' : {'dir': MNT+'Accsession 7/DEG_analysis/csv/',        'condition': 'Neuroinflammation', 'cns_filter': False},
    'Acc8' : {'dir': MNT+'Accsession 8/DEG_analysis/csv/',        'condition': 'Aging',             'cns_filter': True},  # CNS tissues only
    'Acc10': {'dir': MNT+'Accsession 10 /DEG_analysis/work/results/', 'condition': 'TBI',           'cns_filter': False},
    'Acc11': {'dir': MNT+'Accsession 11/DEG_analysis/work/results/', 'condition': 'TBI',            'cns_filter': False},
    'Acc12': {'dir': MNT+'Accsession 12/DEG_analysis/work/results/', 'condition': 'Neurodegeneration', 'cns_filter': False},
}

ACC_ORDER = ['Acc3','Acc4','Acc5','Acc6','Acc7','Acc8','Acc10','Acc11','Acc12']
ACC_LABEL = {
    'Acc3':'Acc3\nGLP-1R Ag.','Acc4':'Acc4\nObesity/RYGB','Acc5':'Acc5\nPYY+GLP-1',
    'Acc6':'Acc6\nGLP-1R Ag.','Acc7':'Acc7\nEAE/MS','Acc8':'Acc8\nAging CNS',
    'Acc10':'Acc10\nTBI Blast','Acc11':'Acc11\nmTBI','Acc12':'Acc12\nHD',
}
ACC_COLOR = {
    'Acc3':'#E07B39','Acc4':'#F4A460','Acc5':'#DAA520','Acc6':'#CD853F',
    'Acc7':'#8E44AD','Acc8':'#27AE60','Acc10':'#2980B9','Acc11':'#1A5276','Acc12':'#C0392B',
}
COND_COLOR = {
    'Obesity/GLP-1R':'#E07B39','Neuroinflammation':'#8E44AD',
    'Aging':'#27AE60','TBI':'#2980B9','Neurodegeneration':'#C0392B',
}
SIG_PADJ = 0.05
SIG_LFC  = 0.2

def download_ilmn_map():
    """Download GPL6885 probe → mouse gene symbol map from GEO."""
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6885/annot/GPL6885.annot.gz'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    data = urllib.request.urlopen(req, timeout=60).read()
    probe_map = {}
    with gzip.open(io.BytesIO(data), 'rt') as f:
        hdr = False; headers = None
        for line in f:
            if line.startswith('!platform_table_begin'): hdr = True; continue
            if line.startswith('!platform_table_end'): break
            if hdr:
                cols = line.strip().split('\t')
                if headers is None: headers = cols; continue
                if len(cols) > 2 and cols[0].startswith('ILMN_') and cols[2].strip():
                    probe_map[cols[0]] = cols[2].strip()
    return probe_map

def fetch_ensm_map(ensm_ids, batch_size=200, sleep=0.15):
    """Map ENSMUSG IDs to gene symbols via Ensembl REST API."""
    result = {}
    ids = sorted(ensm_ids)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        data = json.dumps({'ids': batch}).encode()
        req = urllib.request.Request('https://rest.ensembl.org/lookup/id', data=data,
            headers={'Content-Type':'application/json','Accept':'application/json'})
        try:
            rd = json.loads(urllib.request.urlopen(req, timeout=30).read().decode())
            for eid, info in rd.items():
                if info and 'display_name' in info: result[eid] = info['display_name']
        except Exception as e: print(f"  Ensembl batch error: {e}")
        time.sleep(sleep)
    return result

# Resolve one DEG row to a mouse gene symbol, mapping ILMN_ probes and ENSMUSG IDs.
def std_gene(row, ilmn_map, ensm_map):
    gid = str(row.get('gene_id','')); gnm = str(row.get('gene_name',''))
    if gid.startswith('ILMN_'): return ilmn_map.get(gid, None)
    if gnm.startswith('ILMN_'): return ilmn_map.get(gnm, None)
    if gid.startswith('ENSMUSG'): return ensm_map.get(gid, None)
    if gnm.startswith('ENSMUSG'): return ensm_map.get(gnm, None)
    if gnm and gnm not in ('nan','None',''): return gnm
    return None

# Load and concatenate every accession's DEG CSVs with standardized gene symbols.
def load_all_data(ilmn_map, ensm_map):
    all_data = {}
    for acc_id, cfg in ACC_CFG.items():
        if not os.path.isdir(cfg['dir']): print(f"  {acc_id}: missing dir"); continue
        dfs = []
        for fp in sorted(glob.glob(cfg['dir']+'*.csv')):
            fname = os.path.basename(fp).lower()
            if any(x in fname for x in ['significant','summary']): continue
            if cfg['cns_filter'] and not any(kw in fname for kw in CNS_KW): continue
            try:
                raw = pd.read_csv(fp)
                if raw.empty or 'log2FoldChange' not in raw.columns: continue
                rows = []
                for _, row in raw.iterrows():
                    sym = std_gene(row, ilmn_map, ensm_map)
                    if not sym or sym in ('nan','None',''): continue
                    rows.append({'gene': sym,
                                 'lfc':  float(row['log2FoldChange']) if pd.notna(row.get('log2FoldChange')) else 0.0,
                                 'padj': float(row['padj'])           if pd.notna(row.get('padj'))           else 1.0,
                                 'group': os.path.basename(fp).replace('.csv','').replace('DEG_','')})
                if rows: dfs.append(pd.DataFrame(rows))
            except Exception: pass
        if dfs:
            df = pd.concat(dfs, ignore_index=True); all_data[acc_id] = df
            print(f"  {acc_id}: {df['gene'].nunique():,} genes | {df['group'].nunique()} groups")
    return all_data

# Build per-accession DEG sets (padj < 0.05, |log2FC| > 0.2), ranked lists, and background.
def build_gene_sets(all_data):
    acc_deg = {}; acc_rnk = {}; bg_all = set()
    for acc_id, df in all_data.items():
        sig = df[(df['padj'] < SIG_PADJ) & (df['lfc'].abs() > SIG_LFC)]
        acc_deg[acc_id] = set(sig['gene'].unique())
        bg_all |= set(df['gene'].unique())
        idx_max = df.groupby('gene')['lfc'].apply(lambda x: x.abs().idxmax())
        agg = df.loc[idx_max.values].drop_duplicates('gene').set_index('gene')['lfc']
        acc_rnk[acc_id] = agg[~agg.index.duplicated(keep='first')].sort_values(ascending=False)
    cond_deg = defaultdict(set)
    for acc_id, genes in acc_deg.items(): cond_deg[ACC_CFG[acc_id]['condition']] |= genes
    return acc_deg, acc_rnk, bg_all, dict(cond_deg)

# Run gseapy prerank GSEA against the Hallmark and KEGG gene sets.
def run_gsea(acc_rnk, hallmark_gmt, kegg_upper):
    gsea_hm = {}; gsea_kegg = {}
    for acc_id, rnk in acc_rnk.items():
        if len(rnk) < 10: continue
        print(f"  GSEA {acc_id}...", end=' ', flush=True)
        try:
            pre = gp.prerank(rnk=rnk, gene_sets=hallmark_gmt, min_size=5, max_size=500,
                             permutation_num=100, seed=42, outdir=None, verbose=False, threads=2)
            gsea_hm[acc_id] = pre.res2d
        except Exception as e: gsea_hm[acc_id] = pd.DataFrame(); print(f"[HM:{e}]",end='')
        rnk_u = rnk.copy(); rnk_u.index = rnk_u.index.str.upper(); rnk_u = rnk_u.groupby(level=0).max()
        try:
            pre_k = gp.prerank(rnk=rnk_u, gene_sets=kegg_upper, min_size=5, max_size=500,
                               permutation_num=100, seed=42, outdir=None, verbose=False, threads=2)
            gsea_kegg[acc_id] = pre_k.res2d
        except Exception as e: gsea_kegg[acc_id] = pd.DataFrame(); print(f"[KEGG:{e}]",end='')
        print("done")
    return gsea_hm, gsea_kegg

# Benjamini-Hochberg FDR correction.
def bh_fdr(pvals):
    pv=np.asarray(pvals,float); n=len(pv)
    if n==0: return np.array([])
    r=stats.rankdata(pv,'ordinal'); fdr=np.minimum(pv*n/r,1.0)
    idx=np.argsort(pv)[::-1]; fdr[idx]=np.minimum.accumulate(fdr[idx])
    return fdr

# Hypergeometric over-representation test of a gene list against gene sets.
def local_ora(gene_list, gene_sets, background, min_size=5, max_size=500):
    bg=set(background); gl=set(gene_list)&bg; N,K=len(bg),len(gl)
    if K==0: return pd.DataFrame()
    rows=[]
    for name,gs in gene_sets.items():
        s=set(gs)&bg; M=len(s)
        if not(min_size<=M<=max_size): continue
        k=len(gl&s)
        if k==0: continue
        pval=hypergeom.sf(k-1,N,M,K); OR=(k/(K-k+.1))/(M/(N-M+.1))
        rows.append({'Term':name,'k':k,'M':M,'K':K,'N':N,'pval':pval,'odds_ratio':OR,
                     'genes':','.join(sorted(gl&s))})
    if not rows: return pd.DataFrame()
    df=pd.DataFrame(rows); df['padj']=bh_fdr(df['pval'].values)
    df['neg_log10_pval']=-np.log10(df['pval']+1e-300)
    df['combined_score']=df['neg_log10_pval']*np.log1p(df['odds_ratio'])
    return df.sort_values('pval').reset_index(drop=True)

if __name__ == '__main__':
    print("Cross-Accession Pathway Analysis")
    print("="*50)
    print("\nStep 1: Building gene ID maps...")
    ilmn_map = download_ilmn_map()
    print(f"  ILMN map: {len(ilmn_map):,} probes")

    # Collect the ENSMUSG IDs that need mapping
    all_ensm = set()
    for cfg in ACC_CFG.values():
        for fp in glob.glob(cfg['dir']+'*.csv'):
            try:
                df_tmp = pd.read_csv(fp, usecols=['gene_id'])
                all_ensm.update(df_tmp['gene_id'][df_tmp['gene_id'].str.startswith('ENSMUSG',na=False)])
            except: pass
    ensm_map = fetch_ensm_map(all_ensm) if all_ensm else {}
    print(f"  ENSMUSG map: {len(ensm_map):,} genes")

    print("\nStep 2: Loading DEG data...")
    all_data = load_all_data(ilmn_map, ensm_map)

    print("\nStep 3: Building gene sets...")
    acc_deg, acc_rnk, bg_all, cond_deg = build_gene_sets(all_data)

    print("\nStep 4: Running GSEA...")
    hallmark_gmt_path = '/tmp/mh_hallmark_mouse.gmt'   # download from MSigDB
    kegg_gmt_path     = '/tmp/kegg_mouse.gmt'           # from gseapy Enrichr
    # See README for how to obtain these GMT files.
    def read_gmt(fp):
        gs={}
        with open(fp) as f:
            for line in f:
                p=line.strip().split('\t')
                if len(p)>=3: gs[p[0]]=[g for g in p[2:] if g]
        return gs
    hallmark_gs   = read_gmt(hallmark_gmt_path)
    kegg_gs_upper = {k:[g.upper() for g in v] for k,v in read_gmt(kegg_gmt_path).items()}
    gsea_hm, gsea_kegg = run_gsea(acc_rnk, hallmark_gs, kegg_gs_upper)

    print("\nStep 5: Running ORA...")
    bg_upper = {g.upper() for g in bg_all}
    ora_hm = {a: local_ora(list(g), hallmark_gs, bg_all) for a,g in acc_deg.items() if len(g)>=5}
    ora_kegg = {a: local_ora([x.upper() for x in g], kegg_gs_upper, bg_upper) for a,g in acc_deg.items() if len(g)>=5}

    print("\nStep 6: Generating figures...")
    print(f"  Output: {OUT}")
    # Figures 01-08 are produced by step4_plots.py / GenerateEnrichmentPlots.py.
    print("\nComplete!")
