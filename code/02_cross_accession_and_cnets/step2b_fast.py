# Builds the |log2FC| >= 0 consensus from the full DEG files (padj < 0.05 only).
# Per-accession DEG_*.csv -> final_consensus_LFC0.csv, final_master_deg_LFC0.csv

import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import pandas as pd, numpy as np, os, glob, gzip
import mygene

BASE = BASE + "/mnt/Bulk RNA sequencing"
OUT = BASE + "/mnt/outputs"

# GPL6885 (Illumina MouseRef-8 v2) probe-to-symbol map.
bgx_path = os.path.join(BASE, "GPL6885_MouseRef-8_V2_0_R0_11278551_A.bgx.gz")
probe_map = {}
with gzip.open(bgx_path, 'rt') as f:
    in_probes = False; header = None
    for line in f:
        line = line.strip()
        if line == '[Probes]': in_probes = True; continue
        if in_probes and header is None:
            header = line.split('\t')
            probe_idx = header.index('Probe_Id'); sym_idx = header.index('Symbol'); continue
        if in_probes:
            if line.startswith('['): break
            parts = line.split('\t')
            if len(parts) > max(probe_idx, sym_idx):
                pid, sym = parts[probe_idx], parts[sym_idx]
                if pid and sym and sym.strip(): probe_map[pid] = sym.strip()

# Ensembl IDs used by accessions 3, 5, and 6.
ens_ids = set()
for acc in [3, 5, 6]:
    acc_dir = os.path.join(BASE, f"Accsession {acc}")
    for f in glob.glob(os.path.join(acc_dir, "DEG_analysis/work/results/DEG_*.csv")):
        if 'summary' in f: continue
        try:
            col = pd.read_csv(f, usecols=['gene_id'], nrows=5000)
            ens_ids.update([x for x in col['gene_id'].dropna().unique() if str(x).startswith('ENSMUSG')])
        except: pass

mg = mygene.MyGeneInfo()
ens_map = {}
for i in range(0, len(list(ens_ids)), 1000):
    batch = list(ens_ids)[i:i+1000]
    res = mg.querymany(batch, scopes='ensembl.gene', fields='symbol', species='mouse', returnall=True)
    for hit in res['out']:
        if 'symbol' in hit and 'query' in hit: ens_map[hit['query']] = hit['symbol']
print(f"Maps: {len(probe_map)} probes, {len(ens_map)} ensembl")

acc_info = {3:'ensembl',4:'rat_symbol',5:'ensembl',6:'ensembl',7:'symbol',8:'symbol',10:'illumina',11:'illumina',12:'illumina'}

treatment_classes = {
    'GLP_1':'GLP1_Mono','GLP_1_MK_801':'Conjugate','conjugate_2':'Conjugate',
    'Lira':'GLP1_Mono','PYY_Lira':'Dual_Agonist',
    'IP118_PY115':'Dual_Agonist',
    'AC3174':'GLP1_Mono','AC710222':'Dual_Agonist','COMBO':'Combination',
    'Ex':'GLP1_Mono','Ex_KD':'KD_Experiment','ctrl_KD':'KD_Experiment',
    'Exendin4':'GLP1_Mono','Ex4_PreTBI':'GLP1_Mono','Ex4_PostTBI':'GLP1_Mono',
    'Ex4_Pre_TBI':'GLP1_Mono','Ex4_Po_t_TBI':'GLP1_Mono','Ex4_Po_tTBI':'GLP1_Mono',
    'Ex4':'GLP1_Mono','Ex4_mTBI':'GLP1_Mono','LIRA_EAE':'GLP1_Mono',
    'WT_Ex_4':'GLP1_Mono','WT_GLP_1':'GLP1_Mono','HD_Ex_4':'GLP1_Mono','HD_GLP_1':'GLP1_Mono',
    'combined_WT_Ex_4':'GLP1_Mono','combined_WT_GLP_1':'GLP1_Mono',
    'combined_HD_Ex_4':'GLP1_Mono','combined_HD_GLP_1':'GLP1_Mono',
}

known_regions = ['hippocampus','hypothalamus','accumbens','brainstem','nts','pvn','arc','dvc','mbh',
                 'lumbar_spinal_cord','frontalcortex','kidney','lung','heart','liver','spleen','colon',
                 'adipose','skeletalmuscle','wbcs']

# Map one accession's native gene ID to a mouse symbol.
def get_symbol(gid, acc):
    t = acc_info[acc]
    if t == 'ensembl': return ens_map.get(str(gid))
    if t == 'illumina': return probe_map.get(str(gid))
    if t == 'rat_symbol':
        s = str(gid)
        return s[0].upper() + s[1:].lower() if len(s)>1 else s.upper()
    return str(gid)

# Extract region and treatment from a DEG filename of the form region_treatment_vs_control.
def parse_comp(filename):
    bn = os.path.basename(filename).replace('DEG_','').replace('_significant.csv','').replace('.csv','')
    parts = bn.split('_vs_')
    if len(parts)!=2: return None, None
    prefix = parts[0]
    for r in known_regions:
        if r in prefix:
            idx = prefix.find(r)+len(r)
            return prefix[:idx], prefix[idx:].lstrip('_')
    return prefix.split('_')[0], '_'.join(prefix.split('_')[1:])

records = []  # (symbol, group, log2FC)
for acc_num in sorted(acc_info.keys()):
    acc_dirs = [os.path.join(BASE, f"Accsession {acc_num}"), os.path.join(BASE, f"Accsession {acc_num} ")]
    acc_dir = next((d for d in acc_dirs if os.path.isdir(d)), None)
    if not acc_dir: continue
    
    results_dir = os.path.join(acc_dir, "DEG_analysis/work/results")
    full_files = [f for f in glob.glob(os.path.join(results_dir, "DEG_*.csv")) 
                  if 'significant' not in f and 'summary' not in f]
    if not full_files:
        full_files = glob.glob(os.path.join(results_dir, "*_significant.csv"))
    
    n_processed = 0
    for fpath in full_files:
        region, treatment = parse_comp(fpath)
        if not treatment: continue
        tclass = treatment_classes.get(treatment)
        if not tclass:
            for k,v in treatment_classes.items():
                if k in treatment: tclass=v; break
        if not tclass: continue
        
        group = f"{tclass}|{region}"
        
        try:
            df = pd.read_csv(fpath, usecols=['gene_id','log2FoldChange','padj'])
        except:
            continue
        df = df[df['padj']<0.05].copy()
        if df.empty: continue
        
        df['sym'] = df['gene_id'].apply(lambda x: get_symbol(x, acc_num))
        df = df.dropna(subset=['sym'])
        df['sym'] = df['sym'].apply(lambda x: x[0].upper()+x[1:] if len(str(x))>1 else str(x).upper())
        
        for _, row in df.iterrows():
            records.append((row['sym'], group, row['log2FoldChange']))
        n_processed += 1
    
    print(f"  Acc {acc_num}: {n_processed} files processed")

print(f"\nTotal LFC0 records: {len(records)}")
df_all = pd.DataFrame(records, columns=['symbol','group','log2FC'])

# Consensus requires a gene in at least 2 treatment x region groups.
gene_groups = df_all.groupby('symbol')['group'].nunique().reset_index()
gene_groups.columns = ['symbol','n_groups']
cons = gene_groups[gene_groups['n_groups']>=2].sort_values('n_groups', ascending=False)

stats = df_all.groupby('symbol').agg(
    mean_lfc=('log2FC','mean'), n_up=('log2FC', lambda x: (x>0).sum()),
    n_down=('log2FC', lambda x: (x<0).sum())
).reset_index()
cons = cons.merge(stats, on='symbol')
cons['predominant_direction'] = np.where(cons['n_up']>=cons['n_down'],'UP','DOWN')

print(f"LFC0 consensus (≥2 groups): {len(cons)} genes")
cons.to_csv(os.path.join(OUT, 'final_consensus_LFC0.csv'), index=False)
df_all.to_csv(os.path.join(OUT, 'final_master_deg_LFC0.csv'), index=False)
print("Done!")
