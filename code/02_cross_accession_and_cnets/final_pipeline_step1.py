"""
STEP 1: Load all DEG files, standardize gene IDs, build master catalog
"""
import pandas as pd
import numpy as np
import os, glob, gzip, json

BASE = "/sessions/practical-ecstatic-mendel/mnt/Bulk RNA sequencing"
OUT = "/sessions/practical-ecstatic-mendel/mnt/outputs"

# ============================================================
# 1A. Build GPL6885 probe-to-symbol mapping
# ============================================================
print("=== Building GPL6885 probe-to-symbol map ===")
bgx_path = os.path.join(BASE, "GPL6885_MouseRef-8_V2_0_R0_11278551_A.bgx.gz")
probe_map = {}
with gzip.open(bgx_path, 'rt') as f:
    in_probes = False
    header = None
    for line in f:
        line = line.strip()
        if line == '[Probes]':
            in_probes = True
            continue
        if in_probes and header is None:
            header = line.split('\t')
            probe_idx = header.index('Probe_Id')
            sym_idx = header.index('Symbol')
            continue
        if in_probes:
            if line.startswith('['):
                break
            parts = line.split('\t')
            if len(parts) > max(probe_idx, sym_idx):
                pid = parts[probe_idx]
                sym = parts[sym_idx]
                if pid and sym and sym.strip():
                    probe_map[pid] = sym.strip()

print(f"  Mapped {len(probe_map)} probes to symbols")

# ============================================================
# 1B. Build Ensembl-to-symbol mapping via mygene
# ============================================================
print("\n=== Building Ensembl-to-symbol map ===")
# First collect all Ensembl IDs from Acc 3, 5, 6
ensembl_ids = set()
for acc in [3, 5, 6]:
    acc_dir = os.path.join(BASE, f"Accsession {acc}")
    if not os.path.isdir(acc_dir):
        continue
    for csv_file in glob.glob(os.path.join(acc_dir, "DEG_analysis/work/results/*_significant.csv")):
        df = pd.read_csv(csv_file, nrows=0)
        df_full = pd.read_csv(csv_file)
        ids = df_full['gene_id'].dropna().unique()
        ensembl_ids.update([x for x in ids if str(x).startswith('ENSMUSG')])

print(f"  Found {len(ensembl_ids)} unique Ensembl IDs")

import mygene
mg = mygene.MyGeneInfo()
ens_list = list(ensembl_ids)

# Query in batches
ens_map = {}
batch_size = 1000
for i in range(0, len(ens_list), batch_size):
    batch = ens_list[i:i+batch_size]
    results = mg.querymany(batch, scopes='ensembl.gene', fields='symbol', species='mouse', returnall=True)
    for hit in results['out']:
        if 'symbol' in hit and 'query' in hit:
            ens_map[hit['query']] = hit['symbol']
    print(f"  Batch {i//batch_size + 1}: mapped {len(ens_map)} so far")

print(f"  Total Ensembl mapped: {len(ens_map)}")

# ============================================================
# 1C. Load all DEG files with gene ID standardization
# ============================================================
print("\n=== Loading all DEG files ===")

# Accession metadata
acc_info = {
    3: {'species': 'mouse', 'id_type': 'ensembl', 'gse': 'GSE162614'},
    4: {'species': 'rat', 'id_type': 'symbol', 'gse': 'GSE190218'},
    5: {'species': 'mouse', 'id_type': 'ensembl', 'gse': 'GSE155178'},
    6: {'species': 'mouse', 'id_type': 'ensembl', 'gse': 'GSE184435'},
    7: {'species': 'mouse', 'id_type': 'symbol', 'gse': 'GSE106543'},
    8: {'species': 'mouse', 'id_type': 'symbol', 'gse': 'GSE113071'},
    10: {'species': 'mouse', 'id_type': 'illumina', 'gse': 'GSE71850'},
    11: {'species': 'mouse', 'id_type': 'illumina', 'gse': 'GSE41345'},
    12: {'species': 'mouse', 'id_type': 'illumina', 'gse': 'GSE39586'},
}

# Treatment classification
treatment_classes = {
    # Acc 3
    'GLP_1': 'GLP1_Mono', 'GLP_1_MK_801': 'Conjugate', 'conjugate_2': 'Conjugate', 'MK_801': 'Other',
    # Acc 4
    'Lira': 'GLP1_Mono', 'PYY_Lira': 'Dual_Agonist', 'PYY': 'Other',
    'BWM': 'Other', 'Bypass': 'Other', 'Lean': 'Other', 'Sham': 'Other', 'Saline': 'Other',
    # Acc 5
    'IP118_PY115': 'Dual_Agonist',
    # Acc 6
    'AC3174': 'GLP1_Mono', 'AC710222': 'Dual_Agonist', 'COMBO': 'Combination',
    # Acc 7
    'LIRA_EAE': 'GLP1_Mono',
    # Acc 8
    'Ex': 'GLP1_Mono', 'rapa': 'Other', 'Ex_KD': 'KD_Experiment', 'ctrl_KD': 'KD_Experiment',
    # Acc 10
    'Exendin4': 'GLP1_Mono', 'Ex4_PreTBI': 'GLP1_Mono', 'Ex4_PostTBI': 'GLP1_Mono',
    'Ex4_Pre_TBI': 'GLP1_Mono', 'Ex4_Po_t_TBI': 'GLP1_Mono', 'Ex4_Po_tTBI': 'GLP1_Mono',
    'TBI': 'Other', 'Bla_t_mTBI': 'Other', 'blast_TBI': 'Other', 'mTBI': 'Other',
    # Acc 11
    'Ex4': 'GLP1_Mono', 'Ex4_mTBI': 'GLP1_Mono',
    # Acc 12
    'WT_Ex_4': 'GLP1_Mono', 'WT_GLP_1': 'GLP1_Mono', 'HD_Ex_4': 'GLP1_Mono', 'HD_GLP_1': 'GLP1_Mono',
    'WT_INS': 'Other', 'HD_INS': 'Other', 'HD_CTR': 'Other',
    'combined_WT_Ex_4': 'GLP1_Mono', 'combined_WT_GLP_1': 'GLP1_Mono',
    'combined_HD_Ex_4': 'GLP1_Mono', 'combined_HD_GLP_1': 'GLP1_Mono',
    'combined_WT_INS': 'Other', 'combined_HD_INS': 'Other', 'combined_HD_CTR': 'Other',
}

# Rat-to-mouse ortholog conversion (capitalize first letter for mouse convention)
def rat_to_mouse_symbol(sym):
    if sym and isinstance(sym, str) and len(sym) > 0:
        return sym[0].upper() + sym[1:].lower() if len(sym) > 1 else sym.upper()
    return sym

def standardize_gene(gene_id, acc_num):
    info = acc_info[acc_num]
    if info['id_type'] == 'ensembl':
        return ens_map.get(str(gene_id), None)
    elif info['id_type'] == 'illumina':
        return probe_map.get(str(gene_id), None)
    elif info['id_type'] == 'symbol':
        if info['species'] == 'rat':
            return rat_to_mouse_symbol(str(gene_id))
        return str(gene_id)
    return None

def parse_comparison(filename, acc_num):
    """Extract region, treatment, control from filename"""
    basename = os.path.basename(filename).replace('DEG_', '').replace('_significant.csv', '')
    # Pattern: region_treatment_vs_control
    parts = basename.split('_vs_')
    if len(parts) != 2:
        return None, None, None
    control = parts[1]
    prefix = parts[0]
    # Region is typically the first part(s)
    # For acc 10, format is: hippocampus_blast_Treatment
    # For acc 8: agedst_hippocampus_Treatment or young_hippocampus_Treatment
    # General approach: known regions
    known_regions = ['hippocampus', 'hypothalamus', 'accumbens', 'brainstem', 
                     'nts', 'pvn', 'arc', 'dvc', 'mbh', 'lumbar_spinal_cord',
                     'frontalcortex', 'kidney', 'lung', 'heart', 'liver',
                     'spleen', 'colon', 'adipose', 'skeletalmuscle', 'wbcs']
    
    region = None
    treatment = prefix
    
    for r in known_regions:
        if r in prefix:
            # Find where region ends and treatment begins
            idx = prefix.find(r) + len(r)
            region_part = prefix[:idx]
            treatment = prefix[idx:].lstrip('_')
            region = region_part
            break
    
    if region is None:
        # Try splitting by first underscore after known patterns
        region = prefix.split('_')[0]
        treatment = '_'.join(prefix.split('_')[1:])
    
    return region, treatment, control

all_data = []
file_catalog = []

for acc_num in sorted(acc_info.keys()):
    # Find DEG files - prefer work/results folder
    acc_dir_candidates = [
        os.path.join(BASE, f"Accsession {acc_num}"),
        os.path.join(BASE, f"Accsession {acc_num} "),  # trailing space
    ]
    acc_dir = None
    for c in acc_dir_candidates:
        if os.path.isdir(c):
            acc_dir = c
            break
    if acc_dir is None:
        print(f"  WARNING: No directory found for Acc {acc_num}")
        continue
    
    results_dir = os.path.join(acc_dir, "DEG_analysis/work/results")
    if not os.path.isdir(results_dir):
        results_dir = os.path.join(acc_dir, "DEG_analysis/csv/significant_only")
    
    files = glob.glob(os.path.join(results_dir, "*_significant.csv"))
    print(f"  Acc {acc_num} ({acc_info[acc_num]['gse']}): {len(files)} comparisons in {results_dir}")
    
    for f in files:
        region, treatment, control = parse_comparison(f, acc_num)
        
        # Classify treatment
        tclass = treatment_classes.get(treatment, None)
        if tclass is None:
            # Try partial matching
            for k, v in treatment_classes.items():
                if k in treatment:
                    tclass = v
                    break
        if tclass is None:
            tclass = 'Other'
        
        # Skip non-GLP1R comparisons (controls only)
        if tclass == 'Other':
            continue
        
        df = pd.read_csv(f)
        if df.empty:
            continue
        
        # Standardize gene names
        df['symbol'] = df['gene_id'].apply(lambda x: standardize_gene(x, acc_num))
        df = df.dropna(subset=['symbol'])
        df = df[df['symbol'] != '']
        
        if df.empty:
            continue
        
        # Capitalize first letter for consistency (mouse convention)
        df['symbol'] = df['symbol'].apply(lambda x: x[0].upper() + x[1:] if len(x) > 1 else x.upper())
        
        comparison = os.path.basename(f).replace('DEG_', '').replace('_significant.csv', '')
        
        for _, row in df.iterrows():
            all_data.append({
                'accession': acc_num,
                'gse': acc_info[acc_num]['gse'],
                'comparison': comparison,
                'region': region,
                'treatment': treatment,
                'treatment_class': tclass,
                'control': control,
                'symbol': row['symbol'],
                'log2FC': row['log2FoldChange'],
                'padj': row['padj'],
                'direction': row['direction'],
            })
        
        file_catalog.append({
            'accession': acc_num,
            'gse': acc_info[acc_num]['gse'],
            'file': os.path.basename(f),
            'region': region,
            'treatment': treatment,
            'treatment_class': tclass,
            'n_genes_mapped': len(df),
            'n_genes_original': len(pd.read_csv(f)),
        })

master_df = pd.DataFrame(all_data)
catalog_df = pd.DataFrame(file_catalog)

print(f"\n=== SUMMARY ===")
print(f"Total records: {len(master_df)}")
print(f"Unique genes: {master_df['symbol'].nunique()}")
print(f"Accessions: {sorted(master_df['accession'].unique())}")
print(f"Treatment classes: {master_df['treatment_class'].value_counts().to_dict()}")
print(f"\nComparisons per accession:")
print(catalog_df.groupby('accession')[['file']].count())

# Save
master_df.to_csv(os.path.join(OUT, 'final_master_deg.csv'), index=False)
catalog_df.to_csv(os.path.join(OUT, 'final_file_catalog.csv'), index=False)
print(f"\nSaved master_deg ({len(master_df)} rows) and catalog ({len(catalog_df)} entries)")
