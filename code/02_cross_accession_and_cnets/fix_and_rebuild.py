"""
Fix two issues:
1. Use MASTER data mean LFC for cnet gene coloring (not just consensus)
2. Remove groups with <5 significant DEGs from overlap analysis
"""
import pandas as pd
import numpy as np
import os

OUT = "/sessions/practical-ecstatic-mendel/mnt/outputs"

# ============================================================
# FIX 1: Remove low-DEG groups and rebuild consensus
# ============================================================
print("=== Fixing group overlap: removing groups with <5 DEGs ===")

master = pd.read_csv(os.path.join(OUT, 'final_master_deg.csv'))
master['group'] = master['treatment_class'] + '|' + master['region'].fillna('unknown')

# Count DEGs per group
group_counts = master.groupby('group')['symbol'].nunique()
valid_groups = group_counts[group_counts >= 5].index.tolist()
removed_groups = group_counts[group_counts < 5].index.tolist()

print(f"Total groups: {len(group_counts)}")
print(f"Valid groups (≥5 DEGs): {len(valid_groups)}")
print(f"Removed groups (<5 DEGs): {removed_groups}")

# Filter master to valid groups only
master_filtered = master[master['group'].isin(valid_groups)]
print(f"Records after filtering: {len(master_filtered)} (from {len(master)})")

# Rebuild consensus
gene_groups = master_filtered.groupby('symbol')['group'].nunique().reset_index()
gene_groups.columns = ['symbol', 'n_groups']
cons_new = gene_groups[gene_groups['n_groups'] >= 2].sort_values('n_groups', ascending=False)

# Add stats using ALL master data (for accurate LFC)
gene_stats = master_filtered.groupby('symbol').agg(
    mean_lfc=('log2FC', 'mean'),
    median_lfc=('log2FC', 'median'),
    n_comparisons=('comparison', 'nunique'),
    n_accessions=('accession', 'nunique'),
    n_up=('direction', lambda x: (x == 'UP').sum()),
    n_down=('direction', lambda x: (x == 'DOWN').sum()),
).reset_index()

cons_new = cons_new.merge(gene_stats, on='symbol')
cons_new['predominant_direction'] = np.where(cons_new['n_up'] >= cons_new['n_down'], 'UP', 'DOWN')

print(f"New consensus (LFC≥0.5, ≥2 valid groups): {len(cons_new)} genes")
cons_new.to_csv(os.path.join(OUT, 'final_consensus_LFC05.csv'), index=False)
master_filtered.to_csv(os.path.join(OUT, 'final_master_deg.csv'), index=False)

# ============================================================
# FIX 2: Build a complete LFC lookup from ALL data (for cnet coloring)
# ============================================================
print("\n=== Building complete gene LFC lookup for cnet coloring ===")

# Use ALL records in master (not just consensus) to get accurate mean LFC per gene
all_gene_lfc = master_filtered.groupby('symbol')['log2FC'].mean().to_dict()
print(f"LFC lookup: {len(all_gene_lfc)} genes")

# Save this lookup
lfc_df = pd.DataFrame(list(all_gene_lfc.items()), columns=['symbol', 'mean_lfc'])
lfc_df.to_csv(os.path.join(OUT, 'final_gene_lfc_lookup.csv'), index=False)

# Also fix LFC0
print("\n=== Fixing LFC0 consensus ===")
master_lfc0 = pd.read_csv(os.path.join(OUT, 'final_master_deg_LFC0.csv'))
# Remove same invalid groups
master_lfc0_valid = master_lfc0[master_lfc0['group'].isin(
    master_lfc0.groupby('group')['symbol'].nunique().pipe(lambda x: x[x>=5]).index
)]
gene_groups_0 = master_lfc0_valid.groupby('symbol')['group'].nunique().reset_index()
gene_groups_0.columns = ['symbol', 'n_groups']
cons0_new = gene_groups_0[gene_groups_0['n_groups'] >= 2].sort_values('n_groups', ascending=False)

stats0 = master_lfc0_valid.groupby('symbol').agg(
    mean_lfc=('log2FC', 'mean'),
    n_up=('log2FC', lambda x: (x>0).sum()),
    n_down=('log2FC', lambda x: (x<0).sum())
).reset_index()
cons0_new = cons0_new.merge(stats0, on='symbol')
cons0_new['predominant_direction'] = np.where(cons0_new['n_up']>=cons0_new['n_down'],'UP','DOWN')
print(f"New LFC0 consensus: {len(cons0_new)} genes")
cons0_new.to_csv(os.path.join(OUT, 'final_consensus_LFC0.csv'), index=False)

# LFC0 lookup
all_gene_lfc0 = master_lfc0_valid.groupby('symbol')['log2FC'].mean().to_dict()
lfc0_df = pd.DataFrame(list(all_gene_lfc0.items()), columns=['symbol', 'mean_lfc'])
lfc0_df.to_csv(os.path.join(OUT, 'final_gene_lfc_lookup_LFC0.csv'), index=False)

print("\nDone! Ready to rebuild cnets with corrected LFC coloring.")
