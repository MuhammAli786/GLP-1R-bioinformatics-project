"""
Rebuilds the consensus tables after dropping groups with fewer than 5 significant
DEGs, and writes a whole-master mean log2FC lookup used to colour cnet genes.
final_master_deg*.csv -> final_consensus_*.csv, final_gene_lfc_lookup*.csv
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import pandas as pd
import numpy as np
import os

OUT = BASE + "/mnt/outputs"

# Drop groups with fewer than 5 DEGs, then rebuild the |LFC| >= 0.5 consensus.
print("=== Fixing group overlap: removing groups with <5 DEGs ===")

master = pd.read_csv(os.path.join(OUT, 'final_master_deg.csv'))
master['group'] = master['treatment_class'] + '|' + master['region'].fillna('unknown')

group_counts = master.groupby('group')['symbol'].nunique()
valid_groups = group_counts[group_counts >= 5].index.tolist()
removed_groups = group_counts[group_counts < 5].index.tolist()

print(f"Total groups: {len(group_counts)}")
print(f"Valid groups (≥5 DEGs): {len(valid_groups)}")
print(f"Removed groups (<5 DEGs): {removed_groups}")

master_filtered = master[master['group'].isin(valid_groups)]
print(f"Records after filtering: {len(master_filtered)} (from {len(master)})")

# Consensus requires a gene in at least 2 valid groups.
gene_groups = master_filtered.groupby('symbol')['group'].nunique().reset_index()
gene_groups.columns = ['symbol', 'n_groups']
cons_new = gene_groups[gene_groups['n_groups'] >= 2].sort_values('n_groups', ascending=False)

# Summary statistics over all filtered records, not only the consensus genes.
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

# Gene-level mean log2FC lookup used for cnet node colouring.
print("\n=== Building complete gene LFC lookup for cnet coloring ===")

# Averaged over all master records, not only consensus genes.
all_gene_lfc = master_filtered.groupby('symbol')['log2FC'].mean().to_dict()
print(f"LFC lookup: {len(all_gene_lfc)} genes")

lfc_df = pd.DataFrame(list(all_gene_lfc.items()), columns=['symbol', 'mean_lfc'])
lfc_df.to_csv(os.path.join(OUT, 'final_gene_lfc_lookup.csv'), index=False)

# Same treatment for the |LFC| >= 0 consensus.
print("\n=== Fixing LFC0 consensus ===")
master_lfc0 = pd.read_csv(os.path.join(OUT, 'final_master_deg_LFC0.csv'))
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

all_gene_lfc0 = master_lfc0_valid.groupby('symbol')['log2FC'].mean().to_dict()
lfc0_df = pd.DataFrame(list(all_gene_lfc0.items()), columns=['symbol', 'mean_lfc'])
lfc0_df.to_csv(os.path.join(OUT, 'final_gene_lfc_lookup_LFC0.csv'), index=False)

print("\nDone! Ready to rebuild cnets with corrected LFC coloring.")
