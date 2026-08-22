# Command-line driver that generates every cnet plot (consensus, pathway, and
# specific-pathway sets) at the LFC0/LFC02/LFC05/LFC1 thresholds.

import os, sys, csv, argparse

# Make the sibling cnet_* modules importable when run from the Scripts folder.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from cnet_style import THRESHOLD_LABELS
from cnet_gene_lists import (
    BBB_GENES, JAK_STAT3_GENES, PI3K_AKT_GENES,
    ION_CHANNEL_BASE_GENES, KEYWORD_FILTERS,
)
from cnet_engine import build_cnet, load_consensus, load_enrichment, load_lfc

def get_ion_genes(workspace_root):
    """Merge the base ion channel list with Data_Files/Ion_Channel_Consensus_Genes.csv."""
    extra = []
    csv_path = os.path.join(workspace_root, 'Data_Files',
                            'Ion_Channel_Consensus_Genes.csv')
    if os.path.exists(csv_path):
        try:
            rows = list(csv.DictReader(open(csv_path)))
            extra = [r['gene_symbol'] for r in rows]
            print(f"  Loaded {len(extra)} genes from Ion_Channel_Consensus_Genes.csv")
        except Exception as e:
            print(f"  Warning: could not read {csv_path}: {e}")
    return list(set(ION_CHANNEL_BASE_GENES + extra))



THRESHOLDS = ['LFC0', 'LFC02', 'LFC05', 'LFC1']

def generate_consensus_cnets(data_dir, out_base, lfc_map):
    """Build consensus-gene cnet plots at every threshold."""
    out_dir = os.path.join(out_base, 'Plots', 'Cnet_Plots_Consensus_Based')
    print("\n" + "=" * 60)
    print("CONSENSUS-BASED CNET PLOTS")
    print("=" * 60)

    for thr in THRESHOLDS:
        tl = THRESHOLD_LABELS[thr]
        cons_map = load_consensus(data_dir, thr)
        enr_rows = load_enrichment(data_dir, thr)
        print(f"\n--- {thr} ({len(cons_map)} consensus genes, "
              f"{len(enr_rows)} sig terms) ---")

        build_cnet(
            ref_genes=list(cons_map.values()),
            term_filter=None,
            title=f'Consensus Gene Cnet ({tl})\n{len(cons_map)} genes',
            fname=f'Cnet_Consensus_{thr}',
            out_dir=out_dir,
            lfc_map=lfc_map,
            cons_map=cons_map,
            enr_rows=enr_rows,
        )


def generate_pathway_cnets(data_dir, out_base, lfc_map):
    """Build pathway cnet plots (BBB/MMP, JAK-STAT3 inflammatory, PI3K/Akt survival)."""
    out_dir = os.path.join(out_base, 'Plots', 'Cnet_Plots_Pathway_Based')
    print("\n" + "=" * 60)
    print("PATHWAY-BASED CNET PLOTS")
    print("=" * 60)

    pathways = [
        (BBB_GENES,       KEYWORD_FILTERS['BBB'],
         'MMP-2/9 & Blood-Brain Barrier',     'Cnet_BBB_MMP'),
        (JAK_STAT3_GENES, KEYWORD_FILTERS['INFLAMMATORY'],
         'JAK-STAT3 / Inflammatory Pathway',   'Cnet_JAKSTAT3_Inflammatory'),
        (PI3K_AKT_GENES,  KEYWORD_FILTERS['SURVIVAL'],
         'PI3K/Akt Pro-Survival Pathway',      'Cnet_AktPI3K_ProSurvival'),
    ]

    for thr in THRESHOLDS:
        tl = THRESHOLD_LABELS[thr]
        cons_map = load_consensus(data_dir, thr)
        enr_rows = load_enrichment(data_dir, thr)
        print(f"\n--- {thr} ---")

        for genes, kw_filter, label, prefix in pathways:
            build_cnet(
                ref_genes=genes,
                term_filter=kw_filter,
                title=f'{label} – Concept Network ({tl})',
                fname=f'{prefix}_{thr}',
                out_dir=out_dir,
                lfc_map=lfc_map,
                cons_map=cons_map,
                enr_rows=enr_rows,
            )


def generate_specific_cnets(data_dir, out_base, lfc_map, ion_genes):
    """Build the specific-pathway cnet plots (BBB, inflammatory, survival, ion channel)."""
    out_dir = os.path.join(out_base, 'Specific Pathway Cnets')
    print("\n" + "=" * 60)
    print("SPECIFIC PATHWAY CNET PLOTS")
    print("=" * 60)

    pathways = [
        (BBB_GENES,       KEYWORD_FILTERS['BBB'],
         'BBB Pathway',           'Cnet_BBB'),
        (JAK_STAT3_GENES, KEYWORD_FILTERS['INFLAMMATORY'],
         'JAK-STAT3 Inflammatory', 'Cnet_Inflammatory'),
        (PI3K_AKT_GENES,  KEYWORD_FILTERS['SURVIVAL'],
         'PI3K/Akt Survival',      'Cnet_Survival'),
        (ion_genes,        KEYWORD_FILTERS['ION_CHANNEL'],
         'Ion Channel',            'Cnet_IonChannel'),
    ]

    for thr in THRESHOLDS:
        tl = THRESHOLD_LABELS[thr]
        cons_map = load_consensus(data_dir, thr)
        enr_rows = load_enrichment(data_dir, thr)
        print(f"\n--- {thr} ---")

        for genes, kw_filter, label, prefix in pathways:
            build_cnet(
                ref_genes=genes,
                term_filter=kw_filter,
                title=f'{label} – Concept Network ({tl})',
                fname=f'{prefix}_{thr}',
                out_dir=out_dir,
                lfc_map=lfc_map,
                cons_map=cons_map,
                enr_rows=enr_rows,
            )


def main():
    parser = argparse.ArgumentParser(
        description='Generate all Cnet plots for the GLP-1R CNS meta-analysis')
    parser.add_argument(
        '--data-dir', default=None,
        help='Directory containing final_consensus_*.csv, '
             'final_enrichment_*.csv, and final_gene_lfc_COMPREHENSIVE.csv')
    parser.add_argument(
        '--out-dir', default=None,
        help='Root output directory (Final analysis folder)')
    parser.add_argument(
        '--workspace', default=None,
        help='Workspace root (parent of Data_Files/). '
             'Used to find Ion_Channel_Consensus_Genes.csv')
    args = parser.parse_args()

    # Fall back to the usual data locations relative to the Scripts folder.
    if args.data_dir is None:
        candidates = [
            os.path.join(SCRIPT_DIR, '..', '..', 'outputs'),
            os.path.join(SCRIPT_DIR, '..', 'Data_Files'),
        ]
        for c in candidates:
            if os.path.exists(os.path.join(c, 'final_gene_lfc_COMPREHENSIVE.csv')):
                args.data_dir = os.path.abspath(c)
                break
        if args.data_dir is None:
            print("ERROR: Could not find data directory. Use --data-dir.")
            sys.exit(1)

    if args.out_dir is None:
        args.out_dir = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))

    if args.workspace is None:
        args.workspace = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))

    print(f"Data dir:   {args.data_dir}")
    print(f"Output dir: {args.out_dir}")
    print(f"Workspace:  {args.workspace}")

    print("\nLoading LFC map...")
    lfc_map = load_lfc(args.data_dir)
    print(f"  {len(lfc_map)} genes loaded")

    ion_genes = get_ion_genes(args.workspace)
    print(f"  {len(ion_genes)} ion channel genes total")

    generate_consensus_cnets(args.data_dir, args.out_dir, lfc_map)
    generate_pathway_cnets(args.data_dir, args.out_dir, lfc_map)
    generate_specific_cnets(args.data_dir, args.out_dir, lfc_map, ion_genes)

    print("\n" + "=" * 60)
    print("ALL DONE")
    print("=" * 60)


if __name__ == '__main__':
    main()
