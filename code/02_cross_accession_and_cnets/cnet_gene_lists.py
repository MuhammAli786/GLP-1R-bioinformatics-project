# Curated pathway gene lists and enrichment-term keyword filters used by the cnet scripts.

# Blood-brain barrier (44 genes): tight and adherens junctions, transporters,
# MMPs, basement membrane, endothelial signalling.
BBB_GENES = [
    'Adam10', 'Adam9', 'Adam15', 'Cdh5', 'Ctnnb1', 'Ctnnd1',
    'Vegfa', 'Hif1a', 'Aqp4', 'Slc2a1', 'Lrp1', 'Slc7a5',
    'Slc16a1', 'Tfrc', 'Abcb1b', 'Nid1', 'Lamc1', 'Fn1',
    'Spp1', 'Kdr', 'Cav1', 'Vwf', 'Pecam1', 'Eng', 'Nos3',
    'Mmp2', 'Mmp9', 'Mmp16', 'Mmp14', 'Mmp25',
    'Nfkbia', 'Nfkb2', 'Mapk1', 'Rock2', 'Fyn', 'Rock1', 'Akt1',
    'Cldn5', 'Tjp1', 'Jam2', 'F11r', 'Cldn12', 'Ocln', 'Tjp2',
]

# JAK-STAT3 inflammatory (33 genes): cytokine receptors, JAK/STAT cascade,
# SOCS negative regulators, downstream transcriptional targets.
JAK_STAT3_GENES = [
    'Gfap', 'Il1r1', 'Il6ra', 'Il6st', 'Lifr', 'Osmr',
    'Mapk1', 'Fyn', 'Akt1', 'Jak1', 'Nfkbia', 'Nfkb2',
    'Ikbkg', 'Pias2', 'Pias1', 'Ptpn11', 'Ptprd', 'Ptprt',
    'Egfr', 'Erbb2', 'Fgfr1', 'Socs5', 'Socs3', 'Socs4', 'Socs6',
    'Stat1', 'Cdkn1a', 'Bcl2l1', 'Vegfa', 'Hif1a', 'Mmp2',
    'Mcl1', 'Ccnd1',
]

# PI3K/Akt pro-survival (68 unique genes): PI3K subunits, Akt isoforms, mTOR
# complex, FOXO transcription factors, apoptosis, autophagy, RAS/MAPK crosstalk.
PI3K_AKT_GENES = list(set([
    'Irs1', 'Irs2', 'Sos1', 'Gab1', 'Grb2', 'Pdpk1', 'Akt3', 'Akt1',
    'Bcl2l1', 'Bad', 'Mcl1', 'Bax', 'Xiap', 'Birc2',
    'Ulk2', 'Map1lc3a', 'Becn1', 'Ulk1', 'Creb1', 'Creb5', 'Creb3',
    'Cdkn1a', 'Ccnd1', 'Foxo1', 'Foxo3', 'Foxo4', 'Foxo6',
    'Gsk3a', 'Gsk3b', 'Vegfa', 'Kdr', 'Insr', 'Igf1r',
    'Egfr', 'Erbb2', 'Erbb3', 'Fgfr1', 'Fgfr2', 'Glp1r',
    'Nfkbia', 'Nfkb2', 'Phlpp1', 'Phlpp2',
    'Pik3ca', 'Pik3r1', 'Pik3r3', 'Pik3c2a', 'Pik3cb', 'Pik3cg', 'Pik3c2b',
    'Pten', 'Hras', 'Kras', 'Mapk1', 'Map2k2', 'Nras', 'Raf1', 'Fyn',
    'Mtor', 'Tsc2', 'Rptor', 'Rictor', 'Rheb', 'Tsc1',
    'Eif4g1', 'Eif4ebp1', 'Rps6kb1', 'Rps6kb2',
]))

# Ion channel (66+ genes): voltage-gated K+/Na+/Ca2+ channels, glutamate and GABA
# receptors, TRP, HCN, and chloride channels. This base list is augmented at
# runtime from Ion_Channel_Consensus_Genes.csv.
ION_CHANNEL_BASE_GENES = [
    'Scn2a', 'Scn1a', 'Kcnma1', 'Kcnj11', 'Cacna1c', 'Cacna2d1',
    'Gria1', 'Gria2', 'Grin1', 'Grin2a', 'Gabra1', 'Gabrb3',
    'Kcnb1', 'Kcnk2', 'Kcnk3', 'Hcn1', 'Hcn2', 'Trpv1', 'Trpm8', 'Trpc6',
    'Cacnb4', 'Scn3a', 'Scn9a', 'Kcnj4', 'Kcnj16', 'Clcn5',
    'Grik2', 'Grik1', 'Grm1', 'Grm5',
]


# Regex keyword filters used to select enrichment terms per pathway.
KEYWORD_FILTERS = {
    'BBB': (
        r'(?i)blood.brain|BBB|tight.junction|endotheli|barrier|claudin|occludin|'
        r'cell.junction|adherens|basement.membrane|extracellular.matrix|focal.adhesion|'
        r'vascular|angiogen|matrix.metallo|MMP|collagen|laminin|integrin|cadherin|'
        r'VEGF|Rho.GTPase|actin.cytoskeleton|leukocyte.transendothelial|cell.adhesion|'
        r'blood.vessel|endothelial|cell.migrat|tube.morpho|wound.heal|permeab|caveol|'
        r'efflux|transporter|aquaporin'
    ),
    'INFLAMMATORY': (
        r'(?i)JAK|STAT|interleukin|cytokine|inflam|interferon|NF.?kB|tumor.necrosis|'
        r'TNF|chemokine|toll.like|innate.immune|adaptive.immune|complement|apoptosis|'
        r'programmed.cell.death|caspase|cell.death|MAPK|signal.transduct|'
        r'growth.factor.recep|phosphatase|kinase.activ|receptor.signal|Src.family|immune'
    ),
    'SURVIVAL': (
        r'(?i)PI3K|AKT|mTOR|phosphoinositide|insulin|IGF|growth.factor|FOXO|GSK|'
        r'autophagy|apoptosis|pro.?surviv|MAPK|RAS|RAF|MEK|ERK|PTEN|phosphoryl|'
        r'kinase.signal|RTK|receptor.tyrosine|CREB|translation|EIF4|ribosom|'
        r'cell.cycle|proliferat|survival|ubiquitin|proteolys|longevity|PIP3'
    ),
    'ION_CHANNEL': (
        r'(?i)ion.channel|potassium|sodium|calcium|chloride|GABA|glutamate|'
        r'synap|postsynap|presynap|dendrit|axon|neurotransmit|action.potential|'
        r'membrane.potential|voltage.gated|ligand.gated|ion.transport|channel.activity|'
        r'receptor.channel|excitatory|inhibitory|cholinerg|dopamin|seroton|'
        r'neuropeptide|synaptic.vesicle|synaptic.transmission'
    ),
}
