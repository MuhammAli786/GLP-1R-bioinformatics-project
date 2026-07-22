import os, re, textwrap
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import sys
sys.path.insert(0,"/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/scripts")
from cnet_style import *
from cnet_gene_lists import PI3K_AKT_GENES, KEYWORD_FILTERS
D="/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/Data"
OUT="/sessions/amazing-zen-bardeen/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/Plots/Cnet plots/Survival/LFC0.2"

cons=pd.read_csv(f"{D}/consensus_LFC02.csv")
cons_map={r.gene_symbol.upper():r.gene_symbol for r in cons.itertuples()}
comp=pd.read_csv(f"{D}/gene_lfc_comprehensive_LFC02.csv")
lfc_map={r.gene_symbol.upper():float(r.mean_log2FC) for r in comp.itertuples()}
input_upper={g.upper() for g in PI3K_AKT_GENES if g.upper() in cons_map}
enr=[r for _,r in pd.read_csv(f"{D}/enrichment_LFC02.csv").iterrows() if r['Adjusted P-value']<0.05]
pat=re.compile(KEYWORD_FILTERS["SURVIVAL"])
ft=[r for r in enr if pat.search(str(r['Term'])) or len(set(str(r['Genes']).upper().split(';'))&input_upper)>=2]
def clean(t): t=re.sub(r"\s*\(GO:\d+\)","",str(t)); return re.sub(r"\s*R-HSA-\d+","",t)[:52]
def wrap(t): return "\n".join(textwrap.wrap(t, TERM_WRAP_WIDTH))

G=nx.Graph()
for r in ft:                       # NO LIMIT: include every filtered term
    tg=set(str(r['Genes']).upper().split(';'))&input_upper
    if not tg: continue
    t=clean(r['Term'])
    G.add_node(("T",t),ntype="term",lib=r['Gene_set'])
    for gu in sorted(tg,key=lambda g:-abs(lfc_map.get(g,0)))[:MAX_GENES_PER_TERM]:
        gn=cons_map.get(gu,gu)
        if ("G",gn) not in G: G.add_node(("G",gn),ntype="gene",lfc=lfc_map.get(gu,0))
        G.add_edge(("T",t),("G",gn))
G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n)==0])
term_nodes=[n for n,d in G.nodes(data=True) if d['ntype']=="term"]
gene_nodes=[n for n,d in G.nodes(data=True) if d['ntype']=="gene"]
print(f"NO-LIMIT Survival: {len(term_nodes)} terms, {len(gene_nodes)} genes, {G.number_of_edges()} edges")
pos=nx.spring_layout(G,k=K_SPRING,iterations=ITERATIONS,seed=SEED)
fig,ax=plt.subplots(figsize=(46,40)); fig.patch.set_alpha(0); ax.set_facecolor("none")
nx.draw_networkx_edges(G,pos,ax=ax,alpha=EDGE_ALPHA,width=EDGE_WIDTH,edge_color=EDGE_COLOR)
gsz=[max(GENE_MIN_SIZE,G.degree(g)*GENE_DEGREE_SCALE) for g in gene_nodes]
gxy=np.array([pos[g] for g in gene_nodes])
sc=ax.scatter(gxy[:,0],gxy[:,1],c=[np.clip(G.nodes[g]['lfc'],-VMAX,VMAX) for g in gene_nodes],
              cmap="RdBu_r",vmin=-VMAX,vmax=VMAX,s=gsz,zorder=3,edgecolors=GENE_EDGE_COLOR,linewidths=GENE_EDGE_WIDTH)
cb=plt.colorbar(sc,ax=ax,shrink=0.3,pad=0.02); cb.set_label("mean log2FC",fontsize=24)
for t in term_nodes:
    col=LIB_COLORS.get(G.nodes[t]['lib'],"#888"); xy=pos[t]
    ax.scatter(xy[0],xy[1],s=TERM_SQUARE_SIZE,c=col,zorder=4,edgecolors=TERM_EDGE_COLOR,linewidths=TERM_EDGE_WIDTH,marker="s")
    ax.text(xy[0],xy[1]+TERM_LABEL_OFFSET,wrap(t[1]),fontsize=11,fontweight="bold",color="white",ha="center",va="bottom",zorder=5,
            bbox=dict(boxstyle="round,pad=0.2",facecolor=col,alpha=TERM_BOX_ALPHA,edgecolor="none"))
for g in gene_nodes:
    xy=pos[g]; ax.text(xy[0],xy[1]-GENE_LABEL_OFFSET,g[1],fontsize=16,fontweight="bold",color="white",ha="center",va="top",zorder=6,
            bbox=dict(boxstyle="round,pad=0.2",facecolor=GENE_BOX_COLOR,alpha=GENE_BOX_ALPHA,edgecolor="none"))
libs={G.nodes[t]['lib'] for t in term_nodes}
leg=[Line2D([0],[0],marker="s",color="w",markerfacecolor=c,markersize=14,label=LIB_LABELS.get(l,l)) for l,c in LIB_COLORS.items() if l in libs]
ax.legend(handles=leg,loc="lower left",fontsize=22,framealpha=0.8,facecolor="white")
ax.set_title("PI3K/Akt Pro-Survival Cnet - ALL filtered terms (no limit)  (LFC 0.2)",fontsize=30,fontweight="bold",pad=20); ax.axis("off")
for fmt in ("PNG","PDF"):
    d=os.path.join(OUT,fmt); os.makedirs(d,exist_ok=True)
    fig.savefig(os.path.join(d,f"Cnet_Survival_NoLimit_LFC0.2.{fmt.lower()}"),dpi=120,bbox_inches="tight",transparent=True)
plt.close(fig); print("saved")
