"""Montage of one concept network per group, using each group's own LFC0.2 DEGs with no gene-list restriction.

work/group_genes.json and the per-group Enrichr results in work/group_enr/ -> work/montage/Cnet_PerGroup_AllGroups_LFC0.2.png and .pdf.
"""
import os as _os
BASE = _os.environ.get("GLP1R_BASE")
if not BASE:
    raise SystemExit(
        "Set GLP1R_BASE to the directory holding the analysis data tree, e.g.\n"
        "  export GLP1R_BASE=/path/to/workspace"
    )

import os, re, csv, math, json, time, textwrap
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import sys
sys.path.insert(0,BASE + "/mnt/Bulk RNA sequencing/Finalized Bioinformatics Workflow/scripts")
from cnet_style import LIB_COLORS, LIB_LABELS, VMAX
DATA_LOCAL=BASE + "/work"
ENRDIR=BASE + "/work/group_enr"
OUT_LOCAL=BASE + "/work/montage"; os.makedirs(OUT_LOCAL,exist_ok=True)
MAXT=6; MAXG=8

groups=json.load(open(f"{DATA_LOCAL}/group_genes.json"))
import pandas as pd
cat=pd.read_csv(f"{DATA_LOCAL}/data_local/group_catalog.csv") if os.path.exists(f"{DATA_LOCAL}/data_local/group_catalog.csv") else None
order=sorted(groups.keys())
# The catalog also carries groups with no DEGs, which get an empty panel.
all_groups=sorted(cat["group"].tolist()) if cat is not None else order

def clean(t):
    """Strip GO and Reactome accessions from a term."""
    t=re.sub(r'\s*\(GO:\d+\)','',str(t)); return re.sub(r'\s*R-HSA-\d+','',t).strip()
def load_enr(g):
    """Return one group's enrichment rows with adjusted p-value < 0.05, or [] if absent."""
    safe="".join(c if c.isalnum() else "_" for c in g)
    p=f"{ENRDIR}/{safe}.csv"
    if not os.path.exists(p): return []
    rows=[r for r in csv.DictReader(open(p)) if r.get("Adjusted P-value") and r["Adjusted P-value"] not in ("","nan")]
    return [r for r in rows if float(r["Adjusted P-value"])<0.05]

def build_graph(genes, lfc, enr):
    """Build one group's gene-term graph by greedy set cover, or None if too small."""
    input_upper={g.upper() for g in genes}
    if len(genes)<3 or not enr: return None
    enr=sorted(enr, key=lambda r:float(r["Adjusted P-value"]))
    covered=set(); selected=[]; remaining=list(range(len(enr)))
    while len(selected)<MAXT and remaining:
        best=max(remaining,key=lambda i:(len(set(enr[i]["Genes"].upper().split(";"))&input_upper-covered),
                  float(enr[i].get("Combined Score",0) or 0)))
        tg=set(enr[best]["Genes"].upper().split(";"))&input_upper
        if not (tg-covered) and len(selected)>=3: break
        selected.append(best); covered|=tg; remaining.remove(best)
    G=nx.Graph()
    for idx in selected:
        r=enr[idx]; t=clean(r["Term"])[:40]
        G.add_node(("T",t),ntype="term",lib=r["Gene_set"])
        tg=set(r["Genes"].upper().split(";"))&input_upper
        for gu in sorted(tg,key=lambda x:-abs(lfc.get(x,0)))[:MAXG]:
            gn=next((s for s in genes if s.upper()==gu),gu)
            if ("G",gn) not in G: G.add_node(("G",gn),ntype="gene",lfc=lfc.get(gu,0))
            G.add_edge(("T",t),("G",gn))
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n)==0])
    genes_n=[n for n,d in G.nodes(data=True) if d["ntype"]=="gene"]
    if len(genes_n)<2: return None
    return G

def draw(ax,G,title):
    """Draw one montage panel, or a placeholder when the group has too few DEGs."""
    ax.set_title(title,fontsize=8,fontweight="bold")
    ax.axis("off")
    if G is None:
        ax.text(0.5,0.5,"too few DEGs",ha="center",va="center",fontsize=8,color="#999999",transform=ax.transAxes); return
    pos=nx.spring_layout(G,k=1.0,iterations=80,seed=42)
    term_n=[n for n,d in G.nodes(data=True) if d["ntype"]=="term"]
    gene_n=[n for n,d in G.nodes(data=True) if d["ntype"]=="gene"]
    for e in G.edges(): ax.plot([pos[e[0]][0],pos[e[1]][0]],[pos[e[0]][1],pos[e[1]][1]],color="#888",lw=0.8,alpha=0.5,zorder=1)
    gx=np.array([pos[g] for g in gene_n])
    ax.scatter(gx[:,0],gx[:,1],c=[np.clip(G.nodes[g]["lfc"],-VMAX,VMAX) for g in gene_n],
               cmap="RdBu_r",vmin=-VMAX,vmax=VMAX,s=90,zorder=3,edgecolors="white",linewidths=0.4)
    for g in gene_n:
        ax.text(pos[g][0],pos[g][1]-0.04,g[1],fontsize=5,ha="center",va="top",color="white",zorder=5,
                bbox=dict(boxstyle="round,pad=0.1",facecolor="#333",alpha=0.85,edgecolor="none"))
    for t in term_n:
        col=LIB_COLORS.get(G.nodes[t]["lib"],"#888"); xy=pos[t]
        ax.scatter(xy[0],xy[1],s=120,c=col,marker="s",zorder=4,edgecolors="white",linewidths=0.8)
        ax.text(xy[0],xy[1]+0.05,"\n".join(textwrap.wrap(t[1],18)),fontsize=5,ha="center",va="bottom",
                color="white",fontweight="bold",zorder=5,bbox=dict(boxstyle="round,pad=0.15",facecolor=col,alpha=0.9,edgecolor="none"))

N=len(all_groups); ncols=5; nrows=math.ceil(N/ncols)
fig,axes=plt.subplots(nrows,ncols,figsize=(ncols*5.0,nrows*4.4))
fig.patch.set_alpha(0)
axes=axes.flatten()
import matplotlib as mpl
for i,g in enumerate(all_groups):
    d=groups.get(g); 
    if d is None: draw(axes[i],None,f"{g}\n(0 DEGs)"); continue
    G=build_graph(d["genes"],{k.upper():v for k,v in d["lfc"].items()},load_enr(g))
    draw(axes[i],G,f"{g}\n({len(d['genes'])} DEGs, LFC0.2)")
for j in range(N,len(axes)): axes[j].axis("off")
libs=[lib for lib in LIB_COLORS if lib in {"GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023","KEGG_2021_Human","Reactome_2022"}]
handles=[Line2D([0],[0],marker="s",color="w",markerfacecolor=LIB_COLORS[l],markersize=10,label=LIB_LABELS[l]) for l in libs]
fig.legend(handles=handles,loc="lower center",ncol=5,fontsize=11,frameon=False,bbox_to_anchor=(0.5,0.005))
sm=mpl.cm.ScalarMappable(cmap="RdBu_r",norm=mpl.colors.Normalize(-VMAX,VMAX)); sm.set_array([])
cax=fig.add_axes([0.92,0.4,0.012,0.2]); fig.colorbar(sm,cax=cax,label="mean log2FC")
fig.suptitle("Per-group concept networks — each group's LFC0.2 DEGs (no gene restriction)",fontsize=16,fontweight="bold",y=0.997)
fig.tight_layout(rect=[0,0.02,0.91,0.99])
for ext,dpi in [("png",150),("pdf",None)]:
    fig.savefig(f"{OUT_LOCAL}/Cnet_PerGroup_AllGroups_LFC0.2.{ext}",dpi=dpi,bbox_inches="tight",transparent=True)
plt.close(fig)
print("montage saved", N, "panels")
