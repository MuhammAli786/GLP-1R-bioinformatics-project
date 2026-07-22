import os, json, time, sys
import gseapy as gp
DB=['GO_Biological_Process_2023','GO_Molecular_Function_2023','GO_Cellular_Component_2023','KEGG_2021_Human','Reactome_2022']
groups=json.load(open("group_genes.json"))
items=sorted(groups.items(), key=lambda kv: len(kv[1]["genes"]))  # small first (fast wins)
deadline=time.time()+30
done=skip=0
for g,d in items:
    safe="".join(ch if ch.isalnum() else "_" for ch in g)
    out=f"group_enr/{safe}.csv"
    if os.path.exists(out): continue
    genes=d["genes"]
    if len(genes)<3:
        open(out,"w").write("Gene_set,Term,Overlap,P-value,Adjusted P-value,Odds Ratio,Combined Score,Genes\n")
        skip+=1; continue
    if time.time()>deadline: break
    try:
        res=gp.enrichr(gene_list=genes, gene_sets=DB, organism='mouse', no_plot=True).results
        res.to_csv(out, index=False); done+=1
        print(f"  {g}: {len(genes)} genes -> {len(res)} terms")
    except Exception as e:
        print(f"  {g}: ERR {str(e)[:80]}"); time.sleep(2)
print(f"done={done} skip={skip} remaining={sum(1 for g,_ in items if not os.path.exists('group_enr/'+''.join(c if c.isalnum() else '_' for c in g)+'.csv'))}")
