import numpy as np, os
CO="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
tags=["1901","1902","1903","1904","1905","1906"]
def load(yr):
    m=f"{CO}/submesh_{yr}-12-31T00:00:00"
    return (np.loadtxt(f"{m}/map_nod.out").astype(int),
            np.loadtxt(f"{m}/elem2d.out",skiprows=1).astype(int)-1,
            np.loadtxt(f"{m}/cavity_nlvls.out").astype(int),
            np.loadtxt(f"{m}/cavity_elvls.out").astype(int),
            np.loadtxt(f"{m}/nod2d.out",skiprows=1))
# For each mesh: for every OPEN node (cavn==1), count surrounding elements that are
# open at the surface (cave==1). A "pinch" = open node with very few open elements
# while it HAS many total elements (surrounded by cavity).
def pinch_stats(yr):
    mp,el,cavn,cave,nod=load(yr)
    N=cavn.size
    tot=np.zeros(N,int); opn=np.zeros(N,int)
    for e in range(el.shape[0]):
        c1=(cave[e]==1)
        for k in el[e]:
            tot[k]+=1
            if c1: opn[k]+=1
    openmask=(cavn==1)
    # pinched open nodes: open node, has >=4 elems, but <=2 open-surface elems
    pinched=np.where(openmask & (tot>=4) & (opn<=2))[0]
    return mp,pinched,tot,opn,nod,cavn
print(f"{'mesh':>6} {'#open nodes':>11} {'#pinched(open,tot>=4,opn<=2)':>28} {'worst opn/tot':>14}")
crashid=19690
for yr in tags:
    mp,pinched,tot,opn,nod,cavn=pinch_stats(yr)
    nopen=int(np.sum(cavn==1))
    # worst pinch severity among open nodes with tot>=4
    om=(cavn==1)&(tot>=4)
    sev=opn[om]
    worst=sev.min() if sev.size else -1
    # crash node status
    idx=np.where(mp==crashid)[0]
    cs=""
    if idx.size:
        i=idx[0]; cs=f"  crash19690: open?={cavn[i]==1} tot_elems={tot[i]} open_elems={opn[i]}"
    print(f"{yr:>6} {nopen:11d} {pinched.size:28d} {worst:14d}{cs}")
