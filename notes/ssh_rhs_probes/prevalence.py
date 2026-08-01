import numpy as np, os
CO="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
tags=["1901","1902","1903","1904","1905","1906"]
crashid=19690
def load(yr):
    m=f"{CO}/submesh_{yr}-12-31T00:00:00"
    return (np.loadtxt(f"{m}/map_nod.out").astype(int),
            np.loadtxt(f"{m}/elem2d.out",skiprows=1).astype(int)-1,
            np.loadtxt(f"{m}/cavity_nlvls.out").astype(int),
            np.loadtxt(f"{m}/cavity_elvls.out").astype(int),
            np.loadtxt(f"{m}/nlvls.out").astype(int))
def perlevel_all(mp,el,cavn,cave,nln,maxlev=10):
    # returns dict maxid -> (ul, array of wet-element-count per level 1..maxlev), only open nodes ul==1
    N=cavn.size
    # build node->elem adjacency counts per level
    cnt=np.zeros((N,maxlev),int)
    for e in range(el.shape[0]):
        ce=cave[e]
        lo=ce-1
        if lo<0: lo=0
        for k in el[e]:
            # element wet at level lvl (1-based) if cave<=lvl
            for lvl in range(max(ce,1),maxlev+1):
                cnt[k,lvl-1]+=1
    return cnt
data={yr:load(yr) for yr in tags}
cnts={}
for yr in tags:
    mp,el,cavn,cave,nln=data[yr]
    cnts[yr]=(mp,cavn,perlevel_all(mp,el,cavn,cave,nln))
print("Per cycle-transition: open nodes that GAINED newly-trapped upper levels (conn >2 -> <=2)")
print(f"{'transition':>12} {'#open both':>10} {'#nodes+1lvl':>12} {'#nodes>=3lvl':>13} {'#nodes>=4lvl':>13} {'maxlvls trapped':>15} {'crash node trapd':>16}")
for a,b in zip(tags[:-1],tags[1:]):
    mpa,cavna,ca=cnts[a]; mpb,cavnb,cb=cnts[b]
    da={mid:i for i,mid in enumerate(mpa)}
    n1=n3=n4=0; maxtr=0; crashtr=None
    for ib,mid in enumerate(mpb):
        if cavnb[ib]!=1: continue        # open in new
        if mid not in da: continue
        ia=da[mid]
        if cavna[ia]!=1: continue        # open in old too
        # levels newly trapped: old conn>2, new conn<=2
        newtrap=int(np.sum((ca[ia]>2)&(cb[ib]<=2)))
        if newtrap>=1: n1+=1
        if newtrap>=3: n3+=1
        if newtrap>=4: n4+=1
        maxtr=max(maxtr,newtrap)
        if mid==crashid: crashtr=newtrap
    print(f"{a+'->'+b:>12} {'':>10} {n1:12d} {n3:13d} {n4:13d} {maxtr:15d} {str(crashtr):>16}")
