import numpy as np, os
CO="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
tags=["1901","1902","1903","1904","1905","1906"]
def load(yr):
    m=f"{CO}/submesh_{yr}-12-31T00:00:00"
    mp=np.loadtxt(os.path.join(m,"map_nod.out")).astype(int)      # maxmesh id per submesh node
    cd=np.loadtxt(os.path.join(m,"cavity_depth@node.out"))
    nl=np.loadtxt(os.path.join(m,"cavity_nlvls.out"))
    nod=np.loadtxt(os.path.join(m,"nod2d.out"),skiprows=1)
    return mp,cd,nl,nod
data={yr:load(yr) for yr in tags}
# crash location maxid
mp6,cd6,nl6,nod6=data["1906"]
LON0,LAT0=-100.382545,-72.8639145
ic=int(np.argmin((nod6[:,1]-LON0)**2+(nod6[:,2]-LAT0)**2))
crash_maxid=mp6[ic]
print(f"crash maxmesh id={crash_maxid} at lon {nod6[ic,1]:.3f} lat {nod6[ic,2]:.3f}\n")
print(f"{'transition':>14} {'common':>8} {'max|dDraft|':>11} {'p99':>7} {'p99.9':>7} {'n|d|>40':>8} {'n|d|>70':>8} {'nlvls chg':>9} {'Amundsen max|d|':>15}")
for a,b in zip(tags[:-1],tags[1:]):
    mpa,cda,nla,noda=data[a]; mpb,cdb,nlb,nodb=data[b]
    da=dict(zip(mpa,cda)); db=dict(zip(mpb,cdb))
    na=dict(zip(mpa,nla)); nb=dict(zip(mpb,nlb))
    lb=dict(zip(mpb,zip(nodb[:,1],nodb[:,2])))
    common=[k for k in db if k in da]
    dd=np.array([db[k]-da[k] for k in common])
    dnl=np.array([nb[k]-na[k] for k in common])
    # Amundsen box (near crash)
    amu=[abs(db[k]-da[k]) for k in common if abs(lb[k][0]-LON0)<3 and abs(lb[k][1]-LAT0)<2]
    amax=max(amu) if amu else 0
    print(f"{a+'->'+b:>14} {len(common):8d} {np.max(np.abs(dd)):11.2f} {np.percentile(np.abs(dd),99):7.2f} {np.percentile(np.abs(dd),99.9):7.2f} {int(np.sum(np.abs(dd)>40)):8d} {int(np.sum(np.abs(dd)>70)):8d} {int(np.sum(dnl!=0)):9d} {amax:15.2f}")
# crash node history
print("\ncrash-node draft & nlvls history (by maxid):")
for yr in tags:
    mp,cd,nl,nod=data[yr]
    idx=np.where(mp==crash_maxid)[0]
    if idx.size: print(f"  {yr}: draft={cd[idx[0]]:8.2f}  nlvls={nl[idx[0]]:.0f}")
    else: print(f"  {yr}: (not in mesh)")
