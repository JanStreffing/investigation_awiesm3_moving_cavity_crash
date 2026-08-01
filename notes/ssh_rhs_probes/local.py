import numpy as np, os, netCDF4 as nc
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
CO="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
R="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work/fesom.1905.oce.restart"
LON0,LAT0=-100.382545,-72.8639145
def load(yr):
    m=f"{CO}/submesh_{yr}-12-31T00:00:00"
    return (np.loadtxt(os.path.join(m,"map_nod.out")).astype(int),
            np.loadtxt(os.path.join(m,"nod2d.out"),skiprows=1),
            np.loadtxt(os.path.join(m,"elem2d.out"),skiprows=1).astype(int)-1,
            np.loadtxt(os.path.join(m,"cavity_nlvls.out")),
            np.loadtxt(os.path.join(m,"cavity_depth@node.out")))
mp5,nod5,el5,nl5,cd5=load("1905")
mp6,nod6,el6,nl6,cd6=load("1906")
set5=set(mp5.tolist())
ic=int(np.argmin((nod6[:,1]-LON0)**2+(nod6[:,2]-LAT0)**2))
# nodes within ~0.8 deg lon / 0.35 lat of crash
box=np.where((np.abs(nod6[:,1]-LON0)<0.9)&(np.abs(nod6[:,2]-LAT0)<0.35))[0]
print(f"crash node submesh-idx={ic+1} maxid={mp6[ic]}  box nodes={box.size}")
print("\nbox nodes: maxid  lon      lat     nlvls  draft   NEW?(absent in1905)  prevDraft prevNlvls")
d5=dict(zip(mp5,cd5)); n5=dict(zip(mp5,nl5))
for n in box[np.argsort(nod6[box,0])] if False else box:
    mid=mp6[n]
    new = mid not in set5
    pd = d5.get(mid,np.nan); pn=n5.get(mid,np.nan)
    tag="  <== crash" if n==ic else ("  NEW" if new else "")
    print(f"  {mid:7d} {nod6[n,1]:8.3f} {nod6[n,2]:8.3f} {nl6[n]:5.0f} {cd6[n]:7.1f}   {'NEW' if new else '   '}    {pd:8.1f} {pn if pn==pn else -1:6.0f}{tag}")
# element quality near crash: area + min interior angle
R_E=6371000.0; rad=np.pi/180
def geom(nod,elem,ndsel):
    x=R_E*rad*(nod[:,1]-LON0)*np.cos(LAT0*rad); y=R_E*rad*(nod[:,2]-LAT0)
    emask=np.where(np.any(np.isin(elem,ndsel),axis=1))[0]
    out=[]
    for e in emask:
        p=elem[e]; X=x[p]; Y=y[p]
        a=np.hypot(X[1]-X[0],Y[1]-Y[0]); b=np.hypot(X[2]-X[1],Y[2]-Y[1]); c=np.hypot(X[0]-X[2],Y[0]-Y[2])
        ar=0.5*abs((X[1]-X[0])*(Y[2]-Y[0])-(X[2]-X[0])*(Y[1]-Y[0]))
        # min angle via law of cosines
        import math
        def ang(o,u,v):
            import numpy as np
            d1=np.array([X[u]-X[o],Y[u]-Y[o]]); d2=np.array([X[v]-X[o],Y[v]-Y[o]])
            cs=np.dot(d1,d2)/(np.linalg.norm(d1)*np.linalg.norm(d2)+1e-30)
            return math.degrees(math.acos(max(-1,min(1,cs))))
        mn=min(ang(0,1,2),ang(1,0,2),ang(2,0,1))
        out.append((e,ar,mn))
    return out
q=geom(nod6,el6,box)
ars=np.array([o[1] for o in q]); mns=np.array([o[2] for o in q])
print(f"\ncrash-region elements: n={len(q)}  area min/med={ars.min():.3e}/{np.median(ars):.3e}  min-angle min={mns.min():.2f} deg (median {np.median(mns):.1f})")
print("worst (smallest min-angle) elems:")
for e,ar,mn in sorted(q,key=lambda o:o[2])[:6]:
    nds=el6[e]; print(f"  elem {e+1} area={ar:.3e} minang={mn:.2f}  nodes(maxid)={[mp6[k] for k in nds]} lon/lat={[(round(nod6[k,1],2),round(nod6[k,2],2)) for k in nds]}")
