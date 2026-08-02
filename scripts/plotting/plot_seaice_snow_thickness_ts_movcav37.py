#!/usr/bin/env python3
# RESTORED 2026-07-22 from session transcript (original was in node-local /tmp, lost on node switch).
# Timeseries domain: open ocean, cavity nodes excluded, AND enclosed PISM lagoons excluded BY AREA
# (small open-ocean components disconnected from the world ocean by the ice-shelf/cavity ring; they
# fill with ~10 m sea ice). Lagoons kept visible in the maps.
import netCDF4 as nc, numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R="/work/ab0246/a270092/runtime/awiesm3-develop-is/movcav37"
POOL="/work/ab0246/a270092/input/fesom2/pism_cavity_ini/submesh"
def mesh_for(y): return POOL if y<=1900 else f"{R}/couple/submesh_{y}-12-31T00:00:00"
def load_mesh(s):
    ao=np.loadtxt(f"{s}/nod2d.out",skiprows=1); lon,lat=ao[:,1],ao[:,2]
    ul=np.loadtxt(f"{s}/cavity_nlvls.out").astype(int)   # >1 = cavity (under ice shelf)
    e=np.loadtxt(f"{s}/elem2d.out",skiprows=1).astype(int)-1
    R0=6371e3; _la=np.radians(lat); _lo=np.radians(lon)
    # spherical excess (Girard): plane formula explodes on dateline-crossing elems
    v=np.stack([np.cos(_la)*np.cos(_lo),np.cos(_la)*np.sin(_lo),np.sin(_la)],axis=1)
    va,vb,vc=v[e[:,0]],v[e[:,1]],v[e[:,2]]
    _num=np.abs(np.einsum('ij,ij->i',va,np.cross(vb,vc)))
    _den=1.0+np.einsum('ij,ij->i',va,vb)+np.einsum('ij,ij->i',vb,vc)+np.einsum('ij,ij->i',vc,va)
    ar=2.0*np.arctan2(_num,_den)*R0*R0
    na=np.zeros(len(ao))
    for k in range(3): np.add.at(na,e[:,k],ar/3.0)
    # exclude lagoons BY AREA: flood-fill open-ocean nodes through the mesh and keep only
    # ocean-scale components (>1e6 km^2). Enclosed pockets (Ross embayment, ~5e3 km^2, cut off
    # from the world ocean by the cavity ring) fall out; the global ocean (8.6e8 km^2) is kept.
    oo=(ul<=1)
    pairs=np.vstack([e[:,[0,1]],e[:,[1,2]],e[:,[0,2]]])
    both=oo[pairs[:,0]]&oo[pairs[:,1]]; pp=pairs[both]
    A=coo_matrix((np.ones(len(pp)),(pp[:,0],pp[:,1])),shape=(len(ao),len(ao)))
    ncomp,lab=connected_components(A,directed=False)
    carea=np.bincount(lab[oo],weights=na[oo],minlength=ncomp)   # m^2 per component
    big=carea>1e12                                             # 1e6 km^2 = ocean-scale
    domain=oo&big[lab]
    return lat,na,domain,int((oo&~big[lab]).sum())
def rd(f,v):
    d=nc.Dataset(f); a=np.ma.filled(np.ma.masked_invalid(d.variables[v][:]),np.nan); d.close()
    a[np.abs(a)>1e29]=np.nan; return a
tt=[]; nlag=0; S={'shI':[],'shS':[],'nhI':[],'nhS':[]}
for y in range(1900,1910):
    s=mesh_for(y); lat,na,domain,nlag_y=load_mesh(s); nlag=max(nlag,nlag_y)
    mi=rd(f"{R}/outdata/fesom/m_ice.fesom.{y}.nc",'m_ice')
    ms=rd(f"{R}/outdata/fesom/m_snow.fesom.{y}.nc",'m_snow')
    ai=rd(f"{R}/outdata/fesom/a_ice.fesom.{y}.nc",'a_ice')
    sh=(lat<-45)&domain; nh=(lat>45)&domain
    for m in range(12):
        tt.append(y+(m+0.5)/12.)
        for tag,msk in (('sh',sh),('nh',nh)):
            cov=msk&(ai[m]>0.15)&np.isfinite(mi[m])&np.isfinite(ms[m])
            A=na[cov]; iarea=np.sum(ai[m,cov]*A)
            if iarea>0:
                S[tag+'I'].append(np.sum(mi[m,cov]*A)/iarea)   # vol/ice-area = mean actual thickness
                S[tag+'S'].append(np.sum(ms[m,cov]*A)/iarea)
            else: S[tag+'I'].append(np.nan); S[tag+'S'].append(np.nan)
tt=np.array(tt)
fig,axs=plt.subplots(2,2,figsize=(13,7),facecolor="white",sharex=True)
def band(ax,lo,hi,txt): ax.axhspan(lo,hi,color="#cfe0d5",alpha=0.7,zorder=0); ax.text(1900.1,(lo+hi)/2,txt,fontsize=8,va="center",color="0.3")
cfg=[(axs[0,0],'shI',"#2a5fb5",(0.5,1.0),"obs 0.5–1.0 m (ASPeCt/ICESat)","(a) SH sea-ice thickness (vol/area)","ice thickness [m]"),
     (axs[0,1],'nhI',"#2a5fb5",(1.5,3.0),"obs ~1.5–3 m (Arctic)","(b) NH sea-ice thickness (vol/area)",None),
     (axs[1,0],'shS',"#7a3fb5",(0.1,0.4),"obs 0.1–0.4 m","(c) SH snow depth on sea ice","snow depth [m]"),
     (axs[1,1],'nhS',"#7a3fb5",(0.2,0.4),"obs 0.2–0.4 m (Warren)","(d) NH snow depth on sea ice",None)]
for ax,k,c,(lo,hi),ob,ti,yl in cfg:
    band(ax,lo,hi,ob); ax.plot(tt,S[k],color=c,lw=1.4); ax.set_title(ti,fontsize=10.5)
    if yl: ax.set_ylabel(yl)
    ax.spines[["top","right"]].set_visible(False); ax.set_xlim(1900,1910); ax.grid(axis="y",color="0.9")
for ax in axs[1]: ax.set_xlabel("year")
fig.suptitle("movcav37: sea-ice & snow THICKNESS vs obs — hemispheric mean (volume/ice-area), monthly 1900–1909\n(world-ocean only; enclosed lagoons excluded by connected-area, cavity nodes excluded — shown only in the maps)",fontsize=11)
fig.tight_layout(rect=(0,0,1,0.96))
out="/home/a/a270092/esm_tools/movcav4_crash_plots/movcav37_seaice_snow_thickness_ts.png"; fig.savefig(out,dpi=145); print("wrote",out,"| lagoon nodes excluded:",nlag)
i5=tt>=1905
for nm,k in (("SH ice",'shI'),("NH ice",'nhI'),("SH snow",'shS'),("NH snow",'nhS')):
    a=np.array(S[k])[i5]; print(f"{nm:8s} last-5yr: winter-max={np.nanmax(a):.2f}  summer-min={np.nanmin(a):.2f}  mean={np.nanmean(a):.2f} m")
