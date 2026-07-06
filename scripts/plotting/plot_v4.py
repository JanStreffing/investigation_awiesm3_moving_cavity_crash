import numpy as np, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation
from scipy.spatial import cKDTree
C="/work/ab0246/a270092/pism_repro/experiments/movcav8/couple"; R=f"{C}/restart_remapped_v4"
OUT="/home/a/a270092/esm_tools/movcav4_crash_plots/initstate"
LON0,LAT0=-62.96,-66.81; KLEV=11; DLON,DLAT=3.0,1.3
ao=np.loadtxt(f"{C}/latest_submesh/nod2d.out",skiprows=1); lon,lat=ao[:,1],ao[:,2]; nn=len(lon)
uln=np.loadtxt(f"{C}/latest_submesh/cavity_nlvls.out").astype(int); nln=np.loadtxt(f"{C}/latest_submesh/nlvls.out").astype(int)
oo=np.loadtxt(f"{C}/mother_as_old/nod2d.out",skiprows=1); ulo=np.loadtxt(f"{C}/mother_as_old/cavity_nlvls.out").astype(int)
io=cKDTree(np.c_[oo[:,1],oo[:,2]]).query(np.c_[lon,lat])[1]; klass=np.where(uln<ulo[io],1,np.where(uln>ulo[io],2,0))
elem=np.loadtxt(f"{C}/latest_submesh/elem2d.out",skiprows=1).astype(int)-1
with open(f"{C}/latest_submesh/aux3d.out") as f:
    nlv=int(f.readline()); zbar=np.array([float(f.readline()) for _ in range(nlv)])
zmid=0.5*(zbar[:-1]+zbar[1:])
def n3(fn,v):
    a=np.squeeze(np.asarray(nc.Dataset(f"{R}/{fn}")[v][:])); return a.T if (a.ndim==2 and a.shape[0]==nn) else a
def e3(fn,v):
    a=np.squeeze(np.asarray(nc.Dataset(f"{R}/{fn}")[v][:])); return a.T if (a.ndim==2 and a.shape[0]==len(elem)) else a
T=n3("temp.nc","temp"); S=n3("salt.nc","salt"); W=n3("w.nc","w"); H=n3("hnode.nc","hnode")
U=e3("u.nc","u"); V=e3("v.nc","v"); SPDe=np.sqrt(U**2+V**2)
cnt=np.bincount(elem.ravel(),minlength=nn); SPDn=np.zeros((SPDe.shape[0],nn))
for k in range(SPDe.shape[0]):
    acc=np.zeros(nn); np.add.at(acc,elem.ravel(),np.repeat(SPDe[k],3)); SPDn[k]=acc/np.maximum(cnt,1)
ev=lon[elem]; lonspan=ev.max(1)-ev.min(1); NOWRAP=lonspan<180        # <-- the fix
ecx=lon[elem].mean(1); ecy=lat[elem].mean(1)
def wet(k): return (uln<=k+1)&(nln-1>=k+1)
def verts(sel): e=elem[sel]; return np.stack([np.column_stack([lon[e[:,i]],lat[e[:,i]]]) for i in range(3)],axis=1)
# cavity-outline triangulation (mask wrapping tris)
trg=Triangulation(lon,lat,elem); trg.set_mask(~NOWRAP)
def cavity_outline(ax):
    try: ax.tricontour(trg,uln.astype(float),levels=[1.5],colors='k',linewidths=1.6)
    except Exception as e: print("outline:",e)
# ===== plan view =====
k=KLEV-1
sel=np.where(NOWRAP & (ecx>LON0-DLON)&(ecx<LON0+DLON)&(ecy>LAT0-DLAT)&(ecy<LAT0+DLAT) & wet(k)[elem].all(axis=1))[0]
Vv=verts(sel)
def tn(f): return f[k][elem[sel]].mean(axis=1)
fig,axs=plt.subplots(2,3,figsize=(17,9)); axs=axs.ravel()
for ax,(ttl,arr,cm,opt) in zip(axs,[("temp (C)",tn(T),'RdYlBu_r',None),("salt (psu)",tn(S),'viridis',None),
        ("speed |u,v| (m/s)",SPDe[k][sel],'magma',None),("w (m/s)",tn(W),'RdBu_r','sym'),("hnode (m)",tn(H),'cividis',None)]):
    pc=PolyCollection(Vv,array=arr,cmap=cm,edgecolors='face',linewidths=0.1)
    if opt=='sym':
        a=np.nanmax(np.abs(arr)) or 1e-6; pc.set_clim(-a,a)
    ax.add_collection(pc); plt.colorbar(pc,ax=ax,shrink=.85); cavity_outline(ax)
    ax.set_title(ttl); ax.plot(LON0,LAT0,'k*',ms=15,mec='w')
    ax.set_xlim(LON0-DLON,LON0+DLON); ax.set_ylim(LAT0-DLAT,LAT0+DLAT)
ax=axs[5]
for kk,cc,lb in [(0,'0.75','unch'),(1,'tab:blue','opened'),(2,'tab:red','advanced')]:
    m=(klass==kk)&(lon>LON0-DLON)&(lon<LON0+DLON)&(lat>LAT0-DLAT)&(lat<LAT0+DLAT); ax.scatter(lon[m],lat[m],s=9,c=cc,label=lb)
cavity_outline(ax); ax.plot(LON0,LAT0,'k*',ms=15,mec='w'); ax.legend(loc='upper left',fontsize=8); ax.set_title("node class + cavity outline")
ax.set_xlim(LON0-DLON,LON0+DLON); ax.set_ylim(LAT0-DLAT,LAT0+DLAT)
for ax in axs: ax.set_xlabel("lon"); ax.set_ylabel("lat")
fig.suptitle(f"Initial state @ level {KLEV} (~{-zmid[k]:.0f} m) — dateline tris filtered — black line = cavity (ice-front) outline — seed ({LON0},{LAT0})")
fig.tight_layout(); fig.savefig(f"{OUT}/planview_lev{KLEV}.png",dpi=115); plt.close(fig); print("planview ok",len(sel),"tris")
# ===== sections with ice-base line =====
def dosection(axis,fname):
    if axis=='lon':
        band=(np.abs(lat-LAT0)<0.10)&(lon>LON0-DLON)&(lon<LON0+DLON); coord=lon; xl=f"lon (lat={LAT0})"
    else:
        band=(np.abs(lon-LON0)<0.15)&(lat>LAT0-DLAT)&(lat<LAT0+DLAT); coord=lat; xl=f"lat (lon={LON0})"
    idx=np.where(band)[0]; idx=idx[np.argsort(coord[idx])]; xc=coord[idx]
    edges=np.concatenate([[xc[0]-(xc[1]-xc[0])/2],(xc[:-1]+xc[1:])/2,[xc[-1]+(xc[-1]-xc[-2])/2]])
    maxdep=zbar[nln[idx]-1].min(); icebase=zbar[uln[idx]-1]  # top wet level depth (=ice base for cavity)
    def secpoly(field):
        polys=[]; vals=[]
        for a,j in enumerate(idx):
            for kk in range(len(zmid)):
                if (uln[j]<=kk+1) and (nln[j]-1>=kk+1):
                    polys.append([[edges[a],zbar[kk]],[edges[a+1],zbar[kk]],[edges[a+1],zbar[kk+1]],[edges[a],zbar[kk+1]]]); vals.append(field[kk,j])
        return polys,np.array(vals)
    fig,axs=plt.subplots(2,2,figsize=(15,9)); axs=axs.ravel()
    for ax,(ttl,fld,cm,opt) in zip(axs,[("temp",T,'RdYlBu_r',None),("salt",S,'viridis',None),("speed",SPDn,'magma',None),("w",W,'RdBu_r','sym')]):
        polys,vals=secpoly(fld); pc=PolyCollection(polys,array=vals,cmap=cm,edgecolors='face',linewidths=0.1)
        if opt=='sym':
            a=np.nanmax(np.abs(vals)) or 1e-6; pc.set_clim(-a,a)
        ax.add_collection(pc); plt.colorbar(pc,ax=ax,shrink=.85)
        ax.plot(xc,icebase,'k-',lw=1.8,drawstyle='steps-mid')   # ice base / cavity outline
        ax.fill_between(xc,icebase,5,step='mid',color='none',hatch='xx',edgecolor='0.4',lw=0)  # ice region hatched
        ax.set_title(ttl+"  (black=ice base)"); ax.set_xlabel(xl); ax.set_ylabel("depth (m)")
        ax.set_xlim(edges[0],edges[-1]); ax.set_ylim(maxdep*1.05,5)
        ax.axvline(LON0 if axis=='lon' else LAT0,color='r',ls=':',lw=.8); ax.axhline(zmid[KLEV-1],color='k',ls='--',lw=.7)
    fig.suptitle(f"Initial state {axis}-depth section; bottom={-maxdep:.0f} m; black line=ice base(cavity), dashed=level {KLEV}")
    fig.tight_layout(); fig.savefig(f"{OUT}/{fname}.png",dpi=115); plt.close(fig); print(fname,"ok")
dosection('lon','section_lon'); dosection('lat','section_lat'); print("DONE")
