import numpy as np, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation
from scipy.spatial import cKDTree
C="/work/ab0246/a270092/pism_repro/experiments/movcav8/couple"; R=f"{C}/restart_remapped_v4"
OUT="/home/a/a270092/esm_tools/movcav4_crash_plots/initstate"
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
ev=lon[elem]; NOWRAP=(ev.max(1)-ev.min(1))<180
ecx=lon[elem].mean(1); ecy=lat[elem].mean(1)
trg=Triangulation(lon,lat,elem); trg.set_mask(~NOWRAP)
def wet(k): return (uln<=k+1)&(nln-1>=k+1)
def verts(sel): e=elem[sel]; return np.stack([np.column_stack([lon[e[:,i]],lat[e[:,i]]]) for i in range(3)],axis=1)
def outline(ax):
    try: ax.tricontour(trg,uln.astype(float),levels=[1.5],colors='k',linewidths=1.6)
    except Exception: pass
def make(lon0,lat0,klev,tag,cfl):
    DLON=3.0/np.cos(np.radians(lat0))*np.cos(np.radians(66)); DLAT=1.3   # keep ~same km
    DLON=min(DLON,6)
    k=klev-1
    sel=np.where(NOWRAP&(ecx>lon0-DLON)&(ecx<lon0+DLON)&(ecy>lat0-DLAT)&(ecy<lat0+DLAT)&wet(k)[elem].all(axis=1))[0]
    if len(sel)<8: print(f"{tag}: only {len(sel)} tris wet at lev{klev} — sparse"); 
    Vv=verts(sel); tn=lambda f:f[k][elem[sel]].mean(1)
    fig,axs=plt.subplots(2,3,figsize=(17,9)); axs=axs.ravel()
    for ax,(ttl,arr,cm,opt) in zip(axs,[("temp (C)",tn(T),'RdYlBu_r',None),("salt (psu)",tn(S),'viridis',None),
            ("speed (m/s)",SPDe[k][sel],'magma',None),("w (m/s)",tn(W),'RdBu_r','sym'),("hnode (m)",tn(H),'cividis',None)]):
        if len(sel):
            pc=PolyCollection(Vv,array=arr,cmap=cm,edgecolors='face',linewidths=0.1)
            if opt=='sym': a=np.nanmax(np.abs(arr)) or 1e-6; pc.set_clim(-a,a)
            ax.add_collection(pc); plt.colorbar(pc,ax=ax,shrink=.85)
        outline(ax); ax.set_title(ttl); ax.plot(lon0,lat0,'k*',ms=15,mec='w')
        ax.set_xlim(lon0-DLON,lon0+DLON); ax.set_ylim(lat0-DLAT,lat0+DLAT)
    ax=axs[5]
    for kk,cc,lb in [(0,'0.75','unch'),(1,'tab:blue','opened'),(2,'tab:red','advanced')]:
        m=(klass==kk)&(lon>lon0-DLON)&(lon<lon0+DLON)&(lat>lat0-DLAT)&(lat<lat0+DLAT); ax.scatter(lon[m],lat[m],s=9,c=cc,label=lb)
    outline(ax); ax.plot(lon0,lat0,'k*',ms=15,mec='w'); ax.legend(loc='upper left',fontsize=8); ax.set_title("node class + cavity outline")
    ax.set_xlim(lon0-DLON,lon0+DLON); ax.set_ylim(lat0-DLAT,lat0+DLAT)
    for ax in axs: ax.set_xlabel("lon"); ax.set_ylabel("lat")
    fig.suptitle(f"{tag}: high-CFLz {cfl} @ ({lon0},{lat0}) level {klev} (~{-zmid[k]:.0f} m) — black=cavity outline")
    fig.tight_layout(); fig.savefig(f"{OUT}/planview_{tag}.png",dpi=110); plt.close(fig)
    # sections
    for axis,fn in [('lon',f'section_lon_{tag}'),('lat',f'section_lat_{tag}')]:
        if axis=='lon': band=(np.abs(lat-lat0)<0.10)&(lon>lon0-DLON)&(lon<lon0+DLON); coord=lon; xl=f"lon (lat={lat0})"
        else: band=(np.abs(lon-lon0)<0.18)&(lat>lat0-DLAT)&(lat<lat0+DLAT); coord=lat; xl=f"lat (lon={lon0})"
        idx=np.where(band)[0]
        if len(idx)<3: continue
        idx=idx[np.argsort(coord[idx])]; xc=coord[idx]
        edges=np.concatenate([[xc[0]-(xc[1]-xc[0])/2],(xc[:-1]+xc[1:])/2,[xc[-1]+(xc[-1]-xc[-2])/2]])
        maxdep=zbar[nln[idx]-1].min(); icebase=zbar[uln[idx]-1]
        def secpoly(field):
            P=[];val=[]
            for a,j in enumerate(idx):
                for kk in range(len(zmid)):
                    if (uln[j]<=kk+1)and(nln[j]-1>=kk+1): P.append([[edges[a],zbar[kk]],[edges[a+1],zbar[kk]],[edges[a+1],zbar[kk+1]],[edges[a],zbar[kk+1]]]); val.append(field[kk,j])
            return P,np.array(val)
        fig,axs=plt.subplots(2,2,figsize=(15,9)); axs=axs.ravel()
        for ax,(ttl,fld,cm,opt) in zip(axs,[("temp",T,'RdYlBu_r',None),("salt",S,'viridis',None),("speed",SPDn,'magma',None),("w",W,'RdBu_r','sym')]):
            P,val=secpoly(fld); 
            if len(P):
                pc=PolyCollection(P,array=val,cmap=cm,edgecolors='face',linewidths=0.1)
                if opt=='sym': a=np.nanmax(np.abs(val)) or 1e-6; pc.set_clim(-a,a)
                ax.add_collection(pc); plt.colorbar(pc,ax=ax,shrink=.85)
            ax.plot(xc,icebase,'k-',lw=1.8,drawstyle='steps-mid')
            ax.set_title(ttl+" (black=ice base)"); ax.set_xlabel(xl); ax.set_ylabel("depth (m)")
            ax.set_xlim(edges[0],edges[-1]); ax.set_ylim(maxdep*1.05,5)
            ax.axhline(zmid[k],color='k',ls='--',lw=.7); ax.axvline(lon0 if axis=='lon' else lat0,color='r',ls=':',lw=.8)
        fig.suptitle(f"{tag}: {axis}-depth section @ ({lon0},{lat0}); bottom={-maxdep:.0f}m; dashed=level {klev}")
        fig.tight_layout(); fig.savefig(f"{OUT}/{fn}.png",dpi=110); plt.close(fig)
    print(f"{tag} done ({len(sel)} tris)")
for tag,lo,la,kl,cf in [("A_90E_L38",90.53,-63.45,38,"8.45"),("B_109E_L11",109.25,-66.07,11,"7.76"),("C_163E_ross_L11",163.44,-69.81,11,"4.21")]:
    make(lo,la,kl,tag,cf)
print("ALL DONE")
