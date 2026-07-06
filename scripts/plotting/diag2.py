import numpy as np, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
C="/work/ab0246/a270092/pism_repro/experiments/movcav8/couple"; R=f"{C}/restart_remapped_v4"
OUT="/home/a/a270092/esm_tools/movcav4_crash_plots/initstate"
LON0,LAT0=-62.96,-66.81; KLEV=11; DLON,DLAT=3.0,1.3
ao=np.loadtxt(f"{C}/latest_submesh/nod2d.out",skiprows=1); lon,lat=ao[:,1],ao[:,2]; nn=len(lon)
uln=np.loadtxt(f"{C}/latest_submesh/cavity_nlvls.out").astype(int); nln=np.loadtxt(f"{C}/latest_submesh/nlvls.out").astype(int)
elem=np.loadtxt(f"{C}/latest_submesh/elem2d.out",skiprows=1).astype(int)-1
print("elem index range:",elem.min(),elem.max(),"nnodes:",nn,"| nod2d.out first idx col?",int(ao[0,0]) if ao.shape[1]>3 else "n/a")
# check: nod2d.out — is col0 an index 1..N (so my cols 1,2 = lon,lat correct)?
print("nod2d cols sample:",ao[0],ao[1])
S=np.squeeze(np.asarray(nc.Dataset(f"{R}/salt.nc")['salt'][:])); S=S.T if S.shape[0]==nn else S
k=KLEV-1
# triangle geometry in box
ecx=lon[elem].mean(1); ecy=lat[elem].mean(1)
inbox=(ecx>LON0-DLON)&(ecx<LON0+DLON)&(ecy>LAT0-DLAT)&(ecy<LAT0+DLAT)
ev=lon[elem]; maxedge=np.maximum.reduce([np.abs(ev[:,0]-ev[:,1]),np.abs(ev[:,1]-ev[:,2]),np.abs(ev[:,0]-ev[:,2])])
print(f"in-box tris: {inbox.sum()}; max lon-span of a tri: {maxedge[inbox].max():.3f} deg; #tris spanning>1deg lon: {np.sum(inbox&(maxedge>1))}")
wet=(uln<=k+1)&(nln-1>=k+1)
def verts(sel): e=elem[sel]; return np.stack([np.column_stack([lon[e[:,i]],lat[e[:,i]]]) for i in range(3)],axis=1)
fig,axs=plt.subplots(2,2,figsize=(15,11)); axs=axs.ravel()
sel=np.where(inbox & wet[elem].all(axis=1))[0]; V=verts(sel)
# (a) SMOOTH TEST: color triangles by centroid latitude — must be a clean N-S gradient
pc=PolyCollection(V,array=ecy[sel],cmap='turbo',edgecolors='k',linewidths=0.4); axs[0].add_collection(pc)
plt.colorbar(pc,ax=axs[0]); axs[0].set_title("SMOOTH TEST: triangle centroid LATITUDE\n(if this stripes -> triangulation/order bug)")
# (b) raw scatter of salt (no triangles)
mm=wet&(lon>LON0-DLON)&(lon<LON0+DLON)&(lat>LAT0-DLAT)&(lat<LAT0+DLAT)
sc=axs[1].scatter(lon[mm],lat[mm],c=S[k,mm],s=60,cmap='viridis',edgecolors='k',linewidths=0.3)
plt.colorbar(sc,ax=axs[1]); axs[1].set_title("RAW node scatter: salt @ lev11 (no triangles)")
# (c) salt PolyCollection WITH visible edges + node index labels-ish
pc=PolyCollection(V,array=S[k][elem[sel]].mean(1),cmap='viridis',edgecolors='k',linewidths=0.5); axs[2].add_collection(pc)
plt.colorbar(pc,ax=axs[2]); axs[2].set_title("salt PolyCollection with mesh edges")
# (d) node index as color (reveals ordering structure)
sc=axs[3].scatter(lon[mm],lat[mm],c=np.where(mm)[0][np.argsort(np.argsort(np.where(mm)[0]))] if False else np.arange(nn)[mm],s=60,cmap='nipy_spectral',edgecolors='k',linewidths=0.3)
plt.colorbar(sc,ax=axs[3]); axs[3].set_title("node INDEX (ordering) — reveals if data order != geometry")
for ax in axs: ax.plot(LON0,LAT0,'r*',ms=14); ax.set_xlim(LON0-DLON,LON0+DLON); ax.set_ylim(LAT0-DLAT,LAT0+DLAT)
fig.tight_layout(); fig.savefig(f"{OUT}/diag_artifact.png",dpi=115); plt.close(fig); print("wrote diag_artifact.png")
