import numpy as np, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
C="/work/ab0246/a270092/pism_repro/experiments/movcav8/couple"; R=f"{C}/restart_remapped_v4"
OUT="/home/a/a270092/esm_tools/movcav4_crash_plots/initstate"
LON0,LAT0=-62.96,-66.81; DLON,DLAT=3.0,1.3; KLEV=11
ao=np.loadtxt(f"{C}/latest_submesh/nod2d.out",skiprows=1); lon,lat=ao[:,1],ao[:,2]
uln=np.loadtxt(f"{C}/latest_submesh/cavity_nlvls.out").astype(int)
nln=np.loadtxt(f"{C}/latest_submesh/nlvls.out").astype(int)
cav=np.loadtxt(f"{C}/latest_submesh/cavity_depth@node.out")  # ice draft
box=(lon>LON0-DLON)&(lon<LON0+DLON)&(lat>LAT0-DLAT)&(lat<LAT0+DLAT)
S=np.squeeze(np.asarray(nc.Dataset(f"{R}/salt.nc")['salt'][:]))
if S.shape[0]==len(lon): S=S.T
k=KLEV-1; wet=(uln<=k+1)&(nln-1>=k+1)
fig,axs=plt.subplots(1,3,figsize=(19,6))
m=box
sc=axs[0].scatter(lon[m],lat[m],c=uln[m],s=14,cmap='turbo',vmin=1,vmax=16)
plt.colorbar(sc,ax=axs[0]); axs[0].set_title("ulevels (cavity_nlvls): 1=open surface, >1=under ice")
sc=axs[1].scatter(lon[m],lat[m],c=cav[m],s=14,cmap='viridis')
plt.colorbar(sc,ax=axs[1]); axs[1].set_title("cavity_depth@node (ice draft, m)")
mm=box&wet
sc=axs[2].scatter(lon[mm],lat[mm],c=S[k,mm],s=14,cmap='viridis')
plt.colorbar(sc,ax=axs[2]); axs[2].set_title(f"salt @ level {KLEV} (only nodes wet at that level)")
for ax in axs: ax.plot(LON0,LAT0,'k*',ms=13,mec='w'); ax.set_xlim(LON0-DLON,LON0+DLON); ax.set_ylim(LAT0-DLAT,LAT0+DLAT); ax.set_xlabel('lon'); ax.set_ylabel('lat')
fig.tight_layout(); fig.savefig(f"{OUT}/diag_mask_stripes.png",dpi=115); plt.close(fig)
# numeric: are there horizontal rows of ulev==1 among ulev>1? check unique lats with many nodes
print("nodes in box:",box.sum(),"| ulev==1:",np.sum(box&(uln==1)),"| ulev>1:",np.sum(box&(uln>1)))
# for the warm/salty stripe: look at nodes that are wet-at-11 & salty (>33.5) in the western half where most is cavity
west=box&(lon<LON0)&wet
salty=west&(S[k]>33.5)
print(f"western wet-at-lev11 nodes: {west.sum()}, of which salty(>33.5): {salty.sum()}")
if salty.sum():
    la=lat[salty]; print("  their latitudes (rounded):",np.round(np.sort(np.unique(la)),3)[:15])
    print("  their ulev:",np.unique(uln[salty]),"  their class vs old? need old")
print("wrote diag_mask_stripes.png")
