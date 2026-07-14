#!/usr/bin/env python
"""West Weddell (Filchner-Ronne) transect: does spp keep the cavity cold?

Same transect Xiaojie used, but drawn cavity-first: the ice shelf is hatched, the
bedrock grey, and the ice base / seabed drawn explicitly -- the same rendering as the
PISM-change cutaway. His plot put a 6000 m ocean on a log axis and squashed the cavity
(the whole point) into the left edge.

Panels: CORE3 spp OFF (our configuration) vs spp ON, plus our own coupled run.
"""
import numpy as np, netCDF4 as nc, matplotlib, os
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

OUT = "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/figures/spp"
os.makedirs(OUT, exist_ok=True)
CORE3 = "/work/ab0246/a270092/input/fesom2/core3"

# ---- Xiaojie's transect: starts UNDER the Filchner-Ronne, runs out to the deep Weddell
TL = [-80,-70,-60,-50,-40,-35]
TA = [-80,-79,-77,-75,-70,-65]

def load_mesh(mp):
    ao  = np.loadtxt(f"{mp}/nod2d.out", skiprows=1)
    uln = np.loadtxt(f"{mp}/cavity_nlvls.out").astype(int)
    nln = np.loadtxt(f"{mp}/nlvls.out").astype(int)
    dr  = np.abs(np.loadtxt(f"{mp}/cavity_depth@node.out"))
    with open(f"{mp}/aux3d.out") as f:
        nlv = int(f.readline()); zbar = np.array([float(f.readline()) for _ in range(nlv)])
    return ao[:,1], ao[:,2], uln, nln, dr, zbar

def transect_idx(lon, lat, npts=340):
    pts=[]
    for i in range(len(TL)-1):
        for f in np.linspace(0,1,npts//(len(TL)-1),endpoint=False):
            pts.append((TL[i]+f*(TL[i+1]-TL[i]), TA[i]+f*(TA[i+1]-TA[i])))
    pts.append((TL[-1],TA[-1])); pts=np.array(pts)
    idx = cKDTree(np.c_[lon,lat]).query(pts)[1]
    R=6371.; d=np.zeros(len(pts))
    for i in range(1,len(pts)):
        dla=np.radians(pts[i,1]-pts[i-1,1]); dlo=np.radians(pts[i,0]-pts[i-1,0])
        a=np.sin(dla/2)**2+np.cos(np.radians(pts[i,1]))*np.cos(np.radians(pts[i-1,1]))*np.sin(dlo/2)**2
        d[i]=d[i-1]+2*R*np.arcsin(np.sqrt(a))
    return idx, d

def get_T(path, var="temp", nnode=None):
    a=np.squeeze(np.asarray(nc.Dataset(path)[var][:]))
    if a.ndim==3: a=a.mean(0)                 # time-mean
    if a.shape[0]==nnode: a=a.T               # -> (nz, nnode)
    return a

def panel(ax, T, idx, d, uln, nln, dr, zbar, title, vmin=-2.2, vmax=1.0):
    nz   = T.shape[0]
    zmid = 0.5*(zbar[:-1]+zbar[1:])
    # build the (depth, distance) temperature section, masked outside the water column
    Z = np.full((nz, len(idx)), np.nan)
    for c,j in enumerate(idx):
        lo, hi = uln[j]-1, min(nln[j]-1, nz)
        if hi>lo: Z[lo:hi, c] = T[lo:hi, j]
    ice = np.where(uln[idx]>1, dr[idx], 0.0)          # ice base (positive depth)
    bed = np.abs(zbar[nln[idx]-1])                     # seabed
    pc = ax.pcolormesh(d, np.abs(zmid[:nz]), Z, cmap="RdYlBu_r", vmin=vmin, vmax=vmax,
                       shading="nearest", zorder=2)
    # ---- draw the cavity exactly like the PISM-change cutaway ----
    ax.fill_between(d, bed, 2000, color="0.55", zorder=3)                      # bedrock
    ax.fill_between(d, 0, ice, where=(ice>0), color="#e8eef2", zorder=4,
                    hatch="///", edgecolor="0.55", lw=0.0)                     # ICE SHELF
    ax.plot(d, ice, "k-", lw=2.0, zorder=6)                                    # ice base
    ax.plot(d, bed, "-", lw=1.4, color="0.2", zorder=6)                        # seabed
    cavend = d[uln[idx]>1].max()
    ax.axvline(cavend, color="k", ls=":", lw=1.2, zorder=7)
    ax.text(cavend*0.5, 1450, "CAVITY (Filchner-Ronne)", ha="center", fontsize=9, weight="bold", zorder=8)
    ax.text(cavend+60, 1450, "open shelf $\\rightarrow$ deep Weddell", fontsize=9, zorder=8)
    ax.set_ylim(1600, -20); ax.set_xlim(0, d.max())
    ax.set_title(title, fontsize=11)
    ax.set_ylabel("depth (m)")
    return pc

lon, lat, uln, nln, dr, zbar = load_mesh(CORE3)
idx, d = transect_idx(lon, lat)
nn = len(lon)

RUNS = [
  ("/work/ab0246/a270234/runtime/fesom-2.7/test001_core3_momixon_sppoff/outdata/fesom/temp.fesom.1981.nc",
   "CORE3, spp OFF  —  OUR CONFIGURATION", CORE3),
  ("/work/ab0246/a270234/runtime/fesom-2.7/test001_core3_momixon_sppon/outdata/fesom/temp.fesom.1981.nc",
   "CORE3, spp ON  —  the fix", CORE3),
]
fig, axs = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
for ax,(p,t,mp) in zip(axs, RUNS):
    T = get_T(p, nnode=nn)
    pc = panel(ax, T, idx, d, uln, nln, dr, zbar, t)
axs[-1].set_xlabel("distance along transect (km)")
cb = fig.colorbar(pc, ax=axs, shrink=0.85, pad=0.02)
cb.set_label("temperature (°C)")
fig.suptitle("Does the salt plume parameterisation keep the cavity cold?\n"
             "West Weddell transect (Xiaojie's), drawn cavity-first", fontsize=13, y=0.97)
fig.savefig(f"{OUT}/spp_transect_cavity.png", dpi=115, bbox_inches="tight")
print(f"  wrote {OUT}/spp_transect_cavity.png")

# ---- numbers: cavity water in each ----
print("\n  CAVITY water along the transect (ice-base level):")
cavj = idx[uln[idx]>1]
a,b,c = -0.0575, 0.0901, -7.61e-4
for p,t,mp in RUNS:
    T = get_T(p, nnode=nn)
    S = get_T(p.replace("temp","salt"), var="salt", nnode=nn)
    Tb=np.array([T[uln[j]-1,j] for j in cavj]); Sb=np.array([S[uln[j]-1,j] for j in cavj])
    drv = Tb-(b+a*Sb+c*dr[cavj])
    print(f"    {t:42s} T={Tb.mean():+6.2f} degC   thermal driving={drv.mean():+5.2f} K")
