#!/usr/bin/env python3
"""Coupled-system evolution timelines over the movcavlinfs (linfs) run.
Same panels as report Fig. 8 (movcav34), for the linfs run, to compare the
ice-sheet change magnitude. Prints movcav34 (zstar) numbers alongside."""
import numpy as np, netCDF4 as nc, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/work/ab0246/a270092/runtime/awiesm3-develop-is/movcavlinfs"
POOL = "/work/ab0246/a270092/input/fesom2/pism_cavity_ini/submesh"

# ---- mesh timeline -----------------------------------------------------------
subs = [POOL] + sorted(glob.glob(f"{E}/couple/submesh_*"))
years = [1900] + [int(s.split("submesh_")[1][:4]) for s in subs[1:]]
ntot, ncav, mapsets = [], [], []
for s in subs:
    with open(f"{s}/nod2d.out") as f: n = int(f.readline().split()[0])
    ul = np.loadtxt(f"{s}/cavity_nlvls.out").astype(int)
    ntot.append(n); ncav.append(int((ul > 1).sum()))
    mapsets.append(set(np.loadtxt(f"{s}/map_nod.out").astype(int).tolist()))
added   = [np.nan] + [len(mapsets[i] - mapsets[i-1]) for i in range(1, len(mapsets))]
removed = [np.nan] + [len(mapsets[i-1] - mapsets[i]) for i in range(1, len(mapsets))]

# ---- PISM timeline -----------------------------------------------------------
CELL = 8000.0**2
rests = sorted(glob.glob(f"{E}/restart/pism/*_pismr_restart_*.nc"))
spin = "/work/ab0246/a270193/esm-experiments/1percCO2_b_2640/restart/pism_sh/1percCO2_b_2640_pismr_restart_627500101-627501231.nc"
pyears = [0] + [10*(i+1) for i in range(len(rests))]
vol, afloat = [], []
for f in [spin] + rests:
    d = nc.Dataset(f)
    thk = np.array(d['thk'][-1] if d['thk'].ndim == 3 else d['thk'][:])
    topg = np.array(d['topg'][-1] if d['topg'].ndim == 3 else d['topg'][:])
    vol.append(thk.sum()*CELL/1e9/1e3)
    flo = (thk > 5) & (topg + thk*0.917 < 0)
    afloat.append(flo.sum()*CELL/1e6/1e3)
    d.close()

# ---- FESOM cavity melt timeline ---------------------------------------------
def node_areas(sub):
    ao = np.loadtxt(f"{sub}/nod2d.out", skiprows=1)
    lon, lat = np.radians(ao[:,1]), np.radians(ao[:,2])
    e = np.loadtxt(f"{sub}/elem2d.out", skiprows=1).astype(int)-1
    R = 6371e3
    # spherical excess (Girard): a plane formula explodes on the ~1000 elements
    # spanning +/-180 deg, inflating the global total from 361 to 860e6 km^2
    v = np.stack([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)], axis=1)
    va, vb, vc = v[e[:,0]], v[e[:,1]], v[e[:,2]]
    _num = np.abs(np.einsum('ij,ij->i', va, np.cross(vb, vc)))
    _den = 1.0 + np.einsum('ij,ij->i',va,vb) + np.einsum('ij,ij->i',vb,vc) + np.einsum('ij,ij->i',vc,va)
    ar = 2.0*np.arctan2(_num, _den)*R*R
    na = np.zeros(len(ao))
    for k in range(3): np.add.at(na, e[:,k], ar/3.0)
    return na

melt_mean, melt_tot, myears = [], [], []
for y in range(1900, 1912):
    f = f"{E}/outdata/fesom/fw.fesom.{y}.nc"
    if not os.path.exists(f): continue
    sub = POOL if y <= 1900 else f"{E}/couple/submesh_{y}-12-31T00:00:00"
    if not os.path.exists(f"{sub}/cavity_nlvls.out"): continue
    ul = np.loadtxt(f"{sub}/cavity_nlvls.out").astype(int); cav = ul > 1
    d = nc.Dataset(f); fw = np.array(d['fw'][:]); fw[np.abs(fw) > 1e30] = np.nan
    ann = np.nanmean(fw, axis=0)
    if ann.shape[0] != len(ul): d.close(); continue
    na = node_areas(sub); m_wy = ann*86400*365; w = na[cav]
    melt_mean.append(-np.nansum(m_wy[cav]*w)/w.sum())
    melt_tot.append(-np.nansum(ann[cav]*na[cav])*86400*365*1000/1e12)
    myears.append(y); d.close()

# ---- figure ------------------------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(11, 10), facecolor="white")
c1, c2 = "#4053d3", "#b51d14"

ax = axs[0]
ax.plot(years, ntot, "-o", color=c1, label="submesh nodes")
ax.plot(years, [n*20 for n in ncav], "-^", color="#00beff", label="cavity nodes ($\\times$20)")
ax2 = ax.twinx()
ax2.bar([y-0.18 for y in years], added, width=0.36, color="#00b25d", alpha=0.8, label="nodes added")
ax2.bar([y+0.18 for y in years], removed, width=0.36, color="#b51d14", alpha=0.8, label="nodes removed")
ax.set_ylabel("submesh nodes", color=c1); ax2.set_ylabel("nodes changed per increment")
ax.set_title("(a) The moving mesh (linfs): %d→%d nodes, cavity %d→%d" % (ntot[0], ntot[-1], ncav[0], ncav[-1]), fontsize=11)
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2, l1+l2, fontsize=8, frameon=False, loc="center left"); ax.spines[["top"]].set_visible(False)

ax = axs[1]
ax.plot(pyears, vol, "-o", color=c1, label="ice volume")
ax.set_ylabel("ice volume [10$^3$ km$^3$]", color=c1)
ax2 = ax.twinx(); ax2.plot(pyears, afloat, "-s", color=c2, label="floating area")
ax2.set_ylabel("floating area [10$^3$ km$^2$]", color=c2)
ax.set_xlabel("PISM year")
ax.set_title("(b) PISM under linfs FESOM melt: volume and floating area", fontsize=11); ax.spines[["top"]].set_visible(False)

ax = axs[2]
ax.plot(myears, melt_mean, "-o", color=c1)
ax.set_ylabel("cavity-mean melt [m w.e. yr$^{-1}$]", color=c1)
ax.axhspan(0.3, 1.5, color="0.9", zorder=0)
ax.text(myears[0], 0.9, "obs circum-Antarctic range (Adusumilli 2020)", fontsize=7.5, color="0.4")
ax2 = ax.twinx(); ax2.plot(myears, melt_tot, "-s", color=c2)
ax2.set_ylabel("total cavity melt [Gt yr$^{-1}$]", color=c2)
ax.set_xlabel("year"); ax.set_title("(c) FESOM cavity melt under linfs (no escalation)", fontsize=11); ax.spines[["top"]].set_visible(False)

fig.suptitle("movcavlinfs: coupled evolution under linfs (compare Fig. 8 / movcav34 zstar)", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.97))
for out in ("/home/a/a270092/esm_tools/movcav4_crash_plots/movcavlinfs_evolution_timelines.png",
            "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/report/movcavlinfs_evolution_timelines.png"):
    fig.savefig(out, dpi=140)
print("=== movcavlinfs (linfs) ===")
print("mesh nodes:", ntot, "| cavity nodes:", ncav, "| removed/incr:", removed)
print("PISM vol[0..]:", [round(v,1) for v in vol], "| floating area[0..]:", [round(a,1) for a in afloat])
print("melt_mean (m/yr):", [round(m,2) for m in melt_mean], "| melt_tot (Gt/yr):", [round(m,1) for m in melt_tot])
print("\n=== movcav34 (zstar) for comparison, from the report ===")
print("cavity nodes collapsed 2647 -> 111 (incl hygiene amputation); melt_mean escalated 0.9 -> 5.5 m/yr")
