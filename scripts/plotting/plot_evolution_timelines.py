#!/usr/bin/env python3
"""Coupled-system evolution timelines over the movcav34 production run."""
import numpy as np, netCDF4 as nc, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/work/ab0246/a270092/runtime/awiesm3-develop-is/movcav34"
POOL = "/work/ab0246/a270092/input/fesom2/pism_cavity_ini/submesh"

# ---- mesh timeline -----------------------------------------------------------
subs = [POOL] + sorted(glob.glob(f"{E}/couple/submesh_*"))
years = [1900] + [int(s.split("submesh_")[1][:4]) for s in subs[1:]]  # submesh_Y is the mesh of year Y
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
rests = sorted(glob.glob(f"{E}/restart/pism/movcav34_pismr_restart_*.nc"))
spin = "/work/ab0246/a270193/esm-experiments/1percCO2_b_2640/restart/pism_sh/1percCO2_b_2640_pismr_restart_627500101-627501231.nc"
pyears = [0] + [10*(i+1) for i in range(len(rests))]
vol, afloat = [], []
for f in [spin] + rests:
    d = nc.Dataset(f)
    thk = np.array(d['thk'][-1] if d['thk'].ndim == 3 else d['thk'][:])
    topg = np.array(d['topg'][-1] if d['topg'].ndim == 3 else d['topg'][:])
    vol.append(thk.sum()*CELL/1e9/1e3)                       # 1000 km^3
    flo = (thk > 5) & (topg + thk*0.917 < 0)                 # floating criterion
    afloat.append(flo.sum()*CELL/1e6/1e3)                    # 1000 km^2
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
    ul = np.loadtxt(f"{sub}/cavity_nlvls.out").astype(int)
    cav = ul > 1
    d = nc.Dataset(f)
    fw = np.array(d['fw'][:]); fw[np.abs(fw) > 1e30] = np.nan
    ann = np.nanmean(fw, axis=0)                              # m/s annual mean
    if ann.shape[0] != len(ul): d.close(); continue
    na = node_areas(sub)
    m_wy = ann*86400*365                                      # m w.e./yr
    w = na[cav]
    melt_mean.append(-np.nansum(m_wy[cav]*w)/w.sum())                 # melt positive (fw is ocean-loss-positive)
    melt_tot.append(-np.nansum(ann[cav]*na[cav])*86400*365*1000/1e12)  # Gt/yr, melt positive
    myears.append(y); d.close()

# ---- figure ------------------------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(11, 10), sharex=False, facecolor="white")
c1, c2 = "#4053d3", "#b51d14"

ax = axs[0]
ax.plot(years, ntot, "-o", color=c1, label="submesh nodes")
ax.plot(years, [n*20 for n in ncav], "-^", color="#00beff", label="cavity nodes ($\\times$20)")
ax2 = ax.twinx()
ax2.bar([y-0.18 for y in years], added, width=0.36, color="#00b25d", alpha=0.8, label="nodes added")
ax2.bar([y+0.18 for y in years], removed, width=0.36, color="#b51d14", alpha=0.8, label="nodes removed")
ax.set_ylabel("submesh nodes", color=c1)
ax2.set_ylabel("nodes changed per increment")
ax.set_title("(a) The moving mesh: 12 increments, %d→%d nodes" % (ntot[0], ntot[-1]), fontsize=11)
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2, l1+l2, fontsize=8, frameon=False, loc="center left")
ax.spines[["top"]].set_visible(False)

ax = axs[1]
ax.plot(pyears, vol, "-o", color=c1, label="ice volume")
ax.set_ylabel("ice volume [10$^3$ km$^3$]", color=c1)
ax2 = ax.twinx()
ax2.plot(pyears, afloat, "-s", color=c2, label="floating area")
ax2.set_ylabel("floating area [10$^3$ km$^2$]", color=c2)
ax.set_xlabel("PISM year")
ax.set_title("(b) PISM under FESOM's own melt: volume and floating area over 120 years", fontsize=11)
ax.spines[["top"]].set_visible(False)

ax = axs[2]
ax.plot(myears, melt_mean, "-o", color=c1)
ax.set_ylabel("cavity-mean melt [m w.e. yr$^{-1}$]", color=c1)
ax.axhspan(0.3, 1.5, color="0.9", zorder=0)
ax.text(myears[0], 0.9, "obs circum-Antarctic range (Adusumilli 2020)", fontsize=7.5, color="0.4")
ax2 = ax.twinx()
ax2.plot(myears, melt_tot, "-s", color=c2)
ax2.axhline(1325, color=c2, lw=0.8, ls="--")
ax2.text(myears[-1]-3, 1370, "Rignot 2013: 1325 Gt/yr", fontsize=7.5, color=c2)
ax2.set_ylabel("total cavity melt [Gt yr$^{-1}$]", color=c2)
ax.set_xlabel("year")
ax.set_yscale("log"); ax.set_title("(c) FESOM cavity melt escalates as the shelves thin (log scale)", fontsize=11)
ax.spines[["top"]].set_visible(False)

for i,axx in enumerate(axs):
    x0 = 1902 if i != 1 else 20
    axx.axvspan(x0, axx.get_xlim()[1], color="#fdecec", zorder=0)
axs[0].text(1907, 150000, "TAINTED: cavities deleted by the\nhygiene sweep (14-111 nodes left)", color="#b51d14", fontsize=10, ha="center")
axs[1].text(70, 25550, "TAINTED", color="#b51d14", fontsize=11)
axs[2].text(1907.5, 2.5, "TAINTED", color="#b51d14", fontsize=11)
fig.suptitle("movcav34: the coupled system — real through chunk 4, cavity-less artifact after (hygiene overreach)", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.97))
for out in ("/home/a/a270092/esm_tools/movcav4_crash_plots/movcav34_evolution_timelines.png",
            "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/report/movcav34_evolution_timelines.png"):
    fig.savefig(out, dpi=140)
print("saved; nodes", ntot[0], "->", ntot[-1], "| melt_mean", [round(m,2) for m in melt_mean])
print("vol[0/last]:", round(vol[0],1), round(vol[-1],1), "| afloat[0/last]:", round(afloat[0],1), round(afloat[-1],1))
