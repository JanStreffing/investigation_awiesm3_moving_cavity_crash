#!/usr/bin/env python3
"""fgdbg02 equilibration: 24 coupled years / 240 PISM years after the Bug VIII fix.

Four single-axis panels (no twin axes -- two y-scales on one frame is the most
misread chart form there is; each measure gets its own frame instead):

  (a) PISM ice volume            -- is the ice sheet approaching a steady state?
  (b) PISM floating area         -- are the shelves settling?
  (c) FESOM cavity-mean melt     -- does the melt stabilise inside the observed band?
  (d) mesh nodes changed / cycle -- does the geometry increment shrink?

Panel (d) is the equilibration diagnostic that matters for the coupling itself:
if added/removed per increment decays, successive submeshes are converging.
"""
import glob, os
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02"
POOL = "/work/ab0246/a270092/input/fesom2/pism_cavity_ini/submesh"
OUT = "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/report"

BLUE, ORANGE = "#2a78d6", "#eb6834"          # validated pair (dataviz slots 1,2)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d8d8d4"

# ---- mesh timeline -----------------------------------------------------------
subs = [POOL] + sorted(glob.glob(f"{E}/couple/submesh_*"))
myears = [1900] + [int(s.split("submesh_")[1][:4]) for s in subs[1:]]
ntot, ncav, mapsets = [], [], []
for s in subs:
    with open(f"{s}/nod2d.out") as f:
        ntot.append(int(f.readline().split()[0]))
    ul = np.loadtxt(f"{s}/cavity_nlvls.out").astype(int)
    ncav.append(int((ul > 1).sum()))
    mapsets.append(set(np.loadtxt(f"{s}/map_nod.out").astype(int).tolist()))
added = [len(mapsets[i] - mapsets[i - 1]) for i in range(1, len(mapsets))]
removed = [len(mapsets[i - 1] - mapsets[i]) for i in range(1, len(mapsets))]
inc_years = myears[1:]

# ---- PISM timeline -----------------------------------------------------------
CELL = 8000.0 ** 2
rests = sorted(glob.glob(f"{E}/restart/pism/*_pismr_restart_*.nc"))
spin = ("/work/ab0246/a270193/esm-experiments/1percCO2_b_2640/restart/pism_sh/"
        "1percCO2_b_2640_pismr_restart_627500101-627501231.nc")
pyears, vol, afloat = [], [], []
for i, f in enumerate([spin] + rests):
    if not os.path.exists(f):
        continue
    d = nc.Dataset(f)
    thk = np.array(d["thk"][-1] if d["thk"].ndim == 3 else d["thk"][:])
    topg = np.array(d["topg"][-1] if d["topg"].ndim == 3 else d["topg"][:])
    vol.append(thk.sum() * CELL / 1e9 / 1e3)
    flo = (thk > 5) & (topg + thk * 0.917 < 0)
    afloat.append(flo.sum() * CELL / 1e6 / 1e3)
    pyears.append(10 * i)
    d.close()

# ---- FESOM cavity melt -------------------------------------------------------
def node_areas(sub):
    ao = np.loadtxt(f"{sub}/nod2d.out", skiprows=1)
    lon, lat = np.radians(ao[:, 1]), np.radians(ao[:, 2])
    e = np.loadtxt(f"{sub}/elem2d.out", skiprows=1).astype(int) - 1
    R = 6371e3
    # Spherical excess (Girard). A plane-projection formula explodes on the ~1000
    # elements spanning +/-180 deg (lon span ~2*pi instead of ~0), which inflates
    # the global total from 361 to 860e6 km^2.
    v = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon),
                  np.sin(lat)], axis=1)
    va, vb, vc = v[e[:, 0]], v[e[:, 1]], v[e[:, 2]]
    _num = np.abs(np.einsum('ij,ij->i', va, np.cross(vb, vc)))
    _den = (1.0 + np.einsum('ij,ij->i', va, vb)
                + np.einsum('ij,ij->i', vb, vc)
                + np.einsum('ij,ij->i', vc, va))
    ar = 2.0 * np.arctan2(_num, _den) * R * R
    na = np.zeros(len(ao))
    for k in range(3):
        np.add.at(na, e[:, k], ar / 3.0)
    return na

melt_mean, melt_years = [], []
for y in range(1900, 1930):
    f = f"{E}/outdata/fesom/fw.fesom.{y}.nc"
    if not os.path.exists(f):
        continue
    sub = POOL if y <= 1900 else f"{E}/couple/submesh_{y}-12-31T00:00:00"
    if not os.path.exists(f"{sub}/cavity_nlvls.out"):
        continue
    ul = np.loadtxt(f"{sub}/cavity_nlvls.out").astype(int)
    cav = ul > 1
    d = nc.Dataset(f)
    fw = np.array(d["fw"][:])
    d.close()
    fw[np.abs(fw) > 1e30] = np.nan
    ann = np.nanmean(fw, axis=0)
    if ann.shape[0] != len(ul):
        continue
    na = node_areas(sub)
    w = na[cav]
    melt_mean.append(-np.nansum(ann[cav] * 86400 * 365 * w) / w.sum())
    melt_years.append(y)

# ---- figure ------------------------------------------------------------------
def tidy(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=9)

fig, axs = plt.subplots(2, 2, figsize=(12.5, 7.6), facecolor="white")
fig.suptitle("fgdbg02 after the Bug VIII fix: 24 coupled years / 240 PISM years, no blow-ups",
             fontsize=13, color=INK, y=0.98)

ax = axs[0, 0]
ax.plot(pyears, vol, "-o", color=BLUE, lw=2, ms=5, zorder=3)
ax.set_ylabel("ice volume [10$^3$ km$^3$]", color=MUTED, fontsize=9.5)
ax.set_xlabel("PISM year", color=MUTED, fontsize=9.5)
ax.set_title("(a) PISM ice volume", fontsize=10.5, color=INK, loc="left")
if len(vol) > 1:
    ax.text(0.98, 0.93, f"{vol[0]:.0f} → {vol[-1]:.0f}  ({vol[-1]-vol[0]:+.0f})",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color=MUTED)
tidy(ax)

ax = axs[0, 1]
ax.plot(pyears, afloat, "-o", color=BLUE, lw=2, ms=5, zorder=3)
ax.set_ylabel("floating area [10$^3$ km$^2$]", color=MUTED, fontsize=9.5)
ax.set_xlabel("PISM year", color=MUTED, fontsize=9.5)
ax.set_title("(b) PISM floating (ice-shelf) area", fontsize=10.5, color=INK, loc="left")
if len(afloat) > 1:
    ax.text(0.98, 0.93, f"{afloat[0]:.0f} → {afloat[-1]:.0f}  ({afloat[-1]-afloat[0]:+.0f})",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color=MUTED)
tidy(ax)

ax = axs[1, 0]
ax.axhspan(0.3, 1.5, color="#eceae4", zorder=0)
ax.plot(melt_years, melt_mean, "-o", color=BLUE, lw=2, ms=5, zorder=3)
ax.set_ylabel("cavity-mean melt [m w.e. yr$^{-1}$]", color=MUTED, fontsize=9.5)
ax.set_xlabel("year", color=MUTED, fontsize=9.5)
ax.set_title("(c) FESOM cavity melt", fontsize=10.5, color=INK, loc="left")
ax.set_ylim(0.0, 1.62)
if melt_years:
    ax.text(0.03, 0.94, "shaded: observed circum-Antarctic range (Adusumilli 2020)",
            transform=ax.transAxes, fontsize=8, color=MUTED, va="top")
tidy(ax)

ax = axs[1, 1]
w = 0.4
ax.bar([y - w / 2 for y in inc_years], added, width=w - 0.04, color=BLUE,
       label="nodes added", zorder=3)
ax.bar([y + w / 2 for y in inc_years], removed, width=w - 0.04, color=ORANGE,
       label="nodes removed", zorder=3)
ax.set_ylabel("mesh nodes changed", color=MUTED, fontsize=9.5)
ax.set_xlabel("year of increment", color=MUTED, fontsize=9.5)
ax.set_title("(d) Mesh increment per coupling cycle", fontsize=10.5, color=INK, loc="left")
ax.legend(fontsize=8.5, frameon=False, labelcolor=MUTED, loc="upper right")
ax.text(0.03, 0.94, f"removed/cycle {removed[0]} → {int(np.mean(removed[-5:]))} (last-5 mean)",
        transform=ax.transAxes, fontsize=8, color=MUTED, va="top")
tidy(ax)

fig.tight_layout(rect=(0, 0, 1, 0.955))
fig.savefig(f"{OUT}/fgdbg02_equilibration.png", dpi=150, facecolor="white")

print("cycles (submeshes):", len(subs), "| mesh nodes:", ntot[0], "->", ntot[-1],
      "| cavity nodes:", ncav[0], "->", ncav[-1])
print("PISM years:", pyears[0], "->", pyears[-1],
      "| volume:", round(vol[0], 1), "->", round(vol[-1], 1),
      "| floating area:", round(afloat[0], 1), "->", round(afloat[-1], 1))
print("melt years:", melt_years[0] if melt_years else None, "->",
      melt_years[-1] if melt_years else None,
      "| melt:", [round(m, 2) for m in melt_mean])
print("nodes added/incr :", added)
print("nodes removed/incr:", removed)
