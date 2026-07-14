#!/usr/bin/env python
"""PISM geometry change over one 10-yr coupling leg, on the FESOM mesh.

Two figures, same idioms as plot_layers.py (PolyCollection on the real mesh
triangles, NOWRAP dateline filter, cavity outline):

  pism_change_planview.png  -- top view: d(ice draft) and d(ulevel), old/new cavity fronts
  pism_change_section.png   -- cutaway from the side through the Ross cluster that rang

OLD mesh = the submesh chunk 1 ran on (pool).  NEW mesh = the submesh built from
PISM's evolved geometry for chunk 3.
"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation

NEW = "/work/ab0246/a270092/pism_repro/experiments/movcav12/couple/latest_submesh"
OLD = "/work/ab0246/a270092/input/fesom2/pism_cavity_ini/submesh"
OUT = "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/figures/pism_change"
import os; os.makedirs(OUT, exist_ok=True)

# ---- meshes (same loading pattern as plot_layers.py) -------------------------
ao   = np.loadtxt(f"{NEW}/nod2d.out", skiprows=1); lon, lat = ao[:,1], ao[:,2]; nn = len(lon)
uln  = np.loadtxt(f"{NEW}/cavity_nlvls.out").astype(int)
nln  = np.loadtxt(f"{NEW}/nlvls.out").astype(int)
drN  = np.loadtxt(f"{NEW}/cavity_depth@node.out")          # ice draft, <=0
mapn = np.loadtxt(f"{NEW}/map_nod.out").astype(int)
elem = np.loadtxt(f"{NEW}/elem2d.out", skiprows=1).astype(int) - 1

uO   = np.loadtxt(f"{OLD}/cavity_nlvls.out").astype(int)
drO_ = np.loadtxt(f"{OLD}/cavity_depth@node.out")
mapo = np.loadtxt(f"{OLD}/map_nod.out").astype(int)
b2o  = {b: i for i, b in enumerate(mapo)}

with open(f"{NEW}/aux3d.out") as f:
    nlv = int(f.readline()); zbar = np.array([float(f.readline()) for _ in range(nlv)])

# map OLD -> NEW node order (via the shared mother-mesh base index)
ulev_old = np.full(nn, -1, int); draft_old = np.full(nn, np.nan)
for j in range(nn):
    i = b2o.get(mapn[j])
    if i is not None:
        ulev_old[j] = uO[i]; draft_old[j] = drO_[i]
isnew  = ulev_old < 0
dulev  = np.where(isnew, 0, uln - ulev_old)                 # +ve: ice grew down; -ve: ice removed
ddraft = np.where(isnew, 0.0, drN - draft_old)              # +ve: draft shallower (ice lost)

# ---- mesh geometry helpers (lifted from plot_layers.py) ----------------------
ev = lon[elem]; NOWRAP = (ev.max(1) - ev.min(1)) < 180
ecx, ecy = lon[elem].mean(1), lat[elem].mean(1)
trg = Triangulation(lon, lat, elem); trg.set_mask(~NOWRAP)
def verts(sel):
    e = elem[sel]
    return np.stack([np.column_stack([lon[e[:,i]], lat[e[:,i]]]) for i in range(3)], axis=1)
def cavity_front(ax, ulev, color, lw, label):
    """the 1.5 contour of ulevels = the cavity front (ulev>1 <=> under ice)"""
    t = Triangulation(lon, lat, elem); t.set_mask(~NOWRAP)
    try:
        cs = ax.tricontour(t, ulev.astype(float), levels=[1.5], colors=color, linewidths=lw)
        cs.collections[0].set_label(label)
    except Exception:
        pass

# ============================================================ FIG 1: top view
SEL = NOWRAP & (ecy < -60)
fig, axs = plt.subplots(1, 2, figsize=(17, 7.4))

# (a) draft change
ax = axs[0]
arr = ddraft[elem[SEL]].mean(1)
pc = PolyCollection(verts(np.where(SEL)[0]), array=arr, cmap="RdBu_r",
                    edgecolors="face", linewidths=0.05)
pc.set_clim(-300, 300); ax.add_collection(pc)
plt.colorbar(pc, ax=ax, shrink=0.8, label="$\\Delta$ ice draft  (m)   [+ = ice lost]")
cavity_front(ax, ulev_old.astype(float), "k",  1.4, "cavity front, before")
cavity_front(ax, uln.astype(float),      "lime", 1.0, "cavity front, after")
ax.set_title("PISM ice-draft change over one 10-yr leg\n(chunk-1 mesh $\\to$ chunk-3 mesh)")

# (b) ulevel change -- what the OCEAN actually sees
ax = axs[1]
arr = dulev[elem[SEL]].mean(1)
pc = PolyCollection(verts(np.where(SEL)[0]), array=arr, cmap="RdBu_r",
                    edgecolors="face", linewidths=0.05)
pc.set_clim(-20, 20); ax.add_collection(pc)
plt.colorbar(pc, ax=ax, shrink=0.8, label="$\\Delta$ ulevel  (levels)   [- = levels opened]")
cavity_front(ax, ulev_old.astype(float), "k",  1.4, "before")
cavity_front(ax, uln.astype(float),      "lime", 1.0, "after")
# mark the nodes that rang (all had ulev 17 -> 1)
rang = (~isnew) & (ulev_old >= 10) & (uln == 1)
ax.plot(lon[rang], lat[rang], "x", color="magenta", ms=5, mew=1.4,
        label=f"fully opened ({rang.sum()} nodes)")
ax.legend(loc="lower left", fontsize=8)
ax.set_title("What the ocean sees: change in top wet level\n"
             "magenta = ice removed ENTIRELY (these are the nodes that rang)")

for ax in axs:
    ax.set_xlim(-180, 180); ax.set_ylim(-90, -60)
    ax.set_xlabel("lon"); ax.set_ylabel("lat"); ax.set_facecolor("0.93")
fig.tight_layout(); fig.savefig(f"{OUT}/pism_change_planview.png", dpi=115); plt.close(fig)
print(f"  wrote {OUT}/pism_change_planview.png")

# ======================================================= FIG 2: side cutaway
# transect through the Ross cluster that rang (~ -166 E, -82 S), along latitude
LON0, BAND = -166.0, 1.2
band = np.abs(((lon - LON0 + 180) % 360) - 180) < BAND
band &= (lat < -76) & (lat > -85)
o = np.argsort(lat[band]); idx = np.where(band)[0][o]

zs = lat[idx]
ice_old = np.where(ulev_old[idx] > 1, draft_old[idx], 0.0)   # ice base before (<=0)
ice_new = np.where(uln[idx]      > 1, drN[idx],       0.0)   # ice base after
bed     = zbar[nln[idx] - 1]                                  # seabed

fig, ax = plt.subplots(figsize=(13.5, 6.2))
ax.fill_between(zs, bed, -6000, color="0.55", zorder=1)                     # bedrock
ax.fill_between(zs, ice_old, bed, color="#cfe6f7", zorder=2, label="ocean, before")
newly = (ulev_old[idx] > 1) & (uln[idx] == 1)
ax.fill_between(zs, 0, ice_old, where=(ice_old < 0), color="#e8eef2",
                zorder=3, hatch="///", edgecolor="0.6", lw=0.0,
                label="ice shelf, before")
ax.plot(zs, ice_old, "k-",  lw=2.0, zorder=6, label="ice base, before")
ax.plot(zs, ice_new, "-",   lw=2.0, color="lime", zorder=6, label="ice base, after")
ax.plot(zs, bed,     "-",   lw=1.4, color="0.25", zorder=6, label="seabed")
# highlight the water column that opened up in ONE leg
ax.fill_between(zs, 0, ice_old, where=newly, color="crimson", alpha=0.35, zorder=5,
                label="water column opened in ONE 10-yr leg")
ax.set_xlim(zs.min(), zs.max()); ax.set_ylim(-1400, 60)
ax.set_xlabel("latitude  (transect at lon $\\approx$ %.0f$^\\circ$E)" % LON0)
ax.set_ylabel("depth (m)")
ax.set_title("Cutaway through the Ross cluster that rang\n"
             "the ice shelf is removed outright: ~16 ocean levels open in a single leg")
ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{OUT}/pism_change_section.png", dpi=115); plt.close(fig)
print(f"  wrote {OUT}/pism_change_section.png")

# ---- numbers for the caption ------------------------------------------------
chg = dulev[(~isnew) & (dulev != 0)]
print(f"\n  changed nodes        : {len(chg)}  ({100*len(chg)/nn:.2f}% of mesh)")
print(f"  |dulev| >= 10        : {np.sum(np.abs(chg) >= 10)}")
print(f"  fully opened (ulev->1 from >=10): {rang.sum()}")
print(f"  draft change: mean|d|={np.nanmean(np.abs(ddraft[~isnew])):.1f} m   max|d|={np.nanmax(np.abs(ddraft[~isnew])):.0f} m")
