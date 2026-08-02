#!/usr/bin/env python3
"""Sub-ice-shelf cavity barotropic current: zstar (movcav34) vs linfs (movcav37).
Same coupled moving-cavity setup, differing only in which_ale. Reuses the
per-year submesh + node_areas + cavity-mask machinery from
plot_evolution_timelines.py. Cavity-mean barotropic speed = area-weighted mean
over cavity nodes of |depth-mean(u,v)|, per year, on that year's moving submesh.
"""
import numpy as np, netCDF4 as nc, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/work/ab0246/a270092/runtime/awiesm3-develop-is"
POOL = "/work/ab0246/a270092/input/fesom2/pism_cavity_ini/submesh"
FILL = 1e30
RUNS = {
    "movcav34": dict(label="zstar (movcav34)",  color="#b51d14"),
    "movcav37": dict(label="linfs (movcav37)",  color="#4053d3"),
}
YEARS = list(range(1900, 1912))

def node_areas(sub):
    ao = np.loadtxt(f"{sub}/nod2d.out", skiprows=1)
    lon, lat = np.radians(ao[:, 1]), np.radians(ao[:, 2])
    e = np.loadtxt(f"{sub}/elem2d.out", skiprows=1).astype(int) - 1
    R = 6371e3
    # spherical excess (Girard): plane formula explodes on dateline-crossing elems
    v = np.stack([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)], axis=1)
    va, vb, vc = v[e[:, 0]], v[e[:, 1]], v[e[:, 2]]
    _num = np.abs(np.einsum('ij,ij->i', va, np.cross(vb, vc)))
    _den = 1.0 + np.einsum('ij,ij->i',va,vb) + np.einsum('ij,ij->i',vb,vc) + np.einsum('ij,ij->i',vc,va)
    ar = 2.0*np.arctan2(_num, _den)*R*R
    na = np.zeros(len(ao))
    for k in range(3):
        np.add.at(na, e[:, k], ar/3.0)
    return na

def submesh(E, y):
    return POOL if y <= 1900 else f"{E}/couple/submesh_{y}-12-31T00:00:00"

def cavity_mean_speed(E, y):
    sub = submesh(E, y)
    if not os.path.exists(f"{sub}/cavity_nlvls.out"):
        return np.nan
    ul = np.loadtxt(f"{sub}/cavity_nlvls.out").astype(int)   # ulevels (top wet)
    nl = np.loadtxt(f"{sub}/nlvls.out").astype(int)          # nlevels (bottom+1)
    cav = ul > 1
    fu = f"{E}/outdata/fesom/unod.fesom.{y}.nc"
    fv = f"{E}/outdata/fesom/vnod.fesom.{y}.nc"
    if not (os.path.exists(fu) and os.path.exists(fv)):
        return np.nan
    du, dv = nc.Dataset(fu), nc.Dataset(fv)
    u = np.array(du["unod"][:], dtype=float); v = np.array(dv["vnod"][:], dtype=float)
    du.close(); dv.close()
    u[np.abs(u) > FILL] = np.nan; v[np.abs(v) > FILL] = np.nan
    uann = np.nanmean(u, axis=0); vann = np.nanmean(v, axis=0)   # (nod2, nz)
    if uann.shape[0] != len(ul):
        return np.nan
    nz = uann.shape[1]
    ki = np.arange(nz)[None, :]
    wet = (ki >= (ul-1)[:, None]) & (ki <= (nl-2)[:, None])       # (nod2, nz) wet layers
    wu = wet & np.isfinite(uann); wv = wet & np.isfinite(vann)
    with np.errstate(invalid="ignore"):
        ubar = np.nansum(np.where(wu, uann, 0), axis=1) / np.maximum(wu.sum(axis=1), 1)
        vbar = np.nansum(np.where(wv, vann, 0), axis=1) / np.maximum(wv.sum(axis=1), 1)
    spd = np.sqrt(ubar**2 + vbar**2) * 100.0                      # cm/s, barotropic
    na = node_areas(sub)
    w = na[cav]; s = spd[cav]
    m = np.isfinite(s) & (w > 0)
    return float(np.nansum(s[m]*w[m]) / np.nansum(w[m])) if m.any() else np.nan

series = {}
for E, meta in RUNS.items():
    ser = [cavity_mean_speed(f"{BASE}/{E}", y) for y in YEARS]
    series[E] = ser
    print(f"{meta['label']:22s} " + " ".join(f"{v:5.1f}" if np.isfinite(v) else "  nan" for v in ser))

fig, ax = plt.subplots(figsize=(9, 5), facecolor="white")
for E, meta in RUNS.items():
    ax.plot(YEARS, series[E], "-o", color=meta["color"], lw=2.2, ms=6, label=meta["label"])
ax.set_xlabel("year"); ax.set_ylabel("cavity-mean barotropic speed  [cm s$^{-1}$]")
ax.set_title("Sub-ice-shelf cavity current: zstar vs linfs\n(same coupled moving-cavity setup, differing only in which_ale)")
ax.grid(alpha=.3); ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.set_ylim(bottom=0)
fig.tight_layout()
for out in ("/home/a/a270092/esm_tools/movcav4_crash_plots/cavity_velocity_zstar34_vs_linfs37.png",
            "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/report/cavity_velocity_zstar34_vs_linfs37.png"):
    fig.savefig(out, dpi=150)
    print("saved", out)
