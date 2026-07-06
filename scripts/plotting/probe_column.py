#!/usr/bin/env python3
"""Vertical-structure probe at the worst cfl_z nodes in the FESOM blowup.
cfl_z = w*dt/dz -> distinguish thin-layer (dz tiny; structural/cavity) from
huge-w (dynamic/forced). Dumps per-level dz(hnode), w, temp, cfl_z."""
import numpy as np, netCDF4 as nc

MESH = "/work/ab0246/a270092/pism_repro/experiments/movcav4/couple/latest_submesh"
BLOW = "/work/ab0246/a270092/pism_repro/experiments/movcav4/run_awiesm3_19010101-19011231/work/fesom.1901.oce.blowup.nc"

with open(f"{MESH}/nod2d.out") as f:
    nn = int(f.readline().split()[0]); nod = np.loadtxt(f, max_rows=nn)
lon, lat = nod[:, 1], nod[:, 2]
cav = np.loadtxt(f"{MESH}/cavity_nlvls.out", skiprows=1).astype(int)  # cavity lvls/node
nl  = np.loadtxt(f"{MESH}/nlvls.out", skiprows=1).astype(int)          # bottom lvl/node

ds = nc.Dataset(BLOW)
def g(n): return np.squeeze(np.asarray(ds.variables[n][:]))
cflz = g("cfl_z"); hnode = g("hnode"); w = g("w"); temp = g("temp")
hflx = g("heat_flux"); aice = g("a_ice")
print("shapes: cfl_z", cflz.shape, "hnode", hnode.shape, "w", w.shape)
dt = 1200.0

# worst nodes by max cfl_z over depth
cmax = np.nanmax(cflz, axis=1)
worst = np.argsort(cmax)[::-1][:8]
print(f"\nglobal cfl_z max = {np.nanmax(cflz):.1f}\n")
for n in worst:
    print(f"=== node {n}  lon={lon[n]:7.2f} lat={lat[n]:6.2f}  "
          f"cav_nlvls={cav[n]} bottom_nlvls={nl[n]}  "
          f"heat_flux={hflx[n]:8.1f}  a_ice={aice[n]:.2f} ===")
    kmax = min(int(nl[n]) + 1, cflz.shape[1])
    print("  lev |   dz(hnode) |      w(m/s) |   cfl_z |  temp")
    for k in range(kmax):
        dz = hnode[n, k] if k < hnode.shape[1] else np.nan
        wk = w[n, k] if k < w.shape[1] else np.nan
        ck = cflz[n, k] if k < cflz.shape[1] else np.nan
        tk = temp[n, k] if k < temp.shape[1] else np.nan
        flag = ""
        if np.isfinite(ck) and ck > 5: flag += " <<CFL"
        if np.isfinite(dz) and dz < 1.0: flag += " <<THIN"
        print(f"  {k:3d} | {dz:11.4f} | {wk:11.4e} | {ck:7.2f} | {tk:7.2f}{flag}")
    print()

# summary: are high-cfl driven by thin dz or big w?
hi = cmax > 10
if hi.sum():
    # for hi nodes, at the level of max cfl, get dz and w
    kmaxlev = np.nanargmax(np.nan_to_num(cflz[hi]), axis=1)
    idx = np.where(hi)[0]
    dzs = np.array([hnode[i, k] for i, k in zip(idx, kmaxlev)])
    ws  = np.array([abs(w[i, k]) for i, k in zip(idx, kmaxlev)])
    print(f"--- {hi.sum()} nodes with cfl_z>10 ---")
    print(f"  dz at max-cfl level:  min={dzs.min():.3f} med={np.median(dzs):.3f} max={dzs.max():.3f} m")
    print(f"  |w| at max-cfl level: min={ws.min():.2e} med={np.median(ws):.2e} max={ws.max():.2e} m/s")
    print(f"  cavity nodes among them: {(cav[idx]>0).sum()}/{hi.sum()}")
    print(f"  a_ice>0.1 among them:    {(aice[idx]>0.1).sum()}/{hi.sum()}")
