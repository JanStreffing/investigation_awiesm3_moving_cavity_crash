#!/usr/bin/env python3
"""Is the Southern-Hemisphere summer melt-back still anomalous?

The report's open item (sec:seaicemelt) was measured in the movcav16/mc25 era:
our runs kept ~7k SH nodes through summer and retreated only ~6x, where the
v3.4 CMIP7 reference retreated ~200x. It was marked "resolved in movcav37" on
the strength of one run. Re-measure on the long records now available:

  repro_movcav45  57 coupled years, correct sub-shelf melt sign throughout
  fgdbg02         24 coupled years (GM fix; wrong melt sign -- ocean side only,
                  so its sea ice is still usable for this question)
  v3.4 reference  AWIESM720_CMIP7_SPINUP_TCO95_CORE3

Metric: SH (lat < -50) sea-ice AREA = sum(a_ice * node_area), monthly
climatology over the last N years of each run; retreat ratio = max/min.
"""
import glob, os
import numpy as np
import netCDF4 as nc

MAX = "/work/ab0246/a270092/input/fesom2/core3"
RUNS = [
    ("repro_movcav45 (57 yr, sign-fixed)",
     "/work/ab0246/a270234/runtime/awiesm3-develop-is/repro_movcav45", None),
    ("fgdbg02 (24 yr)",
     "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02", None),
    ("v3.4 CMIP7 reference",
     "/work/bb1469/a270270/runtime/awiesm3-v3.4/AWIESM720_CMIP7_SPINUP_TCO95_CORE3", None),
]
NLAST = 5


def node_areas(mesh):
    """Spherical-excess triangle areas. NB: a naive equirectangular formula
    explodes on the ~1000 elements spanning +/-180 deg (they get lon spans of
    ~2*pi instead of ~0), inflating the global total from 361 to 860e6 km^2."""
    ao = np.loadtxt(f"{mesh}/nod2d.out", skiprows=1)
    lon, lat = np.radians(ao[:, 1]), np.radians(ao[:, 2])
    e = np.loadtxt(f"{mesh}/elem2d.out", skiprows=1).astype(int) - 1
    R = 6371e3
    # unit vectors -> spherical excess (Girard); wrap-free by construction
    v = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=1)
    a, b, c = v[e[:, 0]], v[e[:, 1]], v[e[:, 2]]
    num = np.abs(np.einsum('ij,ij->i', a, np.cross(b, c)))
    den = (1.0 + np.einsum('ij,ij->i', a, b)
             + np.einsum('ij,ij->i', b, c)
             + np.einsum('ij,ij->i', c, a))
    ar = 2.0 * np.arctan2(num, den) * R * R
    na = np.zeros(len(ao))
    for k in range(3):
        np.add.at(na, e[:, k], ar / 3.0)
    return na, ao[:, 2]


MESHCACHE = {}
base_area, base_lat = node_areas(MAX)
print(f"base mesh {len(base_area)} nodes, SH(<-50) nodes {int((base_lat<-50).sum())}\n")
print(f"{'run':38s} {'years':>12s} {'AREAmax':>8s} {'AREAmin':>8s} | {'EXTmax':>8s} {'EXTmin':>8s}")

for name, root, _ in RUNS:
    fs = sorted(glob.glob(f"{root}/outdata/fesom/a_ice.fesom.*.nc"))
    if not fs:
        print(f"{name:38s}  (no a_ice output)")
        continue
    fs = fs[-NLAST:]
    clim = np.zeros(12)
    climE = np.zeros(12)
    nodes_min = None
    used = []
    for f in fs:
        y = int(f.split("a_ice.fesom.")[1][:4])
        d = nc.Dataset(f)
        a = np.array(d.variables["a_ice"][:])
        d.close()
        a[np.abs(a) > 1e29] = np.nan
        if a.shape[0] != 12:
            continue
        n = a.shape[1]
        # runs on a submesh are a prefix-compatible subset only via map_nod; if the
        # node count differs from the base mesh, fall back to that run's own mesh
        if n == len(base_area):
            na, la = base_area, base_lat
        else:
            # every cycle has its own submesh: pick the one with matching node count
            key = (root, n)
            if key not in MESHCACHE:
                MESHCACHE[key] = None
                for sub in sorted(glob.glob(f"{root}/couple/submesh_*")):
                    try:
                        with open(f"{sub}/nod2d.out") as fh:
                            if int(fh.readline().split()[0]) == n:
                                MESHCACHE[key] = node_areas(sub)
                                break
                    except Exception:
                        pass
            if MESHCACHE[key] is None:
                continue
            na, la = MESHCACHE[key]
        sh = (la < -50)
        aa = np.where(np.isfinite(a[:, sh]), a[:, sh], 0)
        areas = np.nansum(aa * na[sh], axis=1)            # concentration-weighted AREA
        exts = np.nansum((aa > 0.15) * na[sh], axis=1)    # EXTENT (>15% convention)
        clim += areas / 1e12
        climE += exts / 1e12
        used.append(y)
        cnt = np.nansum(a[:, sh] > 0.15, axis=1)
        nodes_min = cnt.min() if nodes_min is None else min(nodes_min, cnt.min())
    if not used:
        print(f"{name:38s}  (no usable years)")
        continue
    clim /= len(used); climE /= len(used)
    mx, mn = clim.max(), clim.min()
    mxE, mnE = climE.max(), climE.min()
    print(f"{name:38s} {min(used)}-{max(used):>5d} {mx:8.2f} {mn:8.2f} | {mxE:8.2f} {mnE:8.2f}")

print(f"{'OBSERVED (NSIDC, approx)':38s} {'':>12s} {15.0:8.1f} {2.0:8.1f} | {18.5:8.1f} {3.0:8.1f}")
print("\nunits 10^6 km^2. AREA = concentration-weighted; EXTENT = cells with a_ice>0.15.")
