"""Interior exact zeros in the velocity field.

area(nz,n) and areasvol(nz,n) came back clean (fesom_meshdiag on the 1906 submesh:
no zeros, nothing below 1 m^2 inside any wet column), so the W==0 holes are not a
control-volume degeneracy.  But W is the cumulative divergence of (u,v)*helem, and
element 39028 at the crash node already showed u = [.0631 .0575 .0081 0. 0. -.0275 ...]
- two exact zeros mid-column.  An exact 0.0 in u is as impossible by chance as one
in W.  Scan every element column for interior zeros in u and v, and see whether the
set explains the five W holes.
"""
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
RUN = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work"
M = f"{CO}/submesh_1906-12-31T00:00:00"

nod = np.loadtxt(f"{M}/nod2d.out", skiprows=1)
elem = np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1
lon, lat = nod[:, 1], nod[:, 2]
ne = elem.shape[0]

# element level bounds straight from the model (fesom_meshdiag dump), not derived:
# min-over-nodes is NOT what FESOM uses and gave 23 780 false positives.
BIN = "/work/ab0246/a270092/tmp/meshdiag_run/out/mesh_areas.bin"
with open(BIN, "rb") as f:
    buf = f.read()
off = 0
_nl, _nn, _ne = (int(x) for x in np.frombuffer(buf, np.int32, 3, off)); off += 12
off += 8 * _nl * _nn * 2                       # area, areasvol
uln = np.frombuffer(buf, np.int32, _nn, off).copy(); off += 4 * _nn
nln = np.frombuffer(buf, np.int32, _nn, off).copy(); off += 4 * _nn
cave = np.frombuffer(buf, np.int32, _ne, off).copy(); off += 4 * _ne
nle = np.frombuffer(buf, np.int32, _ne, off).copy(); off += 4 * _ne
print(f"element ulevels from model: min={cave.min()} max={cave.max()}; "
      f"nlevels min={nle.min()} max={nle.max()}")
print(f"model ulevels(elem) == cavity_elvls.out for "
      f"{int((cave == np.loadtxt(f'{M}/cavity_elvls.out').astype(int)).sum())}/{ne} elements")

ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")

def scan(varname, ul, nl_, n):
    """indices with an exact zero strictly inside [ul, nl-1)"""
    bad = []
    CH = 20000
    var = ds.variables[varname]
    for s in range(0, n, CH):
        e = min(n, s + CH)
        a = np.asarray(var[0, s:e, :])
        for j in range(e - s):
            i = s + j
            lo, hi = ul[i] - 1, nl_[i] - 1
            if hi - lo < 2:
                continue
            if np.any(a[j, lo:hi] == 0):
                bad.append(i)
    return np.array(bad, int)

print("=== interior exact zeros, blowup (post-step-1) ===")
res = {}
for v in ("u", "v", "helem"):
    b = scan(v, cave, nle, ne)
    res[v] = b
    print(f" {v:6s}: {b.size:6d} elements of {ne}")
    if b.size and b.size < 60:
        for i in b:
            c = lon[elem[i]].mean(), lat[elem[i]].mean()
            print(f"    elem {i+1:7d} centroid {c[0]:8.3f} {c[1]:8.3f} ul={cave[i]} nl={nle[i]}")

# same test on the incoming remapped restart, which FESOM had not touched yet
print("\n=== same test on the remapped restart handed to the leg ===")
rst = f"{RUN}/fesom.1905.oce.restart"
class RD:
    def __init__(self, path, var):
        self.d = nc.Dataset(path); self.v = self.d.variables[var]
    def __getitem__(self, k):
        _, sl, _ = k
        return np.squeeze(np.asarray(self.v[:, :, sl])).T
for v in ("u", "v"):
    with nc.Dataset(f"{rst}/{v}.nc") as d:
        a = np.squeeze(np.asarray(d.variables[v][:]))     # (nz, elem)
    bad = []
    for i in range(ne):
        lo, hi = cave[i] - 1, nle[i] - 1
        if hi - lo < 2:
            continue
        if np.any(a[lo:hi, i] == 0):
            bad.append(i)
    print(f" {v:6s}: {len(bad):6d} elements of {ne}")
    res[v + "_rst"] = np.array(bad, int)
    if 0 < len(bad) < 60:
        for i in bad:
            print(f"    elem {i+1:7d} centroid {lon[elem[i]].mean():8.3f} {lat[elem[i]].mean():8.3f} "
                  f"ul={cave[i]} nl={nle[i]}")

# do the holes sit on the five W-hole nodes?
print("\n=== relation to the five W-hole nodes ===")
pts = [(-100.382545, -72.8639145), (-100.112221, -72.7109299),
       (-100.564, -72.592), (-100.719, -72.739), (-100.019, -72.556)]
crash_nodes = [int(np.argmin((lon - a) ** 2 + (lat - b) ** 2)) for a, b in pts]
for key in ("u", "v", "u_rst", "v_rst"):
    b = res.get(key)
    if b is None or b.size == 0:
        continue
    touch = set()
    for i in b:
        for k in elem[i]:
            touch.add(k)
    print(f" {key:6s}: hole elements touch {len(touch)} nodes; "
          f"of the 5 crash nodes: {sum(1 for c in crash_nodes if c in touch)}")
