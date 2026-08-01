"""Read mesh_areas.bin from the instrumented fesom_meshdiag and test the crash columns.

vert_vel_ale does  Wvel(nz,n) = cumulative_divergence(nz,n) / area(nz,n).
An exact 0.0 in W inside the wet column cannot come out of that sum by chance, so
either area(nz,n) is degenerate at those levels or the divergence really is zero.
The tracer volume uses areasvol(nz,n); if that is degenerate too, the near-singular
vertical tridiagonal (temp ~1e250 decaying 1e3/level) has the same origin.

Layout written by dump_areas (Fortran stream, -r8 -i4):
  int32 nl, nod2D, elem2D
  float64 area(nod2D)      x nl
  float64 areasvol(nod2D)  x nl
  int32   ulevels_nod2D(nod2D), nlevels_nod2D(nod2D)
  int32   ulevels(elem2D),      nlevels(elem2D)
  float64 elem_area(elem2D)
"""
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

BIN = "/work/ab0246/a270092/tmp/meshdiag_run/out/mesh_areas.bin"
CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
RUN = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work"
M = f"{CO}/submesh_1906-12-31T00:00:00"

with open(BIN, "rb") as f:
    buf = f.read()
off = 0
nl, nn, ne = np.frombuffer(buf, np.int32, 3, off); off += 12
nl, nn, ne = int(nl), int(nn), int(ne)
print(f"nl={nl} nod2D={nn} elem2D={ne}")
area = np.frombuffer(buf, np.float64, nl * nn, off).reshape(nl, nn); off += 8 * nl * nn
areasvol = np.frombuffer(buf, np.float64, nl * nn, off).reshape(nl, nn); off += 8 * nl * nn
uln = np.frombuffer(buf, np.int32, nn, off).copy(); off += 4 * nn
nln = np.frombuffer(buf, np.int32, nn, off).copy(); off += 4 * nn
ule = np.frombuffer(buf, np.int32, ne, off).copy(); off += 4 * ne
nle = np.frombuffer(buf, np.int32, ne, off).copy(); off += 4 * ne
earea = np.frombuffer(buf, np.float64, ne, off).copy(); off += 8 * ne
print(f"consumed {off} of {len(buf)} bytes")

nod = np.loadtxt(f"{M}/nod2d.out", skiprows=1)
lon, lat = nod[:, 1], nod[:, 2]
cavn = np.loadtxt(f"{M}/cavity_nlvls.out").astype(int)
print(f"ulevels_nod2D vs cavity_nlvls.out: identical for "
      f"{int((uln == cavn).sum())}/{nn}")

ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")

print("\n=== population: degenerate area / areasvol inside the wet column ===")
for name, A in (("area", area), ("areasvol", areasvol)):
    bad_zero, bad_small = [], []
    for i in range(nn):
        u, b = uln[i], nln[i]
        if b - u < 2:
            continue
        seg = A[u - 1:b - 1, i]
        if np.any(seg == 0):
            bad_zero.append(i)
        elif np.any(seg < 1.0):
            bad_small.append(i)
    print(f" {name:9s} == 0 inside column: {len(bad_zero):6d}   "
          f"< 1 m^2 inside column: {len(bad_small):6d}")
    if bad_zero[:10]:
        print(f"   e.g. {[b+1 for b in bad_zero[:10]]}")

print("\n=== the five W-hole columns ===")
pts = [(-100.382545, -72.8639145), (-100.112221, -72.7109299),
       (-100.564, -72.592), (-100.719, -72.739), (-100.019, -72.556)]
np.set_printoptions(linewidth=250)
for LO, LA in pts:
    i = int(np.argmin((lon - LO) ** 2 + (lat - LA) ** 2))
    u, b = int(uln[i]), int(nln[i])
    w = np.asarray(ds.variables["w"][0, i, :b])
    zeros = [k + 1 for k in range(u - 1, b - 1) if w[k] == 0]
    print(f"\n node {i+1} lon {lon[i]:8.3f} lat {lat[i]:8.3f} ul={u} nl={b}  W==0 at {zeros}")
    print(f"   area     [1..{b}] = {np.array2string(area[:b, i], precision=4, formatter={'float_kind':lambda x: f'{x:10.4g}'})}")
    print(f"   areasvol [1..{b}] = {np.array2string(areasvol[:b, i], precision=4, formatter={'float_kind':lambda x: f'{x:10.4g}'})}")
    print(f"   W        [1..{b}] = {np.array2string(w, precision=3)}")

print("\n=== healthy controls in the same region ===")
for LO, LA in ((-101.428, -72.911), (-101.012, -72.238), (-100.3, -73.3)):
    i = int(np.argmin((lon - LO) ** 2 + (lat - LA) ** 2))
    u, b = int(uln[i]), int(nln[i])
    print(f" node {i+1} lon {lon[i]:8.3f} lat {lat[i]:8.3f} ul={u} nl={b}")
    print(f"   area     = {np.array2string(area[:b, i], formatter={'float_kind':lambda x: f'{x:10.4g}'})}")
    print(f"   areasvol = {np.array2string(areasvol[:b, i], formatter={'float_kind':lambda x: f'{x:10.4g}'})}")
