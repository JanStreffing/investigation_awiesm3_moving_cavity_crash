"""What is unique about the two nodes that exploded?

ssh_rhs is NOT it: 506 nodes carry |ssh_rhs| > 1e7 and the crash node's 7.1e6
ranks 143rd of 184 pinched nodes.  u_rhs is not it either: max|u_rhs| ~1.0 m/s is
the median of all 1311 cavity-front elements.  d_eta stays ~1e-3 everywhere.

The temperature profile at the crash node decays geometrically away from a single
level (peak 1e250 at nz=5 for node A, nz=4 for node B) - the signature of a
tridiagonal solve, not of advection.  And W has an exact 0.0 sitting inside the
wet column at that same level.  Look for a vertical-column bookkeeping fault:
hnode==0 inside [ulevel, nlevel), W==0 inside the column, etc.
"""
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
RUN = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work"
M = f"{CO}/submesh_1906-12-31T00:00:00"

nod = np.loadtxt(f"{M}/nod2d.out", skiprows=1)
elem = np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1
cavn = np.loadtxt(f"{M}/cavity_nlvls.out").astype(int)
cave = np.loadtxt(f"{M}/cavity_elvls.out").astype(int)
nln = np.loadtxt(f"{M}/nlvls.out").astype(int)
lon, lat = nod[:, 1], nod[:, 2]
N = cavn.size

ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")
rst = f"{RUN}/fesom.1905.oce.restart"

def rrd(f, v):
    with nc.Dataset(f"{rst}/{f}.nc") as d:
        return np.squeeze(np.asarray(d.variables[v][:]))

i0 = int(np.argmin((lon + 100.382545) ** 2 + (lat + 72.8639145) ** 2))
i1 = int(np.argmin((lon + 100.112221) ** 2 + (lat + 72.7109299) ** 2))

hn_all = ds.variables["hnode"]
w_all = ds.variables["w"]
hn_rst = rrd("hnode", "hnode")     # (nz_1, node) as written by the remap

print("=== the two blown-up nodes: column bookkeeping ===")
for tag, i in (("A", i0), ("B", i1)):
    ul, nl = cavn[i], nln[i]
    hn = np.asarray(hn_all[0, i, :])
    w = np.asarray(w_all[0, i, :])
    hr = hn_rst[:, i] if hn_rst.shape[1] == N else hn_rst[i, :]
    print(f"\n node {tag} idx1={i+1} lon {lon[i]:.3f} lat {lat[i]:.3f}  ul={ul} nl={nl}")
    print(f"   hnode(blowup) [1..{nl}] = {np.array2string(hn[:nl], precision=2, max_line_width=250)}")
    print(f"   hnode(restart)[1..{nl}] = {np.array2string(hr[:nl], precision=2, max_line_width=250)}")
    print(f"   W(blowup)     [1..{nl}] = {np.array2string(w[:nl], precision=3, max_line_width=250)}")
    inner = np.arange(ul - 1, nl - 1)
    print(f"   zeros in hnode inside [ul,nl): {list(inner[hn[inner] == 0] + 1)}")
    print(f"   zeros in W     inside [ul,nl): {list(inner[w[inner] == 0] + 1)}")

# --- population statistics: how many nodes have a zero hole in hnode or W? ---
print("\n=== population scan: interior zeros in hnode / W (all nodes) ===")
CH = 20000
bad_h, bad_w = [], []
for s in range(0, N, CH):
    e = min(N, s + CH)
    hn = np.asarray(hn_all[0, s:e, :])
    w = np.asarray(w_all[0, s:e, :])
    for j in range(e - s):
        i = s + j
        ul, nl = cavn[i], nln[i]
        if nl - ul < 2:
            continue
        inner = np.arange(ul - 1, nl - 1)
        if np.any(hn[j, inner] == 0):
            bad_h.append(i)
        if np.any(w[j, inner] == 0):
            bad_w.append(i)
print(f" nodes with hnode==0 inside the wet column : {len(bad_h)}")
print(f" nodes with W==0     inside the wet column : {len(bad_w)}")
for nm, lst in (("hnode-hole", bad_h), ("W-hole", bad_w)):
    if not lst:
        continue
    a = np.array(lst)
    print(f"  {nm}: crash A in list? {i0 in set(lst)}   crash B in list? {i1 in set(lst)}")
    print(f"  {nm}: lat range {lat[a].min():.1f}..{lat[a].max():.1f}; "
          f"cavity {int((cavn[a]>1).sum())} open {int((cavn[a]==1).sum())}")
    for n in a[:25]:
        print(f"     node {n+1:7d} lon {lon[n]:8.3f} lat {lat[n]:8.3f} ul={cavn[n]} nl={nln[n]}")
