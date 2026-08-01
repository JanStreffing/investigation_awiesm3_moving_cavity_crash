"""Are the crash nodes 'orphan-level' nodes - open at the top, but with every
surrounding element starting below?

Node B's W is exactly 0.0 at levels 1-5, which vert_vel_ale can only produce if NO
edge contribution ever touched those levels, i.e. every element at that node has
ulevels(elem) >= 6.  The submesh's cavity_elvls.out says node B's elements are
[7,1,1,1].  FESOM does not read that file for the element upper level; it derives
ulevels(elem) from the nodes.  Test both conventions.
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

emax = cavn[elem].max(axis=1)
emin = cavn[elem].min(axis=1)
print(f"elements where cavity_elvls.out != max(node ulevels): {int((cave != emax).sum())}")
print(f"elements where cavity_elvls.out != min(node ulevels): {int((cave != emin).sum())}")
mm = np.where(cave != emax)[0]
if mm.size:
    print(f"  mismatch set: cave range {cave[mm].min()}..{cave[mm].max()}, "
          f"max(node) range {emax[mm].min()}..{emax[mm].max()}")

# For each node: the minimum element upper level among its elements, under both
# conventions.  A node is 'orphaned' over [ulevels_nod2D, min_elem_ul) - levels
# that belong to its control volume but receive no element contribution at all.
minel_file = np.full(N, 99, int)
minel_max = np.full(N, 99, int)
nel = np.zeros(N, int)
for e in range(elem.shape[0]):
    for k in elem[e]:
        nel[k] += 1
        if cave[e] < minel_file[k]:
            minel_file[k] = cave[e]
        if emax[e] < minel_max[k]:
            minel_max[k] = emax[e]

for name, minel in (("cavity_elvls.out", minel_file), ("max(node ulevels)", minel_max)):
    orph = np.where((minel > cavn) & (nel > 0))[0]
    print(f"\n convention = {name}: nodes with orphan levels [ul, min_elem_ul): {orph.size}")
    if orph.size and orph.size < 40:
        for n in orph:
            print(f"    node {n+1:7d} lon {lon[n]:8.3f} lat {lat[n]:8.3f} ul_node={cavn[n]} "
                  f"min_elem_ul={minel[n]} nl={nln[n]} nelem={nel[n]} -> orphan levels "
                  f"{list(range(cavn[n], minel[n]))}")
    elif orph.size:
        print(f"    (too many to list) lat range {lat[orph].min():.1f}..{lat[orph].max():.1f}")

# cross-check against the W==0 set found in the blowup
ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")
crash = []
for LO, LA in ((-100.382545, -72.8639145), (-100.112221, -72.7109299),
               (-100.564, -72.592), (-100.719, -72.739), (-100.019, -72.556)):
    crash.append(int(np.argmin((lon - LO) ** 2 + (lat - LA) ** 2)))
print("\n=== the five W-hole nodes ===")
for i in crash:
    w = np.asarray(ds.variables["w"][0, i, :])
    et = np.where(np.any(elem == i, axis=1))[0]
    zeros = [k + 1 for k in range(cavn[i] - 1, nln[i] - 1) if w[k] == 0]
    print(f" node {i+1:7d} lon {lon[i]:8.3f} lat {lat[i]:8.3f} ul={cavn[i]} nl={nln[i]}")
    print(f"   elements       {list(et + 1)}")
    print(f"   cave(file)     {list(cave[et])}")
    print(f"   max(node ul)   {list(emax[et])}")
    print(f"   min over elems: file={minel_file[i]}  max-convention={minel_max[i]}")
    print(f"   W==0 at levels {zeros}")
