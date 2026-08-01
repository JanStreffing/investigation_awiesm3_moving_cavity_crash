"""Pin down the guard predicate and measure its false-positive rate.

Candidate: flag an element where exactly one of u,v is exactly 0.0 at a level
inside the wet range.  Two variants of "inside":
  STRICT  levels ul .. nl-2   (excludes the bottom wet level)
  FULL    levels ul .. nl-1   (the whole wet column)
39030's zero is at its bottom wet level, so STRICT catches only 2 of 3.

A guard is only useful if it never fires on a healthy state.  Run both variants
over every saved restart of the campaign plus the crash step.
"""
import os, glob, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
B = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02"

meshes = {}
for y in ("1901", "1902", "1903", "1904", "1905", "1906"):
    M = f"{CO}/submesh_{y}-12-31T00:00:00"
    meshes[np.loadtxt(f"{M}/cavity_elvls.out").astype(int).size] = (
        y,
        np.loadtxt(f"{M}/cavity_elvls.out").astype(int),
        np.loadtxt(f"{M}/elvls.out").astype(int),
    )
M0 = "/work/ab0246/a270092/input//fesom2/pism_cavity_ini/submesh/"
meshes[np.loadtxt(f"{M0}/cavity_elvls.out").astype(int).size] = (
    "base", np.loadtxt(f"{M0}/cavity_elvls.out").astype(int),
    np.loadtxt(f"{M0}/elvls.out").astype(int))


def test(u, v, ule, nle, full):
    """u,v shaped (elem, nlev). Returns indices where exactly one component is 0."""
    nlev = u.shape[1]
    lv = np.arange(1, nlev + 1)[None, :]
    hi = nle[:, None] - (1 if full else 2)
    interior = (lv >= ule[:, None]) & (lv <= hi)
    zu = (u == 0) & interior
    zv = (v == 0) & interior
    return np.where((zu & ~zv).any(axis=1))[0], np.where((zv & ~zu).any(axis=1))[0]


print(f"{'state':52s} {'STRICT u-only':>14} {'FULL u-only':>12} {'FULL v-only':>12}")
states = []
for d in sorted(glob.glob(f"{B}/run_awiesm3_*/work/fesom.*.oce.restart")):
    leg = d.split("run_awiesm3_")[1].split("/")[0][:4]
    yr = os.path.basename(d).split(".")[1]
    states.append((f"leg {leg} restart {yr}"
                   f" ({'incoming' if int(yr) < int(leg) else 'native  '})", d, None))
for label, d, _ in states:
    with nc.Dataset(f"{d}/u.nc") as f:
        u = np.squeeze(np.asarray(f.variables["u"][:])).T
    with nc.Dataset(f"{d}/v.nc") as f:
        v = np.squeeze(np.asarray(f.variables["v"][:])).T
    if u.shape[0] not in meshes:
        print(f"{label:52s}  (no mesh for elem2D={u.shape[0]})")
        continue
    _, ule, nle = meshes[u.shape[0]]
    s, _ = test(u, v, ule, nle, False)
    fu, fv = test(u, v, ule, nle, True)
    print(f"{label:52s} {s.size:14d} {fu.size:12d} {fv.size:12d}")

RUN = f"{B}/run_awiesm3_19060101-19061231/work"
ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")
u = np.asarray(ds.variables["u"][0, :, :])
v = np.asarray(ds.variables["v"][0, :, :])
_, ule, nle = meshes[u.shape[0]]
s, _ = test(u, v, ule, nle, False)
fu, fv = test(u, v, ule, nle, True)
print(f"{'CRASH STEP 1 (blowup dump)':52s} {s.size:14d} {fu.size:12d} {fv.size:12d}")
print(f"    STRICT hits: {list(s+1)}")
print(f"    FULL   hits: {list(fu+1)}")
