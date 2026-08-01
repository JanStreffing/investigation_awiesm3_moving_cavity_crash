"""Is the u-only interior zero unique in space AND time?

The fesom_meshdiag dump established that cavity_elvls.out == model ulevels(elem) and
elvls.out == model nlevels(elem), exactly, for all 404 607 elements.  So the correct
level bounds are available for every submesh without rerunning the utility.

Two distinct populations to keep apart:
  BOTH  - u and v zero at identical levels (retreating-cavity remap residue; benign,
          the model has carried these for years)
  U-ONLY- u zero where v is not (the step-1 signature at elements 39028-39030)

Scan every restart the campaign wrote: each leg's incoming remapped restart and its
native end-of-year restart.
"""
import os, glob, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
B = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02"

meshes = {}
for y in ("1901", "1902", "1903", "1904", "1905", "1906"):
    M = f"{CO}/submesh_{y}-12-31T00:00:00"
    ule = np.loadtxt(f"{M}/cavity_elvls.out").astype(int)
    nle = np.loadtxt(f"{M}/elvls.out").astype(int)
    nod = np.loadtxt(f"{M}/nod2d.out", skiprows=1)
    el = np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1
    meshes[ule.size] = (y, ule, nle, nod, el)
print("submesh element counts:", {k: v[0] for k, v in meshes.items()})


def scan(u, v, ule, nle):
    both, uonly, vonly = [], [], []
    for i in range(ule.size):
        lo, hi = ule[i] - 1, nle[i] - 1
        if hi - lo < 2:
            continue
        zu = set(np.where(u[lo:hi, i] == 0)[0].tolist())
        zv = set(np.where(v[lo:hi, i] == 0)[0].tolist())
        if not zu and not zv:
            continue
        if zu == zv:
            both.append(i)
        elif zu - zv:
            uonly.append((i, sorted(k + lo + 1 for k in zu - zv)))
        else:
            vonly.append(i)
    return both, uonly, vonly


rows = []
for d in sorted(glob.glob(f"{B}/run_awiesm3_*/work/fesom.*.oce.restart")):
    leg = d.split("run_awiesm3_")[1].split("/")[0]
    yr = os.path.basename(d).split(".")[1]
    kind = "incoming(remapped)" if int(yr) < int(leg[:4]) else "native end-of-leg"
    try:
        with nc.Dataset(f"{d}/u.nc") as f:
            u = np.squeeze(np.asarray(f.variables["u"][:]))
        with nc.Dataset(f"{d}/v.nc") as f:
            v = np.squeeze(np.asarray(f.variables["v"][:]))
    except Exception as e:
        print(f"{leg} {yr}: {e}")
        continue
    ne = u.shape[1]
    if ne not in meshes:
        print(f"leg {leg} {yr}: no mesh with elem2D={ne}")
        continue
    mname, ule, nle, nod, el = meshes[ne]
    both, uonly, vonly = scan(u, v, ule, nle)
    rows.append((leg, yr, kind, mname, len(both), len(uonly), len(vonly), uonly))
    print(f" leg {leg}  {kind:19s} restart {yr}  (mesh {mname})   "
          f"BOTH={len(both):4d}   U-ONLY={len(uonly):3d}   V-ONLY={len(vonly):3d}")
    for i, lv in uonly[:10]:
        print(f"      u-only elem {i+1:7d} lon {nod[el[i],1].mean():8.3f} "
              f"lat {nod[el[i],2].mean():8.3f} ul={ule[i]} nl={nle[i]} levels {lv}")

print("\n=== and the crash step itself (blowup, post-step-1, mesh 1906) ===")
RUN = f"{B}/run_awiesm3_19060101-19061231/work"
ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")
mname, ule, nle, nod, el = meshes[404607]
ub = np.asarray(ds.variables["u"][0, :, :]).T
vb = np.asarray(ds.variables["v"][0, :, :]).T
both, uonly, vonly = scan(ub, vb, ule, nle)
print(f" blowup step 1: BOTH={len(both)}  U-ONLY={len(uonly)}  V-ONLY={len(vonly)}")
for i, lv in uonly:
    print(f"      u-only elem {i+1:7d} lon {nod[el[i],1].mean():8.3f} "
          f"lat {nod[el[i],2].mean():8.3f} ul={ule[i]} nl={nle[i]} levels {lv}")
