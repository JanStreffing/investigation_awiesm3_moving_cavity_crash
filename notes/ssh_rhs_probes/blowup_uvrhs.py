"""Which term makes ssh_rhs = 7.1e6 at the crash node?

ssh_rhs(n) = SUM_edges alpha*((UV+UV_rhs).n_hat)*helem      + (1-alpha)*ssh_rhs_old
ssh_rhs_old(n) = SUM_edges (UV.n_hat)*helem                  (recomputed AFTER update_vel)

The log prints ssh_rhs = 7.14e6 but ssh_rhs_old = -208 at the same node, i.e. the
bare-velocity divergence is sane and the (UV+UV_rhs) one is not.  Here we open the
2 GB blowup dump (post-step-1 state, all fields) and look at the crash node's own
six elements: per-level u, v, u_rhs, v_rhs, urhs_AB, helem, plus the depth integral
of each.  Compare against the restart (pre-step) values.
"""
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

RUN = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work"
BU = f"{RUN}/fesom.1906.oce.blowup.nc"
RST = f"{RUN}/fesom.1905.oce.restart"
M = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple/submesh_1906-12-31T00:00:00"
LON0, LAT0 = -100.382545, -72.8639145
LON1, LAT1 = -100.112221, -72.7109299          # the second blown-up node

nod = np.loadtxt(f"{M}/nod2d.out", skiprows=1)
elem = np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1
cavn = np.loadtxt(f"{M}/cavity_nlvls.out").astype(int)   # node upper level
cave = np.loadtxt(f"{M}/cavity_elvls.out").astype(int)   # elem upper level
nln = np.loadtxt(f"{M}/nlvls.out").astype(int)           # node bottom level
mp = np.loadtxt(f"{M}/map_nod.out").astype(int)
lon, lat = nod[:, 1], nod[:, 2]

# elem bottom level = min over nodes
nle = nln[elem].min(axis=1)

i0 = int(np.argmin((lon - LON0) ** 2 + (lat - LAT0) ** 2))
i1 = int(np.argmin((lon - LON1) ** 2 + (lat - LAT1) ** 2))
print(f"crash node A idx1={i0+1} maxid={mp[i0]} ul={cavn[i0]} nl={nln[i0]}")
print(f"crash node B idx1={i1+1} maxid={mp[i1]} ul={cavn[i1]} nl={nln[i1]}")

for tag, i in (("A", i0), ("B", i1)):
    et = np.where(np.any(elem == i, axis=1))[0]
    print(f"\nnode {tag}: {et.size} elements  ul(elem)={cave[et]}  nl(elem)={nle[et]}")

ds = nc.Dataset(BU)
rs = {v: nc.Dataset(f"{RST}/{v}.nc") for v in ("u", "v", "urhs_AB", "ssh_rhs_old", "hnode")}

def bu(v, idx):
    return np.asarray(ds.variables[v][0, idx, :])

print("\n=== ssh_rhs / d_eta at the two nodes (blowup file vs restart ssh_rhs_old) ===")
for tag, i in (("A", i0), ("B", i1)):
    print(f" node {tag}: ssh_rhs={ds.variables['ssh_rhs'][0,i]:.6g} "
          f"ssh_rhs_old(blowup)={ds.variables['ssh_rhs_old'][0,i]:.6g} "
          f"ssh_rhs_old(restart)={rs['ssh_rhs_old'].variables['ssh_rhs_old'][0,i]:.6g} "
          f"d_eta={ds.variables['d_eta'][0,i]:.6g} eta={ds.variables['eta_n'][0,i]:.6g}")

print("\n=== global extremes in the blowup dump ===")
for v in ("ssh_rhs", "d_eta", "eta_n"):
    a = np.asarray(ds.variables[v][0, :])
    print(f" {v:12s} min={np.nanmin(a):.6g} max={np.nanmax(a):.6g} "
          f"argmax|.|={int(np.nanargmax(np.abs(a)))+1} "
          f"n(|.|>1e5)={int(np.sum(np.abs(a)>1e5))}")

print("\n=== per-element columns at crash node A (post-step-1 blowup state) ===")
et = np.where(np.any(elem == i0, axis=1))[0]
for e in et:
    u = bu("u", e); v = bu("v", e); ur = bu("u_rhs", e); vr = bu("v_rhs", e)
    ab = bu("urhs_AB", e); he = bu("helem", e)
    ul, nl = cave[e], nle[e]
    iu = np.trapz  # unused
    intU = np.sum((u[ul-1:nl-1]) * he[ul-1:nl-1])
    intUr = np.sum((ur[ul-1:nl-1]) * he[ul-1:nl-1])
    intV = np.sum((v[ul-1:nl-1]) * he[ul-1:nl-1])
    intVr = np.sum((vr[ul-1:nl-1]) * he[ul-1:nl-1])
    print(f"\n elem {e+1} ul={ul} nl={nl}")
    print(f"   int(u h)={intU:12.4f}  int(u_rhs h)={intUr:14.4f}   ratio={intUr/(intU+1e-30):10.2f}")
    print(f"   int(v h)={intV:12.4f}  int(v_rhs h)={intVr:14.4f}   ratio={intVr/(intV+1e-30):10.2f}")
    print(f"   u    [1:18]={np.array2string(u[:18], precision=4, max_line_width=200)}")
    print(f"   u_rhs[1:18]={np.array2string(ur[:18], precision=4, max_line_width=200)}")
    print(f"   v_rhs[1:18]={np.array2string(vr[:18], precision=4, max_line_width=200)}")
    print(f"   helem[1:18]={np.array2string(he[:18], precision=2, max_line_width=200)}")
    print(f"   urhsAB[1:18]={np.array2string(ab[:18], precision=2, max_line_width=200)}")

print("\n=== how anomalous is u_rhs globally? ===")
# sample the whole u_rhs field level 1 and the depth-integral proxy on a subset
ur1 = np.asarray(ds.variables["u_rhs"][0, :, 0])
vr1 = np.asarray(ds.variables["v_rhs"][0, :, 0])
mag = np.hypot(ur1, vr1)
print(f" |u_rhs| level1: p50={np.nanpercentile(mag,50):.4g} p99.9={np.nanpercentile(mag,99.9):.4g} "
      f"max={np.nanmax(mag):.4g} at elem {int(np.nanargmax(mag))+1}")
print(f" crash-node elems |u_rhs| level1 = {np.array2string(mag[et], precision=5)}")
