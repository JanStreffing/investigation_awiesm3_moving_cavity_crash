"""Census over ALL 1311 cavity-front elements: is the crash element an outlier?

Established so far:
  * eta == 0 at every cavity node (FESOM's own convention, present in native
    restarts too), open ocean south of 60S sits at ~-1.58 m -> every one of the
    1311 mixed cavity/open elements carries a ~1.59 m eta step.
  * -g*grad(eta)*dt over a 12 km edge = ~1.5 m/s; measured u_rhs on the four
    newly-iced elements at the crash node is -0.69 .. -1.25 m/s, i.e. the step is
    only ~20% cancelled by the baroclinic pgf.

Question: do all 1311 front elements carry that, or only some?  And where do the
nodes with |ssh_rhs| > 1e5 sit relative to the front?
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
nle = nln[elem].min(axis=1)
lon, lat = nod[:, 1], nod[:, 2]

cav = cavn > 1
ncav = cav[elem].sum(axis=1)
front = (ncav > 0) & (ncav < 3)                 # mixed element
deepcav = ncav == 3
openel = ncav == 0
print(f"elements: open {openel.sum()}  front(mixed) {front.sum()}  full-cavity {deepcav.sum()}")

ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")
sshr = np.asarray(ds.variables["ssh_rhs"][0, :])
deta = np.asarray(ds.variables["d_eta"][0, :])

print("\n=== where are the big-|ssh_rhs| nodes? ===")
big = np.abs(sshr) > 1e5
print(f" n(|ssh_rhs|>1e5) = {big.sum()}")
# is a node 'front-adjacent'? touches a mixed element
front_node = np.zeros(cavn.size, bool)
front_node[elem[front].ravel()] = True
print(f"   of which front-adjacent: {int((big & front_node).sum())} "
      f"({100*(big & front_node).sum()/max(big.sum(),1):.1f}%)")
print(f"   cavity nodes: {int((big & cav).sum())}   open nodes: {int((big & ~cav).sum())}")
print(f"   lat range: {lat[big].min():.1f} .. {lat[big].max():.1f};  "
      f"n south of 60S: {int((big & (lat<-60)).sum())}")
for thr in (1e4, 1e5, 1e6, 1e7):
    print(f"   |ssh_rhs|>{thr:8.0e}: {int((np.abs(sshr)>thr).sum()):6d}  "
          f"front-adjacent {int(((np.abs(sshr)>thr)&front_node).sum()):6d}")

print("\n=== u_rhs depth-integral per element class (blowup, post-step-1) ===")
# read u_rhs/v_rhs/helem in chunks; only need front + a control sample
def col_int(v_name, idx):
    out = np.zeros(idx.size)
    var = ds.variables[v_name]
    he = ds.variables["helem"]
    for j, e in enumerate(idx):
        a = np.asarray(var[0, e, :]); h = np.asarray(he[0, e, :])
        out[j] = np.sum(a * h)
    return out

rng = np.random.default_rng(0)
ctrl = rng.choice(np.where(openel)[0], 400, replace=False)
fidx = np.where(front)[0]
didx = rng.choice(np.where(deepcav)[0], 400, replace=False)

for name, idx in (("front(mixed)", fidx), ("full-cavity", didx), ("open ctrl", ctrl)):
    ur = np.asarray([np.asarray(ds.variables["u_rhs"][0, e, :]) for e in idx])
    vr = np.asarray([np.asarray(ds.variables["v_rhs"][0, e, :]) for e in idx])
    he = np.asarray([np.asarray(ds.variables["helem"][0, e, :]) for e in idx])
    mag = np.hypot(ur, vr)
    wet = he > 0
    peak = np.array([mag[i][wet[i]].max() if wet[i].any() else 0 for i in range(idx.size)])
    trans = np.abs(np.sum(np.hypot(ur, vr) * he, axis=1))
    print(f" {name:14s} n={idx.size:5d}  max|u_rhs| over wet levels: "
          f"p50={np.percentile(peak,50):8.4f} p90={np.percentile(peak,90):8.4f} "
          f"p99={np.percentile(peak,99):8.4f} max={peak.max():8.4f}")
    print(f" {'':14s}          |int(u_rhs h)|: p50={np.percentile(trans,50):10.3f} "
          f"p99={np.percentile(trans,99):10.3f} max={trans.max():10.3f}")

# the four crash elements, for reference
i0 = int(np.argmin((lon + 100.382545) ** 2 + (lat + 72.8639145) ** 2))
et = np.where(np.any(elem == i0, axis=1))[0]
print("\n=== crash-node elements in that distribution ===")
for e in et:
    ur = np.asarray(ds.variables["u_rhs"][0, e, :]); vr = np.asarray(ds.variables["v_rhs"][0, e, :])
    he = np.asarray(ds.variables["helem"][0, e, :])
    wet = he > 0
    mg = np.hypot(ur, vr)
    print(f"  elem {e+1} ul={cave[e]} front={front[e]} deepcav={deepcav[e]}  "
          f"max|u_rhs|(wet)={mg[wet].max() if wet.any() else 0:.4f}  "
          f"|int(u_rhs h)|={abs(np.sum(np.hypot(ur,vr)*he)):.2f}")

# which of the front elements have an OPEN node with few open elements (the pinch)?
print("\n=== the pinch, quantified over all open nodes ===")
nelem_at = np.zeros(cavn.size, int)
nopen_at = np.zeros(cavn.size, int)
for e in range(elem.shape[0]):
    o = cave[e] == 1
    for k in elem[e]:
        nelem_at[k] += 1
        if o:
            nopen_at[k] += 1
om = (cavn == 1) & (nelem_at >= 4)
pinched = om & (nopen_at <= 2)
print(f" open nodes with >=4 elements: {int(om.sum())};  pinched (<=2 open elems): {int(pinched.sum())}")
pidx = np.where(pinched)[0]
print(f" pinched nodes: |ssh_rhs| p50={np.percentile(np.abs(sshr[pidx]),50):.4g} "
      f"p90={np.percentile(np.abs(sshr[pidx]),90):.4g} max={np.abs(sshr[pidx]).max():.4g}")
print(f" crash node A ssh_rhs={sshr[i0]:.4g}  rank among pinched: "
      f"{int((np.abs(sshr[pidx])>abs(sshr[i0])).sum())+1} of {pidx.size}")
ordr = pidx[np.argsort(-np.abs(sshr[pidx]))][:15]
print(" top-15 pinched nodes by |ssh_rhs|:")
for n in ordr:
    print(f"   node {n+1:7d} lon {lon[n]:8.3f} lat {lat[n]:8.3f} openelems {nopen_at[n]}/{nelem_at[n]} "
          f"ssh_rhs={sshr[n]:14.4g} d_eta={deta[n]:10.5f}")
