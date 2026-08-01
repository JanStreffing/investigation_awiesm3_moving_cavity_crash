"""Is eta==0 under the cavity a normal FESOM state or a remap artefact?

u_rhs on the four newly-iced elements at the crash node is -1.23 m/s, depth-uniform.
dt=1200 s -> accel 1.0e-3 m/s^2 -> g*grad(eta) with grad(eta) ~ 1.05e-4, i.e. a
~1.3 m jump over a ~12 km edge.  The crash node sits at eta=-1.554 and its cavity
neighbours at eta=0.0000 exactly.  Test whether "exactly zero under ice" holds
cavity-wide, and how it compares to the open ocean next to it.
"""
import os, numpy as np
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import netCDF4 as nc

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
RUN = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work"

def rd(path, var):
    with nc.Dataset(path) as ds:
        return np.squeeze(np.asarray(ds.variables[var][:]))

for yr, rst in (("1906", f"{RUN}/fesom.1905.oce.restart"),):
    M = f"{CO}/submesh_{yr}-12-31T00:00:00"
    cavn = np.loadtxt(f"{M}/cavity_nlvls.out").astype(int)
    nod = np.loadtxt(f"{M}/nod2d.out", skiprows=1)
    ssh = rd(f"{rst}/ssh.nc", "ssh")
    cav = cavn > 1
    opn = ~cav
    print(f"\n=== {yr} restart {os.path.basename(rst)} ===")
    print(f" nodes {ssh.size}  cavity {cav.sum()}  open {opn.sum()}")
    print(f" cavity eta: exactly 0 -> {int(np.sum(ssh[cav]==0))}/{cav.sum()}"
          f"   min={ssh[cav].min():.6g} max={ssh[cav].max():.6g} "
          f"mean={ssh[cav].mean():.6g} std={ssh[cav].std():.6g}")
    print(f" open   eta: exactly 0 -> {int(np.sum(ssh[opn]==0))}/{opn.sum()}"
          f"   min={ssh[opn].min():.6g} max={ssh[opn].max():.6g} mean={ssh[opn].mean():.6g}")
    # open nodes south of 60S
    sh = opn & (nod[:, 2] < -60)
    print(f" open eta south of 60S: n={sh.sum()} min={ssh[sh].min():.4g} "
          f"max={ssh[sh].max():.4g} mean={ssh[sh].mean():.4g}")
    # the jump across the cavity front: for each open node adjacent to a cavity node
    elem = np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1
    cave = np.loadtxt(f"{M}/cavity_elvls.out").astype(int)
    # front elements = mixed cavity/open nodes
    mix = (cav[elem].sum(axis=1) > 0) & (cav[elem].sum(axis=1) < 3)
    jump = np.array([ssh[e].max() - ssh[e].min() for e in elem[mix]])
    print(f" cavity-front elements (mixed nodes): n={mix.sum()}  "
          f"|eta jump| median={np.median(jump):.4f} p99={np.percentile(jump,99):.4f} max={jump.max():.4f}")
    # for those front elements, implied g*grad(eta)*dt with 12 km scale
    print(f"   -> implied |UV_rhs| = g*(jump/12km)*1200s : median={9.81*np.median(jump)/12000*1200:.3f} "
          f"max={9.81*jump.max()/12000*1200:.3f} m/s")

# --- how many of those front elements are also ul=7-style pinched (ul>1 but with an open node)?
M = f"{CO}/submesh_1906-12-31T00:00:00"
cavn = np.loadtxt(f"{M}/cavity_nlvls.out").astype(int)
cave = np.loadtxt(f"{M}/cavity_elvls.out").astype(int)
elem = np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1
ssh = rd(f"{RUN}/fesom.1905.oce.restart/ssh.nc", "ssh")
has_open_node = (cavn[elem] == 1).any(axis=1)
print("\n=== elements with ul>1 that still contain an OPEN (ul=1) node ===")
sel = (cave > 1) & has_open_node
print(f" count = {int(sel.sum())} of {elem.shape[0]}")
jump = np.array([ssh[e].max() - ssh[e].min() for e in elem[sel]])
print(f" eta jump across them: median={np.median(jump):.4f} p99={np.percentile(jump,99):.4f} max={jump.max():.4f}")
print(f" n with jump > 1.0 m : {int(np.sum(jump>1.0))}")
