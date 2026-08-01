"""Reconstruct the vertical-velocity column at the crash node from the blowup u,v,helem.

vert_vel_ale accumulates, per level, the flux of (u,v)*helem through the median-dual
cell boundary, then cumulates bottom-up and divides by area(nz,n).  The model's own
W has an exact 0.0 at level 5 of node A (levels 1-5 of node B) - the only 5 such
columns in the mesh, and zero such columns in any of the six healthy restarts.

Here we redo the FVM flux sum with our own geometry: for each element e containing
node n, the contribution is the flux through the two segments joining the midpoints
of the two edges of e at n to e's centroid.  If the true cumulative sum passes
smoothly through level 5, the model's 0.0 is a fault; if it genuinely reverses sign
there, the hole is a coincidence.
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

R = 6371000.0
rad = np.pi / 180

ds = nc.Dataset(f"{RUN}/fesom.1906.oce.blowup.nc")

def column(i):
    lat0 = lat[i]
    x = R * rad * (lon - lon[i]) * np.cos(lat0 * rad)
    y = R * rad * (lat - lat[i])
    et = np.where(np.any(elem == i, axis=1))[0]
    nlev = int(nln[i]) - 1
    div = np.zeros(nlev)
    print(f"\n node idx1={i+1} lon {lon[i]:.3f} lat {lat[i]:.3f} ul={cavn[i]} nl={nln[i]} "
          f"({et.size} elements)")
    for e in et:
        nds = elem[e]
        cx, cy = x[nds].mean(), y[nds].mean()
        others = [k for k in nds if k != i]
        u = np.asarray(ds.variables["u"][0, e, :])
        v = np.asarray(ds.variables["v"][0, e, :])
        he = np.asarray(ds.variables["helem"][0, e, :])
        ul, nl = int(cave[e]), int(nle[e])
        # two sub-segments: edge midpoint (i,other) -> centroid
        tot = np.zeros(nlev)
        for o in others:
            mx, my = 0.5 * (x[i] + x[o]), 0.5 * (y[i] + y[o])
            dx, dy = cx - mx, cy - my
            # outward normal of the segment for node i's cell: rotate (dx,dy)
            # sign fixed below by requiring the far vertex to be on the inside
            nxv, nyv = dy, -dx
            # ensure the normal points away from node i
            if (nxv * (mx - x[i]) + nyv * (my - y[i])) < 0:
                nxv, nyv = -nxv, -nyv
            for k in range(ul - 1, min(nl - 1, nlev)):
                tot[k] += (u[k] * nxv + v[k] * nyv) * he[k]
        div += tot
        print(f"   elem {e+1:7d} ul={ul:2d} |flux| lvl1={abs(tot[0]):10.2f} "
              f"lvl5={abs(tot[4]):10.2f} lvl7={abs(tot[6]):10.2f} sum={tot.sum():12.2f}")
    # cumulative from the bottom, mirroring vert_vel_ale
    cum = np.zeros(nlev + 1)
    for k in range(nlev - 1, -1, -1):
        cum[k] = cum[k + 1] - div[k]
    wmod = np.asarray(ds.variables["w"][0, i, :])
    print("   lvl   div(flux)      cum(flux)     W_model      cum/W  (implied area)")
    for k in range(nlev):
        r = cum[k] / wmod[k] if wmod[k] != 0 else np.nan
        print(f"   {k+1:3d} {div[k]:13.3f} {cum[k]:13.3f} {wmod[k]:13.4e} {r:14.4g}")

for LO, LA in ((-100.382545, -72.8639145), (-100.112221, -72.7109299)):
    i = int(np.argmin((lon - LO) ** 2 + (lat - LA) ** 2))
    column(i)
