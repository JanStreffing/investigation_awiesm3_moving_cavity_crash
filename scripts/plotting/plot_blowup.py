#!/usr/bin/env python3
"""PolyCollection map (ocp-tool style) of the FESOM blowup over the Antarctic.
Triangles from elem2d.out/nod2d.out, filled by field value; South-Polar-Stereo."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import netCDF4 as nc

MESH = "/work/ab0246/a270092/pism_repro/experiments/movcav4/couple/latest_submesh"
BLOW = "/work/ab0246/a270092/pism_repro/experiments/movcav4/run_awiesm3_19010101-19011231/work/fesom.1901.oce.blowup.nc"
OUT  = "/home/a/a270092/esm_tools/fesom_blowup_map.png"

# --- mesh -------------------------------------------------------------------
with open(f"{MESH}/nod2d.out") as f:
    nn = int(f.readline().split()[0])
    nod = np.loadtxt(f, max_rows=nn)
lon = nod[:, 1]; lat = nod[:, 2]                     # degrees
with open(f"{MESH}/elem2d.out") as f:
    ne = int(f.readline().split()[0])
    elem = np.loadtxt(f, max_rows=ne, dtype=int) - 1  # 0-based
cav_nlvls = np.loadtxt(f"{MESH}/cavity_nlvls.out", skiprows=1).astype(int) \
            if True else None
print(f"nodes={nn} elems={ne}")

# --- fields (node) ----------------------------------------------------------
ds = nc.Dataset(BLOW)
def node_field(name, reduce="max"):
    v = ds.variables[name][:]
    v = np.squeeze(np.asarray(v))
    if v.ndim == 2:                 # (node, nz) -> reduce over depth
        v = np.nanmax(np.abs(v), axis=1) if reduce == "max" else v[:, 0]
    return v
cflz = node_field("cfl_z", "max")
hflx = np.squeeze(ds.variables["heat_flux"][:])
temp = np.squeeze(ds.variables["temp"][:])[:, 0]       # surface temp
iscav = (cav_nlvls > 0)                                 # under ice shelf
print(f"cfl_z max={np.nanmax(cflz):.1f}  cav nodes={iscav.sum()}")

# --- triangles: verts + per-elem value, drop dateline wrap ------------------
LAT_CUT = -60.0                                         # crop to Antarctic ring
tlon = lon[elem]; tlat = lat[elem]
verts = np.stack([tlon, tlat], axis=-1)                 # (ne,3,2)
elat = tlat.mean(1)
wrap = (tlon.max(1) - tlon.min(1)) < 180
keep = wrap & (tlat.max(1) < LAT_CUT)                   # only cells fully south
def evals(nodef, red="mean"):
    a = nodef[elem]
    return np.nanmax(a, 1) if red == "max" else np.nanmean(a, 1)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.path as mpath
    proj = ccrs.SouthPolarStereo(); HAVE_CARTO = True
    tr = ccrs.PlateCarree()
except Exception as e:
    print("no cartopy:", e); HAVE_CARTO = False

def _fix_geoaxes(ax):
    # cartopy 0.20 + newer mpl: hold_limits() reads _autoscaleXon/Yon
    if not hasattr(ax, "_autoscaleXon"):
        ax._autoscaleXon = ax.get_autoscalex_on()
        ax._autoscaleYon = ax.get_autoscaley_on()

panels = [("cfl_z (max over depth)", evals(cflz, "max"), "inferno", 0, 30),
          ("heat_flux  W/m2",        evals(hflx, "mean"), "RdBu_r", -600, 600),
          ("surface temp  degC",     evals(temp, "mean"), "viridis", -3, 5)]

fig = plt.figure(figsize=(19, 6.4))
for i, (title, val, cmap, vmn, vmx) in enumerate(panels):
    if HAVE_CARTO:
        ax = fig.add_subplot(1, 3, i + 1, projection=proj)
        _fix_geoaxes(ax)
        # circular boundary + crop to the Antarctic cap
        theta = np.linspace(0, 2 * np.pi, 200)
        circ = mpath.Path(np.stack([np.cos(theta), np.sin(theta)], 1) * 0.5 + 0.5)
        ax.set_boundary(circ, transform=ax.transAxes)
        try:
            ax.set_extent([-180, 180, -90, LAT_CUT], crs=tr)
        except Exception as e:
            print("set_extent fallback:", e)
        pc = PolyCollection(verts[keep], array=val[keep], cmap=cmap,
                            edgecolors="none", transform=tr)
        pc.set_clim(vmn, vmx); ax.add_collection(pc)
        ax.coastlines(linewidth=0.6, color="k", zorder=5)
        ax.gridlines(linewidth=0.3, color="grey")
    else:
        # manual polar-stereo: r=(90+lat), theta=lon
        ax = fig.add_subplot(1, 3, i + 1)
        th = np.deg2rad(tlon); r = (90 + tlat)
        vx = r * np.cos(th); vy = r * np.sin(th)
        vv = np.stack([vx, vy], axis=-1)
        pc = PolyCollection(vv[keep], array=val[keep], cmap=cmap,
                            edgecolors="none")
        pc.set_clim(vmn, vmx); ax.add_collection(pc)
        ax.set_xlim(-40, 40); ax.set_ylim(-40, 40); ax.set_aspect("equal")
    plt.colorbar(pc, ax=ax, shrink=0.7, pad=0.02)
    ax.set_title(title, fontsize=12)

# overlay: cavity ring + high-cfl nodes on panel 1
fig.suptitle("FESOM blowup (mstep 712) — Antarctic; cfl_z / OIFS heat flux / SST",
             fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("saved", OUT)
