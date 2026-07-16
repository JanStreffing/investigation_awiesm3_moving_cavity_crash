#!/usr/bin/env python3
"""movcav25 (implicit ice-T solve) vs movcav8 (old best) vs the 100s-yr-stable
v3.4 reference: SH sea-ice + snow seasonal cycles, and December snow-depth maps
(movcav25 vs reference, same scale). Same cartopy machinery as plot_blowup.py /
plot_snow_extent.py. Reference nodes with >50 m snow (known Ronne artifact
cluster, 26 nodes) are excluded."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc

OUT = "/home/a/a270092/esm_tools/movcav4_crash_plots/snow_ice_mc25_vs_ref.png"
RUNS = {  # label -> (basepath, year, line style)
    "movcav25 (implicit solve)": ("/work/ab0246/a270092/runtime/awiesm3-develop-is/movcav25/outdata/fesom/", "1900", dict(color="tab:blue", lw=2)),
    "movcav8 (old best)":        ("/work/ab0246/a270092/pism_repro/experiments/movcav8/outdata/fesom/",      "1900", dict(color="tab:orange", lw=2)),
    "v3.4 reference (stable 100s yr)": ("/work/bb1469/a270270/runtime/awiesm3-v3.4/AWIESM720_CMIP7_SPINUP_TCO95_CORE3/outdata/fesom/", "1849", dict(color="k", lw=2, ls="--")),
}
LAT_CUT = -55.0

try:
    import cartopy.crs as ccrs
    import matplotlib.path as mpath
    proj = ccrs.SouthPolarStereo(); tr = ccrs.PlateCarree(); HAVE_CARTO = True
except Exception as e:
    print("no cartopy:", e); HAVE_CARTO = False


def _fix_geoaxes(ax):
    if not hasattr(ax, "_autoscaleXon"):
        ax._autoscaleXon = ax.get_autoscalex_on()
        ax._autoscaleYon = ax.get_autoscaley_on()


def load(base, year):
    sn = nc.Dataset(base + f"m_snow.fesom.{year}.nc")
    lat = np.ma.filled(sn["lat"][:], np.nan)
    lon = np.ma.filled(sn["lon"][:], np.nan)
    snow = np.ma.filled(sn["m_snow"][:].astype("f8"), np.nan)
    ai = np.ma.filled(nc.Dataset(base + f"a_ice.fesom.{year}.nc")["a_ice"][:].astype("f8"), np.nan)
    snow = np.where(snow > 50.0, np.nan, snow)          # drop known artifact nodes
    return lat, lon, snow, ai


def cycles(lat, snow, ai):
    """monthly SH ice-node count + SH mean snow depth over icy nodes"""
    sh = lat < 0
    n_ice, dep = [], []
    for t in range(snow.shape[0]):
        m = sh & (ai[t] > 0.15)
        n_ice.append(int(np.sum(m & np.isfinite(snow[t]))))
        d = np.where(ai[t] > 0.05, snow[t] / np.where(ai[t] > 0.05, ai[t], np.nan), np.nan)[m]
        d = d[np.isfinite(d)]
        dep.append(float(np.mean(d)) if d.size else np.nan)
    return np.array(n_ice), np.array(dep)


def polar_panel(fig, pos, lon, lat, val, title, vmax=2.5):
    m = np.isfinite(val) & (lat < LAT_CUT)
    if HAVE_CARTO:
        ax = fig.add_subplot(*pos, projection=proj)
        _fix_geoaxes(ax)
        theta = np.linspace(0, 2 * np.pi, 200)
        circ = mpath.Path(np.stack([np.cos(theta), np.sin(theta)], 1) * 0.5 + 0.5)
        ax.set_boundary(circ, transform=ax.transAxes)
        try:
            ax.set_extent([-180, 180, -90, LAT_CUT], crs=tr)
        except Exception as e:
            print("set_extent fallback:", e)
        sc = ax.scatter(lon[m], lat[m], c=val[m], s=3.5, cmap="YlGnBu", vmin=0, vmax=vmax,
                        transform=tr, edgecolors="none")
        ax.coastlines(linewidth=0.6, color="k", zorder=5)
        ax.gridlines(linewidth=0.3, color="grey")
    else:
        ax = fig.add_subplot(*pos)
        th = np.deg2rad(lon[m]); r = 90 + lat[m]
        sc = ax.scatter(r * np.cos(th), r * np.sin(th), c=val[m], s=3.5, cmap="YlGnBu",
                        vmin=0, vmax=vmax, edgecolors="none")
        ax.set_aspect("equal")
    plt.colorbar(sc, ax=ax, shrink=0.72, pad=0.03, label="snow depth (m)")
    ax.set_title(title, fontsize=11)


mon = np.arange(1, 13)
fig = plt.figure(figsize=(15.5, 10))

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
maps = {}
for label, (base, year, style) in RUNS.items():
    lat, lon, snow, ai = load(base, year)
    n_ice, dep = cycles(lat, snow, ai)
    ax1.plot(mon, n_ice / 1000.0, "-o", ms=4, label=label, **style)
    ax2.plot(mon, dep, "-o", ms=4, label=label, **style)
    # December actual depth for the maps
    d = np.where(ai[-1] > 0.05, snow[-1] / np.where(ai[-1] > 0.05, ai[-1], np.nan), np.nan)
    maps[label] = (lon, lat, d)

ax1.set_title("SH sea-ice extent (nodes with a_ice > 0.15)", fontsize=11)
ax1.set_ylabel("ice nodes (thousands)")
ax1.annotate("reference: near-total\nsummer collapse", xy=(2, 0.4), xytext=(3.2, 12),
             fontsize=8, color="k", arrowprops=dict(arrowstyle="->", lw=0.8))
ax2.set_title("SH mean snow depth on sea ice", fontsize=11)
ax2.set_ylabel("snow depth (m)")
for ax in (ax1, ax2):
    ax.set_xlabel("month"); ax.set_xticks(mon); ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")

l25 = "movcav25 (implicit solve)"
lrf = "v3.4 reference (stable 100s yr)"
polar_panel(fig, (2, 2, 3), *maps[l25], f"movcav25 — December snow depth\n(max {np.nanmax(maps[l25][2][maps[l25][1] < LAT_CUT]):.2f} m)")
polar_panel(fig, (2, 2, 4), *maps[lrf], f"v3.4 reference — December snow depth\n(max {np.nanmax(maps[lrf][2][maps[lrf][1] < LAT_CUT]):.2f} m, artifacts excluded)")

fig.suptitle("Sea ice & snow, year 1: implicit-solve run vs old best vs the stable reference\n"
             "(cycle bounded and crash-free now, but SH summer melt-back still far weaker than reference)",
             fontsize=12.5, y=0.995)
fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("saved", OUT)
