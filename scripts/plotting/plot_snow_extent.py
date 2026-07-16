#!/usr/bin/env python3
"""Extent of the Antarctic snow-on-sea-ice pile-up in the OLD (bounded) movcav
runs. South-Polar-Stereo, same cartopy machinery as plot_blowup.py, but coloured
by ACTUAL snow depth (= effective m_snow / a_ice) on the FESOM nodes. We scatter
on the lat/lon carried inside the output file itself (self-consistent) rather than
triangulating, because each moving-cavity leg regenerated the mesh and
latest_submesh no longer matches the year-1900 output node count.

Two views per run:
  (a) node map of the annual-MAX actual snow depth  -> where the problem lives
  (b) node map of December actual snow depth        -> the winter loading
Plus a seasonal-cycle line (global max eff snow, movcav8 vs movcav16) showing the
autumn melt-back that keeps the old runs bounded (~2-2.6 m)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc

EXPDIR = "/work/ab0246/a270092/pism_repro/experiments"
OUT    = "/home/a/a270092/esm_tools/movcav4_crash_plots/snow_extent_movcav.png"
RUNS   = ["movcav8", "movcav16"]          # old, bounded, completed a year
LAT_CUT = -55.0                            # Antarctic cap

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


def load(exp):
    b = f"{EXPDIR}/{exp}/outdata/fesom/"
    sn = nc.Dataset(b + "m_snow.fesom.1900.nc")
    lat = np.ma.filled(sn["lat"][:], np.nan)
    lon = np.ma.filled(sn["lon"][:], np.nan)
    snow = np.ma.filled(sn["m_snow"][:].astype("f8"), np.nan)     # (t, nod) effective
    ai = np.ma.filled(nc.Dataset(b + "a_ice.fesom.1900.nc")["a_ice"][:].astype("f8"), np.nan)
    depth = np.where(ai > 0.05, snow / np.where(ai > 0.05, ai, np.nan), np.nan)  # actual depth
    return lat, lon, snow, depth


def polar_panel(fig, pos, lon, lat, val, title, vmax, cmap="YlGnBu"):
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
        sc = ax.scatter(lon[m], lat[m], c=val[m], s=3.5, cmap=cmap, vmin=0, vmax=vmax,
                        transform=tr, edgecolors="none")
        ax.coastlines(linewidth=0.6, color="k", zorder=5)
        ax.gridlines(linewidth=0.3, color="grey")
    else:
        ax = fig.add_subplot(*pos)
        th = np.deg2rad(lon[m]); r = 90 + lat[m]
        sc = ax.scatter(r * np.cos(th), r * np.sin(th), c=val[m], s=3.5, cmap=cmap,
                        vmin=0, vmax=vmax, edgecolors="none")
        ax.set_aspect("equal")
    plt.colorbar(sc, ax=ax, shrink=0.72, pad=0.03, label="snow depth (m)")
    ax.set_title(title, fontsize=11)
    return ax


fig = plt.figure(figsize=(15, 10))
seas = {}
for j, exp in enumerate(RUNS):
    lat, lon, snow, depth = load(exp)
    amax = np.nanmax(depth, axis=0)            # annual-max actual depth per node
    dec = depth[-1]                            # December
    seas[exp] = np.array([np.nanmax(np.where(np.isfinite(snow[t]), snow[t], np.nan))
                          for t in range(snow.shape[0])])
    n2 = np.nansum(amax > 2.0); n15 = np.nansum(amax > 1.5)
    polar_panel(fig, (2, 3, 3 * j + 1), lon, lat, amax,
                f"{exp}: annual-MAX snow depth\n(nodes >2 m: {n2}, >1.5 m: {n15})",
                vmax=2.5)
    polar_panel(fig, (2, 3, 3 * j + 2), lon, lat, dec,
                f"{exp}: December snow depth\n(max {np.nanmax(dec):.2f} m)", vmax=2.5)

# seasonal-cycle line panel (spans right column, both runs)
axl = fig.add_subplot(1, 3, 3)
mon = np.arange(1, 13)
for exp, c in zip(RUNS, ["tab:blue", "tab:orange"]):
    axl.plot(mon, seas[exp], "-o", color=c, label=exp)
axl.axvspan(4, 6, color="grey", alpha=0.15)
axl.text(5, axl.get_ylim()[0], "autumn\nmelt-back", ha="center", va="bottom", fontsize=8)
axl.set_xlabel("month (1900)"); axl.set_ylabel("global-max effective snow (m)")
axl.set_title("Seasonal cycle — old runs stay bounded\n(autumn thin-back, then regrow)", fontsize=11)
axl.set_xticks(mon); axl.grid(alpha=0.3); axl.legend(fontsize=9)

fig.suptitle("Antarctic snow-on-sea-ice pile-up at the ice-shelf coast — old movcav runs "
             "(bounded ~2-2.6 m; real ~0.5-1 m)", fontsize=13, y=1.0)
fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("saved", OUT)
