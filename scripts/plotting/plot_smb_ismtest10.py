"""Verification figure for the first OIFS-driven surface mass balance.

Panel (a) the delivered SMB on the PISM grid, (b) precipitation against surface
elevation with the CloudSat and ERA5 reference values, (c) OIFS's own output
against what the ISM-mapper delivered, on the same footprint.
"""

import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pyproj
from matplotlib.colors import TwoSlopeNorm

RUN = "/work/ab0246/a270092/runtime/awiesm3-develop-is/ismtest10"
OUT = "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/report"
SEC = 365 * 86400
RHO_I, RHO_W = 910.0, 1028.0
STEPS = 8760.0 * 1000.0  # OIFS accumulation per step, m -> mm/yr

# CloudSat and ERA5, Roussel et al. 2020, The Cryosphere 14, 2715
CLOUDSAT_CONT, ERA5_CONT = 186.0, 216.0
CLOUDSAT_PLAT, ERA5_PLAT = 29.0, 65.0

p = nc.Dataset(glob.glob(f"{RUN}/run_pism_*/work/*_pismr_input_*.nc")[0])
thk = np.ma.filled(np.squeeze(p["thk"][:]), 0.0)
topg = np.ma.filled(np.squeeze(p["topg"][:]), 0.0)
x, y = p["x"][:], p["y"][:]
flt = topg < -(RHO_I / RHO_W) * thk
usurf = np.where(flt, thk * (1 - RHO_I / RHO_W), topg + thk)
ice = thk > 10.0
grounded = ice & ~flt

s = nc.Dataset(f"{RUN}/outdata/ismm/antar_smb_1900.nc")
ann = lambda v: np.nanmean(np.ma.filled(s[v][:], np.nan), axis=0)
P = (ann("antar_Precip_liquid") + ann("antar_Precip_solid")) * SEC * 1000
E = ann("antar_Evaporation") * SEC
R = ann("antar_Runoff") * SEC * 1000
SMB = P - E - R


def oifs(v):
    d = nc.Dataset(glob.glob(f"{RUN}/outdata/oifs/atm_remapped_1m_{v}_1m_*.nc")[0])
    n = [k for k in d.variables if k not in
         ("lon", "lat", "time", "time_counter", "time_bnds", "time_centered",
          "time_centered_bounds", "time_counter_bounds", "lon_bnds", "lat_bnds",
          "axis_nbounds")][0]
    return d, np.nanmean(np.ma.filled(d[n][:], np.nan), axis=0)


d, lsp = oifs("lsp")
_, cp = oifs("cp")
_, ev = oifs("e")
_, ro = oifs("ro")
LON, LAT = np.meshgrid(d["lon"][:], d["lat"][:])
tr = pyproj.Transformer.from_crs(
    "EPSG:4326",
    "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 "
    "+datum=WGS84 +units=m",
    always_xy=True,
)
X, Y = tr.transform(LON, LAT)
ix = np.round((X - x[0]) / 8000.0).astype(int)
# The data arrays run opposite to the y coordinate variable. Verified against
# the file's own lat/lon: the northernmost ice cell is the Peninsula tip.
iy = np.round((-Y - y[0]) / 8000.0).astype(int)
inside = (LAT < -55) & (ix >= 0) & (ix < len(x)) & (iy >= 0) & (iy < len(y))
wgt = np.cos(np.deg2rad(LAT))


def sample(cond):
    m = np.zeros_like(inside)
    m[inside] = cond[iy[inside], ix[inside]]
    return m


def oifs_mean(a, m):
    return np.nansum(a[m] * wgt[m]) / np.nansum(wgt[m])


fig = plt.figure(figsize=(14.0, 5.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 0.95], wspace=0.42)

# --- (a) SMB map ---------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
fld = np.where(ice, SMB, np.nan)
im = ax.imshow(fld, origin="upper", cmap="YlGnBu",
               norm=matplotlib.colors.PowerNorm(0.55, vmin=0, vmax=600),
               extent=[x[0] / 1e3, x[-1] / 1e3, y[0] / 1e3, y[-1] / 1e3])
neg = np.where(ice & (SMB < 0), 1.0, np.nan)
ax.contourf(neg, levels=[0.5, 1.5], colors=["#d62728"],
            extent=[x[0] / 1e3, x[-1] / 1e3, y[0] / 1e3, y[-1] / 1e3],
            origin="upper")
ax.contour(np.where(ice, usurf, np.nan), levels=[2500], colors="k",
           linewidths=0.8,
           extent=[x[0] / 1e3, x[-1] / 1e3, y[0] / 1e3, y[-1] / 1e3],
           origin="upper")
ax.set_xlim(-2800, 2900)
ax.set_ylim(-2300, 2300)
ax.set_xlabel("x [km]")
ax.set_ylabel("y [km]")
ax.set_title("(a) delivered SMB $= P-E-R$", fontsize=11)
fig.colorbar(im, ax=ax, shrink=0.82, label="mm w.e. yr$^{-1}$")

# --- (b) precipitation against elevation ---------------------------------
ax = fig.add_subplot(gs[0, 1])
edges = np.arange(0, 4250, 250.0)
mid, val, cnt = [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = grounded & (usurf >= lo) & (usurf < hi)
    if m.sum() < 50:
        continue
    mid.append(0.5 * (lo + hi))
    val.append(P[m].mean())
    cnt.append(m.sum())
ax.plot(val, mid, "o-", color="#1f4e79", lw=2, ms=4, label="OIFS via ISM-mapper")
ax.axvline(CLOUDSAT_PLAT, color="#c00000", ls="--", lw=1.4)
ax.axvline(ERA5_PLAT, color="#c00000", ls=":", lw=1.4)
ax.axhline(2500, color="k", lw=0.8)
ax.text(CLOUDSAT_PLAT + 6, 4120, "CloudSat 29", color="#c00000", fontsize=8,
        va="top", rotation=90)
ax.text(ERA5_PLAT + 6, 4120, "ERA5 65", color="#c00000", fontsize=8, va="top",
        rotation=90)
ax.text(430, 2570, "plateau / margin split", fontsize=8, ha="right")
ax.set_xlabel("precipitation [mm w.e. yr$^{-1}$]")
ax.set_ylabel("surface elevation [m]", labelpad=2)
ax.set_title("(b) grounded ice, 250 m bins", fontsize=11)
ax.set_xlim(0, 450)
ax.set_ylim(0, 4200)
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc="lower right")

# --- (c) fidelity of the remap -------------------------------------------
ax = fig.add_subplot(gs[0, 2])
groups = [
    ("all ice", sample(ice), ice),
    ("margin\n$<2500$ m", sample(grounded & (usurf < 2500)), grounded & (usurf < 2500)),
    ("plateau\n$\\geq 2500$ m", sample(grounded & (usurf >= 2500)), grounded & (usurf >= 2500)),
]
pos = np.arange(len(groups))
o = [oifs_mean(lsp + cp, m) * STEPS for _, m, _ in groups]
i = [P[m].mean() for _, _, m in groups]
ax.bar(pos - 0.19, o, 0.36, label="OIFS, native output", color="#7f7f7f")
ax.bar(pos + 0.19, i, 0.36, label="delivered on the PISM grid", color="#1f4e79")
for k, (a, b) in enumerate(zip(o, i)):
    ax.text(k, max(a, b) + 8, f"{100*(b-a)/a:+.1f}%", ha="center", fontsize=8)
ax.axhline(CLOUDSAT_CONT, color="#c00000", ls="--", lw=1.4)
ax.text(2.48, CLOUDSAT_CONT - 4, "CloudSat 186", color="#c00000", fontsize=8,
        ha="right", va="top")
ax.axhline(ERA5_CONT, color="#c00000", ls=":", lw=1.4)
ax.text(2.48, ERA5_CONT + 3, "ERA5 216", color="#c00000", fontsize=8,
        ha="right", va="bottom")
ax.set_xticks(pos)
ax.set_xticklabels([g[0] for g in groups], fontsize=9)
ax.set_ylabel("precipitation [mm w.e. yr$^{-1}$]", labelpad=2)
ax.set_title("(c) same footprint, both sides", fontsize=11)
ax.set_ylim(0, 300)
ax.legend(fontsize=8, loc="upper left")
ax.grid(alpha=0.3, axis="y")

fig.savefig(f"{OUT}/ismtest10_smb_verification.png", dpi=150, bbox_inches="tight")
print("wrote", f"{OUT}/ismtest10_smb_verification.png")
