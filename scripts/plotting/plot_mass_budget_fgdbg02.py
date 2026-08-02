#!/usr/bin/env python3
"""fgdbg02 Antarctic mass budget over the 240 coupled PISM years.

Reads PISM's scalar time series (ts_*.nc, annual) from every run_pism_* leg:
  tendency_of_ice_mass_due_to_surface_mass_flux   SMB accumulation   (+)
  tendency_of_ice_mass_due_to_basal_mass_flux     basal/cavity melt  (-)
  tendency_of_ice_mass_due_to_discharge           calving+frontal    (-)
already in Gt yr^-1, so no unit conversion.

One axis (all four series share Gt/yr). The sum is the net dM/dt and is
cross-checked against the independent d(ice_mass)/dt.
"""
import glob, os, re
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02"
OUT = "/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/report"

# dataviz slots 1-4 (validated; the aqua/yellow contrast WARN is relieved by direct labels)
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d8d8d4"

V_SMB = "tendency_of_ice_mass_due_to_surface_mass_flux"
V_BAS = "tendency_of_ice_mass_due_to_basal_mass_flux"
V_DIS = "tendency_of_ice_mass_due_to_discharge"

files = sorted(glob.glob(f"{E}/run_pism_*/work/*_pismr_ts_*.nc"),
               key=lambda f: int(re.search(r"ts_(\d+)0101", f).group(1)))
yr, smb, bas, dis, mass = [], [], [], [], []
for i, f in enumerate(files):
    d = nc.Dataset(f)
    n = d.variables[V_SMB].shape[0]
    yr.extend([10 * i + k for k in range(n)])
    smb.extend(np.array(d.variables[V_SMB][:]))
    bas.extend(np.array(d.variables[V_BAS][:]))
    dis.extend(np.array(d.variables[V_DIS][:]))
    mass.extend(np.array(d.variables["ice_mass"][:]) / 1e12)   # kg -> Gt
    d.close()

yr = np.array(yr, float)
smb, bas, dis = np.array(smb), np.array(bas), np.array(dis)
net = smb + bas + dis
mass = np.array(mass)

# independent check: d(ice_mass)/dt from the mass series itself
dmdt = np.full_like(net, np.nan)
dmdt[1:] = np.diff(mass) / np.diff(yr)

fig, axs = plt.subplots(2, 1, figsize=(12.0, 8.2), facecolor="white",
                        gridspec_kw={"height_ratios": [2.1, 1]})
fig.suptitle("fgdbg02: Antarctic ice-sheet mass budget over 240 coupled PISM years",
             fontsize=13, color=INK, y=0.97)

def tidy(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=9)

ax = axs[0]
ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
series = [("SMB accumulation", smb, BLUE),
          ("basal mass balance", bas, ORANGE),
          ("calving + frontal melt", dis, AQUA),
          ("net (sum)", net, YELLOW)]
for lab, y, c in series:
    ax.plot(yr, y, lw=2.0 if lab != "net (sum)" else 2.6, color=c, zorder=3,
            label=lab)
# direct labels at the right edge (relieves the contrast WARN; <=4 series)
for lab, y, c in series:
    ax.annotate(f"{lab}: {y[-1]:.0f}", xy=(yr[-1], y[-1]),
                xytext=(8, 0), textcoords="offset points",
                va="center", fontsize=9, color=MUTED)
ax.set_ylabel("mass tendency [Gt yr$^{-1}$]", color=MUTED, fontsize=10)
ax.set_xlabel("PISM year", color=MUTED, fontsize=10)
ax.set_title("(a) Budget terms, PISM sign convention (positive = ice gain). Net stays negative "
             "throughout, but the loss shrinks fourfold as calving relaxes.",
             fontsize=10.5, color=INK, loc="left")
ax.set_xlim(yr[0], yr[-1] + 40)
ax.set_ylim(-16000, 5500)
ax.text(0.015, 0.97, f"year 0 spin-up shock: calving {dis[0]:.0f}, net {net[0]:.0f} Gt/yr (off scale)",
        transform=ax.transAxes, fontsize=8.5, color=MUTED, va="top")
ax.legend(fontsize=9, frameon=False, labelcolor=MUTED, loc="lower left", ncol=2)
tidy(ax)

ax = axs[1]
ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
ax.plot(yr, net, lw=2.2, color=YELLOW, zorder=3, label="net = SMB + basal + calving")
ax.plot(yr, dmdt, lw=1.4, color=MUTED, ls="--", zorder=4,
        label="d(ice_mass)/dt  (independent check)")
ax.set_ylabel("net tendency [Gt yr$^{-1}$]", color=MUTED, fontsize=10)
ax.set_xlabel("PISM year", color=MUTED, fontsize=10)
ax.set_title("(b) Net budget against the mass series it should reproduce",
             fontsize=10.5, color=INK, loc="left")
ax.legend(fontsize=9, frameon=False, labelcolor=MUTED, loc="lower right")
tidy(ax)

fig.tight_layout(rect=(0, 0, 1, 0.945))
fig.savefig(f"{OUT}/fgdbg02_mass_budget.png", dpi=150, facecolor="white")

def stats(name, a):
    return (f"{name:22s} first10={np.mean(a[:10]):9.1f}  last10={np.mean(a[-10:]):9.1f}  "
            f"min={a.min():9.1f} max={a.max():9.1f}")
print("PISM years:", int(yr[0]), "->", int(yr[-1]), f"({len(yr)} annual values, {len(files)} legs)")
print(stats("SMB accumulation", smb))
print(stats("basal (cavity) melt", bas))
print(stats("calving + discharge", dis))
print(stats("net (sum)", net))
print(f"ice_mass {mass[0]:.0f} -> {mass[-1]:.0f} Gt  (change {mass[-1]-mass[0]:+.0f} Gt)")
print(f"mean net {np.nanmean(net):.1f} Gt/yr vs mean d(mass)/dt {np.nanmean(dmdt):.1f} Gt/yr")
