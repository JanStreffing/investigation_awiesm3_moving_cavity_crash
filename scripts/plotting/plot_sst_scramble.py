#!/usr/bin/env python
"""Global January SST as OIFS sees it: healthy (movcav4, pool weights) vs
scrambled (movcav15, ocp-tool weights). If the exchange is tile-translocated,
it should be instantly visible."""
import numpy as np, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT="/work/ab0246/a270092/postprocessing/investigation_awiesm3_moving_cavity_crash/figures/spp"
fig, axs = plt.subplots(2, 1, figsize=(13, 11))
for ax, run, tag in [(axs[0], "movcav4",  "movcav4 — pool weights (healthy)"),
                     (axs[1], "movcav15", "movcav15 — ocp-tool regenerated weights")]:
    d = nc.Dataset(f"/work/ab0246/a270092/pism_repro/experiments/{run}/outdata/oifs/atm_remapped_1m_sst_1m_1900-1900.nc")
    s = np.squeeze(np.asarray(d["sst"][:]))[0]          # January
    la = np.asarray(d["lat"][:]); lo = np.asarray(d["lon"][:])
    s = np.where(s > 1e10, np.nan, s)
    if np.nanmean(s) > 100: s = s - 273.15
    pc = ax.pcolormesh(lo, la, s, cmap="RdYlBu_r", vmin=-2, vmax=30, shading="nearest")
    plt.colorbar(pc, ax=ax, label="SST (°C)")
    ax.set_title(f"OIFS January SST — {tag}", fontsize=12)
    ax.set_xlabel("lon"); ax.set_ylabel("lat")
fig.suptitle("The ocean surface the ATMOSPHERE is coupled to (January, month 1)", fontsize=13)
fig.tight_layout()
fig.savefig(f"{OUT}/sst_scramble_global.png", dpi=110)
print(f"wrote {OUT}/sst_scramble_global.png")
