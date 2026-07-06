#!/usr/bin/env python3
"""
Plot FESOM blowup mixing coefficients (3D fields).
"""

import sys
sys.path.insert(0, '/work/ab0246/a270092/software/pyfesom2')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import netCDF4 as nc
import pyfesom2 as pf

MESH_PATH = '/work/ab0246/a270092/input/fesom2/dars2'
work_dir = '/work/bb1469/a270092/runtime/awiesm3-develop/test-spinup/run_13580101-13591231_crash/work'
blowup_file = f'{work_dir}/fesom.1358.oce.blowup.nc'
output_dir = '/work/bb1469/a270092/runtime/awiesm3-develop/test-spinup/run_13580101-13591231_crash'

crash_lon = 29.7770298093963
crash_lat = 41.2424041260004
lon_min, lon_max = 25, 35
lat_min, lat_max = 38, 44

print("Loading mesh...")
mesh = pf.load_mesh(MESH_PATH)
lons = mesh.x2
lats = mesh.y2

dist = np.sqrt((lons - crash_lon)**2 + (lats - crash_lat)**2)
crash_idx = np.argmin(dist)
print(f"Crash node: {crash_idx}")

region_mask = (lons >= lon_min) & (lons <= lon_max) & (lats >= lat_min) & (lats <= lat_max)
region_lons = lons[region_mask]
region_lats = lats[region_mask]
triang = tri.Triangulation(region_lons, region_lats)

edge_threshold = 0.5
x, y = triang.x, triang.y
triangles = triang.triangles
mask = np.zeros(len(triangles), dtype=bool)
for i, t in enumerate(triangles):
    d01 = np.sqrt((x[t[0]] - x[t[1]])**2 + (y[t[0]] - y[t[1]])**2)
    d12 = np.sqrt((x[t[1]] - x[t[2]])**2 + (y[t[1]] - y[t[2]])**2)
    d20 = np.sqrt((x[t[2]] - x[t[0]])**2 + (y[t[2]] - y[t[0]])**2)
    if max(d01, d12, d20) > edge_threshold:
        mask[i] = True
triang.set_mask(mask)

print("Loading mixing data...")
with nc.Dataset(blowup_file, 'r') as ds:
    Kv = ds.variables['Kv'][0, :, :]  # (node, nz)
    N2 = ds.variables['N2'][0, :, :]
    hnode = ds.variables['hnode'][0, :, :]
    iter_count = ds.variables['iter'][0]

nlevels = Kv.shape[1]
print(f"Vertical levels: {nlevels}")

# Compute depths
depths = np.zeros(nlevels)
for k in range(nlevels):
    if k == 0:
        depths[k] = hnode[crash_idx, min(k, hnode.shape[1]-1)] / 2
    else:
        k_prev = min(k-1, hnode.shape[1]-1)
        k_curr = min(k, hnode.shape[1]-1)
        depths[k] = depths[k-1] + (hnode[crash_idx, k_prev] + hnode[crash_idx, k_curr]) / 2

print(f"\nMixing at crash (top 5 levels):")
for k in range(min(5, nlevels)):
    print(f"  Level {k}: Kv={Kv[crash_idx, k]:.4e}, N2={N2[crash_idx, k]:.4e}")

mixing_vars = {
    'Kv': (Kv, 'Vertical Diffusivity (Kv)', 'm²/s', 'viridis'),
    'N2': (N2, 'Buoyancy Frequency (N²)', '1/s²', 'RdYlBu_r'),
}

print("\nCreating mixing overview...")
fig, axes = plt.subplots(2, 2, figsize=(18, 12))

for row, (varname, (data, title, units, cmap)) in enumerate(mixing_vars.items()):
    # Surface map
    ax = axes[row, 0]
    region_data = data[region_mask, 0]
    finite_data = region_data[np.isfinite(region_data)]
    if len(finite_data) > 0:
        vmin, vmax = np.percentile(finite_data, [2, 98])
    else:
        vmin, vmax = -1, 1
    
    tpc = ax.tripcolor(triang, region_data, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.plot(crash_lon, crash_lat, 'k*', markersize=20)
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect('equal')
    ax.set_title(f'{title} Surface', fontsize=18, fontweight='bold')
    ax.set_xlabel('Longitude [°E]', fontsize=16)
    ax.set_ylabel('Latitude [°N]', fontsize=16)
    ax.tick_params(labelsize=14)
    cbar = plt.colorbar(tpc, ax=ax, shrink=0.8)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(f'{units}', fontsize=16)
    
    # Profile
    ax = axes[row, 1]
    profile = data[crash_idx, :]
    valid = np.isfinite(profile) & (profile != 0)
    ax.plot(profile[valid], np.abs(depths[valid]), 'o-', markersize=6, linewidth=2.5)
    ax.set_xlabel(f'{title} [{units}]', fontsize=16)
    ax.set_ylabel('Depth [m]', fontsize=16)
    ax.set_title(f'{title} Profile at Crash', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.3, linewidth=1)
    ax.invert_yaxis()

fig.suptitle(f'FESOM Blowup - Mixing Coefficients\nIteration {iter_count}',
             fontsize=22, fontweight='bold')
plt.tight_layout()

outfile = f'{output_dir}/mixing_overview.png'
fig.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"Saved: {outfile}")
plt.close(fig)

print("\n" + "="*60)
print("MIXING PLOTS DONE!")
print("="*60)
