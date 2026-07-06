#!/usr/bin/env python3
"""
Plot FESOM blowup layer thickness variables.
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

print("Loading layer thickness data...")
with nc.Dataset(blowup_file, 'r') as ds:
    hnode = ds.variables['hnode'][0, :, :]  # (node, nz_1)
    zbar_n_bot = ds.variables['zbar_n_bot'][0, :]
    bottom_node_thickness = ds.variables['bottom_node_thickness'][0, :]
    iter_count = ds.variables['iter'][0]

nlevels = hnode.shape[1]
print(f"Vertical levels: {nlevels}")

# Compute depths
depths = np.zeros(nlevels)
for k in range(nlevels):
    if k == 0:
        depths[k] = hnode[crash_idx, k] / 2
    else:
        depths[k] = depths[k-1] + (hnode[crash_idx, k-1] + hnode[crash_idx, k]) / 2

print(f"\nLayer thickness at crash:")
print(f"  zbar_n_bot = {zbar_n_bot[crash_idx]:.4e} m")
print(f"  bottom_node_thickness = {bottom_node_thickness[crash_idx]:.4e} m")
print(f"  Top 5 layer thicknesses:")
for k in range(min(5, nlevels)):
    print(f"    Level {k}: hnode = {hnode[crash_idx, k]:.4e} m")

thickness_vars = {
    'hnode_surface': (hnode[:, 0], 'Surface Layer Thickness (hnode[0])', 'm', 'viridis'),
    'zbar_n_bot': (zbar_n_bot, 'Bottom Depth (zbar_n_bot)', 'm', 'RdBu_r'),
    'bottom_node_thickness': (bottom_node_thickness, 'Bottom Layer Thickness', 'm', 'plasma'),
}

print("\nCreating layer thickness overview...")
fig, axes = plt.subplots(2, 2, figsize=(18, 12))
axes = axes.flatten()

for idx, (varname, (data, title, units, cmap)) in enumerate(thickness_vars.items()):
    ax = axes[idx]
    region_data = data[region_mask]
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
    ax.set_title(f'{title}\nValue at crash: {data[crash_idx]:.2e} {units}', 
                 fontsize=16, fontweight='bold')
    ax.set_xlabel('Longitude [°E]', fontsize=14)
    ax.set_ylabel('Latitude [°N]', fontsize=14)
    ax.tick_params(labelsize=12)
    cbar = plt.colorbar(tpc, ax=ax, shrink=0.9)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(f'{units}', fontsize=14)

# Add hnode profile in 4th panel
ax = axes[3]
profile = hnode[crash_idx, :]
valid = np.isfinite(profile) & (profile != 0)
ax.plot(profile[valid], np.abs(depths[valid]), 'o-', markersize=6, linewidth=2.5, color='blue')
ax.set_xlabel('Layer Thickness [m]', fontsize=16)
ax.set_ylabel('Depth [m]', fontsize=16)
ax.set_title('Layer Thickness Profile at Crash', fontsize=18, fontweight='bold')
ax.tick_params(labelsize=14)
ax.grid(True, alpha=0.3, linewidth=1)
ax.invert_yaxis()

fig.suptitle(f'FESOM Blowup - Layer Thickness\nIteration {iter_count}',
             fontsize=22, fontweight='bold')
plt.tight_layout()

outfile = f'{output_dir}/layer_thickness_overview.png'
fig.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"Saved: {outfile}")
plt.close(fig)

print("\n" + "="*60)
print("LAYER THICKNESS PLOTS DONE!")
print("="*60)
