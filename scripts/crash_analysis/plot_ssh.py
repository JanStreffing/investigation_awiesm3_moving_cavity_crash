#!/usr/bin/env python3
"""
Plot FESOM blowup SSH-related variables (2D node fields).
"""

import sys
sys.path.insert(0, '/work/ab0246/a270092/software/pyfesom2')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import pyfesom2 as pf

# Paths
MESH_PATH = '/work/ab0246/a270092/input/fesom2/dars2'
work_dir = '/work/bb1469/a270092/runtime/awiesm3-develop/test-spinup/run_13580101-13591231_crash/work'
blowup_file = f'{work_dir}/fesom.1358.oce.blowup.nc'
output_dir = '/work/bb1469/a270092/runtime/awiesm3-develop/test-spinup/run_13580101-13591231_crash'

# Crash location
crash_lon = 29.7770298093963
crash_lat = 41.2424041260004

# Region to plot
lon_min, lon_max = 25, 35
lat_min, lat_max = 38, 44

print("Loading mesh with pyfesom2...")
mesh = pf.load_mesh(MESH_PATH)
lons = mesh.x2
lats = mesh.y2

print(f"Mesh: {mesh.n2d} nodes")

# Find crash node
dist = np.sqrt((lons - crash_lon)**2 + (lats - crash_lat)**2)
crash_idx = np.argmin(dist)
print(f"Crash node: {crash_idx} at ({lons[crash_idx]:.4f}, {lats[crash_idx]:.4f})")

# Filter to region
region_mask = (lons >= lon_min) & (lons <= lon_max) & (lats >= lat_min) & (lats <= lat_max)
region_indices = np.where(region_mask)[0]
print(f"Nodes in region: {len(region_indices)}")

# Create regional triangulation
region_lons = lons[region_mask]
region_lats = lats[region_mask]
triang = tri.Triangulation(region_lons, region_lats)

# Mask long edges
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

print("Loading SSH data...")
with nc.Dataset(blowup_file, 'r') as ds:
    eta_n = ds.variables['eta_n'][0, :]
    hbar = ds.variables['hbar'][0, :]
    d_eta = ds.variables['d_eta'][0, :]
    ssh_rhs = ds.variables['ssh_rhs'][0, :]
    ssh_rhs_old = ds.variables['ssh_rhs_old'][0, :]
    iter_count = ds.variables['iter'][0]

# Print crash point values
print(f"\n{'='*60}")
print(f"SSH VALUES AT CRASH POINT (node {crash_idx})")
print(f"{'='*60}")
print(f"  eta_n        = {eta_n[crash_idx]:.6e}")
print(f"  hbar         = {hbar[crash_idx]:.6e}")
print(f"  d_eta        = {d_eta[crash_idx]:.6e}")
print(f"  ssh_rhs      = {ssh_rhs[crash_idx]:.6e}")
print(f"  ssh_rhs_old  = {ssh_rhs_old[crash_idx]:.6e}")

ssh_vars = {
    'eta_n': (eta_n, 'Sea Surface Elevation (eta_n)', 'm', 'RdBu_r'),
    'hbar': (hbar, 'ALE Surface Elevation (hbar)', 'm', 'RdBu_r'),
    'd_eta': (d_eta, 'SSH Change from Solver (d_eta)', 'm', 'RdBu_r'),
    'ssh_rhs': (ssh_rhs, 'SSH RHS', '?', 'RdBu_r'),
    'ssh_rhs_old': (ssh_rhs_old, 'SSH RHS Old', '?', 'RdBu_r'),
}

# Create combined overview
print("\nCreating SSH overview...")
fig, axes = plt.subplots(3, 2, figsize=(18, 20))
axes = axes.flatten()

for idx, (varname, (data, title, units, cmap)) in enumerate(ssh_vars.items()):
    if idx >= len(axes):
        break
        
    ax = axes[idx]
    region_data = data[region_mask]
    
    # Robust percentile calc
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

# Hide unused subplot
if len(ssh_vars) < len(axes):
    axes[-1].axis('off')

fig.suptitle(f'FESOM Blowup - SSH Variables\nIteration {iter_count}',
             fontsize=22, fontweight='bold')
plt.tight_layout()

outfile = f'{output_dir}/ssh_overview.png'
fig.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"Saved: {outfile}")
plt.close(fig)

print("\n" + "="*60)
print("SSH PLOTS DONE!")
print("="*60)
