#!/usr/bin/env python3
"""
Plot FESOM blowup ice variables (2D node fields).
"""

import sys
sys.path.insert(0, '/work/ab0246/a270092/software/pyfesom2')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
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

# Region
lon_min, lon_max = 25, 35
lat_min, lat_max = 38, 44

print("Loading mesh...")
mesh = pf.load_mesh(MESH_PATH)
lons = mesh.x2
lats = mesh.y2

# Find crash node
dist = np.sqrt((lons - crash_lon)**2 + (lats - crash_lat)**2)
crash_idx = np.argmin(dist)
print(f"Crash node: {crash_idx}")

# Filter region
region_mask = (lons >= lon_min) & (lons <= lon_max) & (lats >= lat_min) & (lats <= lat_max)
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

print("Loading ice data...")
with nc.Dataset(blowup_file, 'r') as ds:
    a_ice = ds.variables['a_ice'][0, :]
    m_ice = ds.variables['m_ice'][0, :]
    m_snow = ds.variables['m_snow'][0, :]
    u_ice = ds.variables['u_ice'][0, :]
    v_ice = ds.variables['v_ice'][0, :]
    iter_count = ds.variables['iter'][0]

print(f"Ice values at crash:")
print(f"  a_ice  = {a_ice[crash_idx]:.4e}")
print(f"  m_ice  = {m_ice[crash_idx]:.4e}")
print(f"  m_snow = {m_snow[crash_idx]:.4e}")
print(f"  u_ice  = {u_ice[crash_idx]:.4e}")
print(f"  v_ice  = {v_ice[crash_idx]:.4e}")

ice_vars = {
    'a_ice': (a_ice, 'Ice Concentration', '0-1', 'Blues'),
    'm_ice': (m_ice, 'Ice Mass', 'kg/m²', 'viridis'),
    'm_snow': (m_snow, 'Snow Mass', 'kg/m²', 'plasma'),
    'u_ice': (u_ice, 'Ice Velocity U', 'm/s', 'RdBu_r'),
    'v_ice': (v_ice, 'Ice Velocity V', 'm/s', 'RdBu_r'),
}

print("\nCreating ice overview...")
fig, axes = plt.subplots(3, 2, figsize=(18, 20))
axes = axes.flatten()

for idx, (varname, (data, title, units, cmap)) in enumerate(ice_vars.items()):
    if idx >= len(axes):
        break
        
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

if len(ice_vars) < len(axes):
    axes[-1].axis('off')

fig.suptitle(f'FESOM Blowup - Ice Variables\nIteration {iter_count}',
             fontsize=22, fontweight='bold')
plt.tight_layout()

outfile = f'{output_dir}/ice_overview.png'
fig.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"Saved: {outfile}")
plt.close(fig)

print("\n" + "="*60)
print("ICE PLOTS DONE!")
print("="*60)
