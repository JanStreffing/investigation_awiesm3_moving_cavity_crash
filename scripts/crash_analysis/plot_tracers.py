#!/usr/bin/env python3
"""
Plot FESOM blowup tracers (temp, salt) - using pyfesom2 mesh + regional filtering.
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
elem = mesh.elem

print(f"Mesh: {mesh.n2d} nodes, {elem.shape[0]} triangles")

# Find crash node
dist = np.sqrt((lons - crash_lon)**2 + (lats - crash_lat)**2)
crash_idx = np.argmin(dist)
print(f"Crash node: {crash_idx} at ({lons[crash_idx]:.4f}, {lats[crash_idx]:.4f})")

# Filter to region - only include nodes in the box
region_mask = (lons >= lon_min) & (lons <= lon_max) & (lats >= lat_min) & (lats <= lat_max)
region_indices = np.where(region_mask)[0]
print(f"Nodes in region: {len(region_indices)}")

# Create regional triangulation from filtered nodes
region_lons = lons[region_mask]
region_lats = lats[region_mask]
triang = tri.Triangulation(region_lons, region_lats)

# Mask triangles with long edges (avoid artifacts)
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

print("Loading tracer data...")
with nc.Dataset(blowup_file, 'r') as ds:
    temp = ds.variables['temp'][0, :, :]
    salt = ds.variables['salt'][0, :, :]
    hnode = ds.variables['hnode'][0, :, :]
    iter_count = ds.variables['iter'][0]

nlevels = temp.shape[1]
print(f"Vertical levels: {nlevels}")

# Compute depths at crash point
depths = np.zeros(nlevels)
for k in range(nlevels):
    if k == 0:
        depths[k] = hnode[crash_idx, k] / 2
    else:
        depths[k] = depths[k-1] + (hnode[crash_idx, k-1] + hnode[crash_idx, k]) / 2

# Print crash point values
print(f"\n{'='*60}")
print(f"TRACER VALUES AT CRASH POINT (node {crash_idx})")
print(f"{'='*60}")
print(f"{'Level':>5} {'Depth(m)':>10} {'Temp(C)':>12} {'Salt(psu)':>12}")
for k in range(min(10, nlevels)):
    print(f"{k:>5} {depths[k]:>10.1f} {temp[crash_idx, k]:>12.4f} {salt[crash_idx, k]:>12.4f}")
if nlevels > 10:
    print(f"  ... ({nlevels - 10} more levels)")

tracers = {
    'temp': (temp, 'Temperature', '°C', 'RdYlBu_r'),
    'salt': (salt, 'Salinity', 'PSU', 'viridis'),
}

for varname, (data, title, units, cmap) in tracers.items():
    print(f"\nPlotting {varname}...")
    
    # Get regional data
    region_data = data[region_mask, 0]  # surface layer
    
    # ===== Figure 1: Surface map =====
    fig1 = plt.figure(figsize=(12, 10))
    ax1 = fig1.add_subplot(111, projection=ccrs.PlateCarree())
    
    vmin, vmax = np.nanpercentile(region_data[np.isfinite(region_data)], [2, 98])
    
    cf = ax1.tripcolor(triang, region_data, cmap=cmap, vmin=vmin, vmax=vmax,
                       transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cf, ax=ax1, orientation='vertical', shrink=0.7, pad=0.02)
    cbar.set_label(f'{title} [{units}]', fontsize=12)
    
    ax1.plot(crash_lon, crash_lat, 'k*', markersize=20, transform=ccrs.PlateCarree())
    ax1.plot(crash_lon, crash_lat, 'r*', markersize=15, transform=ccrs.PlateCarree())
    
    ax1.coastlines(resolution='10m', linewidth=1, color='black')
    ax1.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.5)
    ax1.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    gl = ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    
    ax1.set_title(f'Surface {title}\nCrash location: ({crash_lon:.2f}°E, {crash_lat:.2f}°N)\n'
                  f'Value at crash: {data[crash_idx, 0]:.4f} {units}', fontsize=14, fontweight='bold')
    
    outfile = f'{output_dir}/tracers_{varname}_surface.png'
    fig1.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close(fig1)
    
    # ===== Figure 2: Vertical profile at crash point =====
    fig2, ax2 = plt.subplots(figsize=(8, 10))
    
    profile = data[crash_idx, :]
    valid = np.isfinite(profile) & (profile != 0)
    
    # Use absolute depths for y-axis
    ax2.plot(profile[valid], np.abs(depths[valid]), 'o-', linewidth=2, markersize=4, color='blue')
    ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax2.set_xlabel(f'{title} [{units}]', fontsize=12)
    ax2.set_ylabel('Depth [m]', fontsize=12)
    ax2.set_title(f'{title} Profile at Crash Point\n({crash_lon:.2f}°E, {crash_lat:.2f}°N)', 
                  fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.invert_yaxis()  # 0 at top, increasing downward
    
    outfile = f'{output_dir}/tracers_{varname}_profile.png'
    fig2.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close(fig2)
    
    # ===== Figure 3: Zonal section =====
    fig3, ax3 = plt.subplots(figsize=(14, 8))
    
    lat_band = np.abs(lats - crash_lat) < 0.1
    section_mask = lat_band & (lons >= lon_min) & (lons <= lon_max)
    section_idx = np.where(section_mask)[0]
    
    if len(section_idx) > 10:
        sort_order = np.argsort(lons[section_idx])
        section_idx = section_idx[sort_order]
        section_lons = lons[section_idx]
        section_data = data[section_idx, :]
        section_hnode = hnode[section_idx, :]
        
        section_depths = np.zeros_like(section_data)
        for i in range(len(section_idx)):
            for k in range(nlevels):
                if k == 0:
                    section_depths[i, k] = section_hnode[i, k] / 2
                else:
                    section_depths[i, k] = section_depths[i, k-1] + (section_hnode[i, k-1] + section_hnode[i, k]) / 2
        
        LON = np.tile(section_lons[:, np.newaxis], (1, nlevels))
        vmin_s, vmax_s = np.nanpercentile(section_data[section_data != 0], [2, 98])
        section_data_masked = np.ma.masked_where((section_data == 0) | ~np.isfinite(section_data), section_data)
        
        cf = ax3.pcolormesh(LON, -section_depths, section_data_masked, cmap=cmap, 
                            vmin=vmin_s, vmax=vmax_s, shading='auto')
        plt.colorbar(cf, ax=ax3, shrink=0.8, pad=0.02, label=f'{title} [{units}]')
        ax3.axvline(crash_lon, color='red', linestyle='--', linewidth=2, label='Crash')
        ax3.set_xlabel('Longitude [°E]')
        ax3.set_ylabel('Depth [m]')
        ax3.set_title(f'{title} - Zonal Section at {crash_lat:.2f}°N', fontsize=14, fontweight='bold')
        ax3.legend()
        
        outfile = f'{output_dir}/tracers_{varname}_section_zonal.png'
        fig3.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f"  Saved: {outfile}")
    plt.close(fig3)
    
    # ===== Figure 4: Meridional section =====
    fig4, ax4 = plt.subplots(figsize=(14, 8))
    
    lon_band = np.abs(lons - crash_lon) < 0.1
    section_mask = lon_band & (lats >= lat_min) & (lats <= lat_max)
    section_idx = np.where(section_mask)[0]
    
    if len(section_idx) > 10:
        sort_order = np.argsort(lats[section_idx])
        section_idx = section_idx[sort_order]
        section_lats = lats[section_idx]
        section_data = data[section_idx, :]
        section_hnode = hnode[section_idx, :]
        
        section_depths = np.zeros_like(section_data)
        for i in range(len(section_idx)):
            for k in range(nlevels):
                if k == 0:
                    section_depths[i, k] = section_hnode[i, k] / 2
                else:
                    section_depths[i, k] = section_depths[i, k-1] + (section_hnode[i, k-1] + section_hnode[i, k]) / 2
        
        LAT = np.tile(section_lats[:, np.newaxis], (1, nlevels))
        vmin_s, vmax_s = np.nanpercentile(section_data[section_data != 0], [2, 98])
        section_data_masked = np.ma.masked_where((section_data == 0) | ~np.isfinite(section_data), section_data)
        
        cf = ax4.pcolormesh(LAT, -section_depths, section_data_masked, cmap=cmap, 
                            vmin=vmin_s, vmax=vmax_s, shading='auto')
        plt.colorbar(cf, ax=ax4, shrink=0.8, pad=0.02, label=f'{title} [{units}]')
        ax4.axvline(crash_lat, color='red', linestyle='--', linewidth=2, label='Crash')
        ax4.set_xlabel('Latitude [°N]')
        ax4.set_ylabel('Depth [m]')
        ax4.set_title(f'{title} - Meridional Section at {crash_lon:.2f}°E', fontsize=14, fontweight='bold')
        ax4.legend()
        
        outfile = f'{output_dir}/tracers_{varname}_section_merid.png'
        fig4.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f"  Saved: {outfile}")
    plt.close(fig4)

# ===== Combined overview =====
print("\nCreating combined overview...")
fig, axes = plt.subplots(2, 2, figsize=(18, 12))

for row, (varname, (data, title, units, cmap)) in enumerate(tracers.items()):
    region_data = data[region_mask, 0]
    vmin, vmax = np.nanpercentile(region_data[np.isfinite(region_data)], [2, 98])
    
    # Surface map
    ax = axes[row, 0]
    tpc = ax.tripcolor(triang, region_data, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.plot(crash_lon, crash_lat, 'r*', markersize=20)
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

fig.suptitle(f'FESOM Blowup - Tracers Overview',
             fontsize=22, fontweight='bold')
plt.tight_layout()

outfile = f'{output_dir}/tracers_overview.png'
fig.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"Saved: {outfile}")
plt.close(fig)

print("\n" + "="*60)
print("TRACER PLOTS DONE!")
print("="*60)
