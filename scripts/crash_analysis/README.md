# FESOM2 Crash Analysis Tool

Automated plotting toolkit for analyzing FESOM2 model crashes from blowup files.

## Overview

This tool generates comprehensive visualizations of FESOM2 model state at crash time, including:
- **Tracers**: Temperature and salinity (surface maps + vertical profiles)
- **SSH Variables**: Sea surface elevation and related diagnostics
- **Ice Variables**: Ice concentration, mass, snow, and velocities
- **Surface Fluxes**: Heat and water fluxes
- **Vertical Velocity**: w, w_expl, w_impl, and CFL numbers
- **Mixing Coefficients**: Vertical diffusivity and buoyancy frequency
- **Layer Thickness**: Vertical grid structure diagnostics

## Quick Start

```bash
# 1. Edit configuration
vim config.yaml

# 2. Run all categories
python run_analysis.py --all

# 3. Run specific categories
python run_analysis.py --categories tracers,ssh,ice

# 4. Interactive mode
python run_analysis.py
```

## Configuration

Edit `config.yaml` to set:

```yaml
# Essential paths
mesh:
  path: /path/to/mesh

data:
  work_dir: /path/to/run/work
  blowup_file: fesom.XXXX.oce.blowup.nc
  output_dir: /path/to/output

# Crash location (from error log)
crash:
  lon: 29.777
  lat: 41.242

# Region bounds
region:
  lon_min: 25
  lon_max: 35
  lat_min: 38
  lat_max: 44
```

## Usage

### List Available Categories
```bash
python run_analysis.py --list
```

### Run All Enabled Categories
```bash
python run_analysis.py --all
```

### Run Specific Categories
```bash
python run_analysis.py --categories tracers
python run_analysis.py --categories ssh,ice,fluxes
```

### Interactive Mode
```bash
python run_analysis.py
```

### Use Custom Config
```bash
python run_analysis.py --config my_config.yaml --all
```

## Output

Plots are saved to the directory specified in `config.yaml` as:
- `tracers_overview.png`
- `ssh_overview.png`
- `ice_overview.png`
- `fluxes_overview.png`
- `vertical_velocity_overview.png`
- `mixing_overview.png`
- `layer_thickness_overview.png`

Plus individual plots for each variable (surface, profile, sections).

## File Structure

```
fesom2_crash_analysis/
├── config.yaml                    # Configuration file
├── run_analysis.py                # Main runner script
├── README.md                      # This file
├── plot_tracers.py               # Tracer plotting
├── plot_ssh.py                   # SSH variable plotting
├── plot_ice.py                   # Ice variable plotting
├── plot_fluxes.py                # Surface flux plotting
├── plot_vertical_velocity.py     # Vertical velocity plotting
├── plot_mixing.py                # Mixing coefficient plotting
└── plot_layer_thickness.py       # Layer thickness plotting
```

## Requirements

- Python 3.x
- numpy
- matplotlib
- cartopy
- netCDF4
- pyfesom2
- pyyaml

## Customization

### Plot Settings

Modify `config.yaml` plotting section:
```yaml
plotting:
  dpi: 150
  edge_threshold: 0.5
  fonts:
    title: 22
    subtitle: 18
    axis_label: 16
```

### Disable Categories

Set to `false` in config.yaml:
```yaml
categories:
  tracers: true
  ssh: false      # Skip SSH plots
  ice: true
```

## Troubleshooting

**Config not found:**
```bash
python run_analysis.py --config /full/path/to/config.yaml
```

**Script not found:**
Ensure all `plot_*.py` scripts are in the same directory as `run_analysis.py`.

**Permission denied:**
```bash
chmod +x run_analysis.py
./run_analysis.py --all
```

## Example Workflow

```bash
# 1. Copy your crash run to analysis directory
cp -r /path/to/crash_run /path/to/analysis

# 2. Edit config.yaml with correct paths
vim config.yaml

# 3. Run analysis
python run_analysis.py --all

# 4. View results
ls -lh /path/to/output/*.png
```

## Notes

- Node-based variables only (element-based fields like u, v on elements not yet supported)
- Optimized for regional plotting around crash location
- Triangulation edge masking prevents plotting artifacts
- All plots use robust percentile-based color scaling (2-98%)

## Support

For issues or questions, refer to the FESOM2 documentation or contact your local FESOM expert.
