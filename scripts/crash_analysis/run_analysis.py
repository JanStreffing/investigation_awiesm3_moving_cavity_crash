#!/usr/bin/env python3
"""
FESOM2 Crash Analysis Runner
=============================
Orchestrates plotting of blowup file variables based on config.yaml

Usage:
    python run_analysis.py [--categories cat1,cat2,...]
    python run_analysis.py --all
    python run_analysis.py --list
"""

import sys
import os
import yaml
import argparse
import subprocess
from pathlib import Path

# Available plotting scripts
PLOT_SCRIPTS = {
    'tracers': 'plot_tracers.py',
    'ssh': 'plot_ssh.py',
    'ice': 'plot_ice.py',
    'fluxes': 'plot_fluxes.py',
    'vertical_velocity': 'plot_vertical_velocity.py',
    'mixing': 'plot_mixing.py',
    'layer_thickness': 'plot_layer_thickness.py',
}

def load_config(config_path='config.yaml'):
    """Load configuration from YAML file."""
    if not os.path.exists(config_path):
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

def update_script_config(script_path, config):
    """
    Update a plotting script with configuration parameters.
    Reads the script, replaces hardcoded values with config values,
    and writes to a temporary script.
    """
    with open(script_path, 'r') as f:
        script_content = f.read()
    
    # Create temporary script with updated config
    replacements = {
        "MESH_PATH = '/work/ab0246/a270092/input/fesom2/dars2'": 
            f"MESH_PATH = '{config['mesh']['path']}'",
        
        "work_dir = '/work/bb1469/a270092/runtime/awiesm3-develop/test-spinup/run_13580101-13591231_crash/work'":
            f"work_dir = '{config['data']['work_dir']}'",
        
        "blowup_file = f'{work_dir}/fesom.1358.oce.blowup.nc'":
            f"blowup_file = f\"{{work_dir}}/{config['data']['blowup_file']}\"",
        
        "output_dir = '/work/bb1469/a270092/runtime/awiesm3-develop/test-spinup/run_13580101-13591231_crash'":
            f"output_dir = '{config['data']['output_dir']}'",
        
        "crash_lon = 29.7770298093963":
            f"crash_lon = {config['crash']['lon']}",
        
        "crash_lat = 41.2424041260004":
            f"crash_lat = {config['crash']['lat']}",
        
        "lon_min, lon_max = 25, 35":
            f"lon_min, lon_max = {config['region']['lon_min']}, {config['region']['lon_max']}",
        
        "lat_min, lat_max = 38, 44":
            f"lat_min, lat_max = {config['region']['lat_min']}, {config['region']['lat_max']}",
        
        "sys.path.insert(0, '/work/ab0246/a270092/software/pyfesom2')":
            f"sys.path.insert(0, '{config['pyfesom2']['path']}')",
        
        "dpi=150":
            f"dpi={config['plotting']['dpi']}",
        
        "edge_threshold = 0.5":
            f"edge_threshold = {config['plotting']['edge_threshold']}",
    }
    
    for old, new in replacements.items():
        script_content = script_content.replace(old, new)
    
    # Write temporary script
    temp_script = script_path.replace('.py', '_temp.py')
    with open(temp_script, 'w') as f:
        f.write(script_content)
    
    return temp_script

def run_category(category, config, script_dir):
    """Run plotting script for a specific category."""
    if category not in PLOT_SCRIPTS:
        print(f"WARNING: Unknown category '{category}', skipping")
        return False
    
    script_name = PLOT_SCRIPTS[category]
    script_path = os.path.join(script_dir, script_name)
    
    if not os.path.exists(script_path):
        print(f"WARNING: Script not found: {script_path}")
        return False
    
    print(f"\n{'='*60}")
    print(f"Running: {category}")
    print(f"{'='*60}")
    
    # Update script with config and run
    temp_script = update_script_config(script_path, config)
    
    try:
        result = subprocess.run(
            ['python3', temp_script],
            cwd=script_dir,
            capture_output=False,
            text=True
        )
        
        # Clean up temp script
        if os.path.exists(temp_script):
            os.remove(temp_script)
        
        if result.returncode == 0:
            print(f"✓ {category} completed successfully")
            return True
        else:
            print(f"✗ {category} failed with exit code {result.returncode}")
            return False
            
    except Exception as e:
        print(f"✗ {category} failed with error: {e}")
        if os.path.exists(temp_script):
            os.remove(temp_script)
        return False

def print_config_summary(config):
    """Print a summary of the configuration."""
    print("\n" + "="*60)
    print("CONFIGURATION SUMMARY")
    print("="*60)
    print(f"Mesh:        {config['mesh']['path']}")
    print(f"Work Dir:    {config['data']['work_dir']}")
    print(f"Blowup File: {config['data']['blowup_file']}")
    print(f"Output Dir:  {config['data']['output_dir']}")
    print(f"Crash:       ({config['crash']['lon']:.4f}°E, {config['crash']['lat']:.4f}°N)")
    print(f"Region:      lon=[{config['region']['lon_min']}, {config['region']['lon_max']}], "
          f"lat=[{config['region']['lat_min']}, {config['region']['lat_max']}]")
    print("="*60 + "\n")

def main():
    parser = argparse.ArgumentParser(
        description='Run FESOM2 crash analysis plotting scripts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_analysis.py --all                    # Run all categories
  python run_analysis.py --categories tracers     # Run only tracers
  python run_analysis.py --categories ssh,ice     # Run SSH and ice plots
  python run_analysis.py --list                   # List available categories
        """
    )
    
    parser.add_argument('--config', default='config.yaml',
                        help='Path to config YAML file (default: config.yaml)')
    parser.add_argument('--categories', type=str,
                        help='Comma-separated list of categories to plot')
    parser.add_argument('--all', action='store_true',
                        help='Run all enabled categories from config')
    parser.add_argument('--list', action='store_true',
                        help='List available categories and exit')
    
    args = parser.parse_args()
    
    # List categories
    if args.list:
        print("\nAvailable categories:")
        for cat in PLOT_SCRIPTS.keys():
            print(f"  - {cat}")
        print()
        return 0
    
    # Load configuration
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, args.config)
    config = load_config(config_path)
    
    # Print config summary
    print_config_summary(config)
    
    # Determine which categories to run
    if args.all:
        categories = [cat for cat, enabled in config['categories'].items() if enabled]
        print(f"Running all enabled categories: {', '.join(categories)}\n")
    elif args.categories:
        categories = [cat.strip() for cat in args.categories.split(',')]
        print(f"Running specified categories: {', '.join(categories)}\n")
    else:
        # Interactive mode - ask user
        print("No categories specified. Choose an option:")
        print("  1. Run all enabled categories")
        print("  2. Select specific categories")
        print("  3. Exit")
        
        choice = input("\nEnter choice (1-3): ").strip()
        
        if choice == '1':
            categories = [cat for cat, enabled in config['categories'].items() if enabled]
        elif choice == '2':
            print("\nAvailable categories:")
            for i, cat in enumerate(PLOT_SCRIPTS.keys(), 1):
                enabled = "✓" if config['categories'].get(cat, False) else " "
                print(f"  {i}. [{enabled}] {cat}")
            
            selection = input("\nEnter category numbers (comma-separated) or names: ").strip()
            
            # Parse selection
            if selection.replace(',', '').replace(' ', '').isdigit():
                # Numeric selection
                indices = [int(x.strip())-1 for x in selection.split(',')]
                cat_list = list(PLOT_SCRIPTS.keys())
                categories = [cat_list[i] for i in indices if 0 <= i < len(cat_list)]
            else:
                # Name selection
                categories = [cat.strip() for cat in selection.split(',')]
        else:
            print("Exiting.")
            return 0
    
    if not categories:
        print("No categories to run. Exiting.")
        return 0
    
    # Verify output directory exists
    output_dir = config['data']['output_dir']
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    # Run selected categories
    results = {}
    for category in categories:
        success = run_category(category, config, script_dir)
        results[category] = success
    
    # Print summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    
    successful = [cat for cat, success in results.items() if success]
    failed = [cat for cat, success in results.items() if not success]
    
    if successful:
        print(f"✓ Successful ({len(successful)}): {', '.join(successful)}")
    if failed:
        print(f"✗ Failed ({len(failed)}): {', '.join(failed)}")
    
    print(f"\nOutput directory: {output_dir}")
    print("="*60 + "\n")
    
    return 0 if not failed else 1

if __name__ == '__main__':
    sys.exit(main())
