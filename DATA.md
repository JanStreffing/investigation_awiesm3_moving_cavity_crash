# Data paths, experiments, and how to reproduce

All on **DKRZ Levante**, user `a270092`, project `ab0246`.

## Driver / how the runs are launched
- Runscripts: `/work/ab0246/a270092/pism_repro/scripts_is/`
  - `spinup_coupleAnt.yaml` — the iterative-coupling **driver** (model1 = awiesm3, 1-yr chunks;
    model2 = pism, 10-yr chunks).
  - `awiesm3_dyn_ocean_core3.yaml` — the awiesm3 (FESOM CORE3 + OIFS) leg runscript.
  - `spinup_pismAnt.yaml` — the Antarctic PISM (08 km) leg; `compute_time` here must cover
    PISM **+** the next chunk's `couple_in` mesh regen (raised 40 min → 2 h during this work).
- Launch / resume (from `scripts_is/`):
  `esm_runscripts spinup_coupleAnt.yaml -e <expid> --open-run`
  (`--open-run` = reuse the current dev `esm_tools` checkout; **not** a fresh install.)
- Experiments base: `/work/ab0246/a270092/pism_repro/experiments/`
- The mesh-change leg of interest is **chunk-3 = awiesm3 1901**, dir
  `run_awiesm3_19010101-19011231/`.

## Experiments used

| exp | role | key outcome |
|---|---|---|
| **movcav4** | the **evolution / hourly-movie** run (Blocker II first characterised) | FESOM CFLz blowup ~day 5–10; hourly T/S/w output → the movies |
| **movcav6** | earlier native de-containerised pipeline | OIFS **`voskin`** floating-overflow crash (Blocker I variant) |
| **movcav7** | first fully-clean native-pipeline run | chunk-3 FESOM CFLz blowup at `mstep≈205` (no voskin) |
| **movcav8** | **main step-1 analysis vehicle**; T/S-fix + 60 s tests | see below — all `restart_remapped_v{2,3,4}`, blowup file, 60 s run live here |
| **movcav10** | fresh clean run to reconfirm | OIFS **`voskin`** again at chunk-3 (voskin is state-dependent, recurs) |

### movcav8 — the analysis experiment (most files referenced in the report)
Base: `/work/ab0246/a270092/pism_repro/experiments/movcav8/`

- **Sub-meshes** (`couple/`):
  - `mother_as_old/` — pre-leg mesh (full CORE3, **211 567** nodes) = remap source geometry.
  - `latest_submesh/` — carved chunk-3 sub-mesh (**208 676** nodes) = remap target geometry.
    Files: `nod2d.out`, `elem2d.out`, `aux3d.out` (z-levels), `cavity_nlvls.out` (ulevels /
    top wet level), `nlvls.out` (bottom level), `cavity_depth@node.out` (ice draft).
- **Remapped FESOM restart** (`couple/`):
  - `restart_remapped/`      — as produced in the run (per-field `.nc`: temp, salt, w, w_impl,
    w_expl, u, v, ssh, hbar, hnode, *_AB, *_M1, ice fields).
  - `restart_remapped_v2/`   — ssh open-ocean donor + partial T/S donor.
  - `restart_remapped_v3/`   — full-column donor for fully-opened cells.
  - `restart_remapped_v4/`   — **+ shelf-base hold for still-sub-ice retreated cells → all
    static density inversions removed** (this is what the `figures/initstate/*` plots use).
  - regenerated on a **compute node** via `couple/regen_v{2,3,4}/run_remap.sbatch`.
- **What the model reads / writes** (`run_awiesm3_19010101-19011231/work/`):
  - `fesom.1900.oce.restart/` — the staged FESOM restart the compute leg ingests (the vN
    restarts were staged over this for the isolation tests).
  - `fesom.1901.oce.blowup.nc` — **the blowup snapshot** (eta_n, u, v, temp, salt, hnode, …);
    localises the 60 s `eta_n→NaN` crash.
  - `ICMGGawi3INIT` — the per-leg regenerated OIFS init (land–sea mask etc.); voskin driver.
  - `namelist.config` (`step_per_day`: 72 = 1200 s; set to **1440 = 60 s** for the timestep test).

#### movcav8 SLURM jobs (logs in `movcav8/log/`)
| job | what |
|---|---|
| `26081445` | original chunk-3 compute (pre-T/S-fix restart) — FESOM CFLz blowup |
| `26087097` | **v4 restart, 1200 s** — blew `mstep≈164`, `CFLz≈8.45` (the fast blowup; seed 90.5/−63.5) |
| `26086009` / `26086904` | v2 / v3 isolation tests — still blew (`mstep≈162/160`) |
| `26087452` | **60 s** (`step_per_day=1440`) — 0 CFLz warnings, ran to `mstep=14 629` (~10 d), then **`eta_n→NaN`** at 112.8/−65.9 |
| `26085827/26086774/26087053` | the v2/v3/v4 **remap regens** on a compute node |

### movcav4 — the movie/evolution run
Base: `/work/ab0246/a270092/pism_repro/experiments/movcav4/run_awiesm3_19010101-19011231/`.
Blowup runs: `26055476` (1200 s, m696), `26060316` (ice-seeded, m311), `26061443` (600 s, m845),
`26061910` (1200 s + hourly output, m608), `26063909` (diagnostic: hourly fw/fh/N²/Kv).
Hourly, hourly-**split** T/S/w files (`temp.1h.fesom_*`, `salt.1h`, `w.1h`) in `work/` — the
enabler for `figures/movies/*.mp4` (grab before the next leg overwrites them).

## The fix code (in `scripts/remap_and_coupling/`, from the `esm_tools` tree)
- `mo_remap_fields.F90` — the F90 restart remap. Fix commit **`ad4e84b2`** on branch
  `development_interactiv_mesh_awiesm3` (coherent-column T/S fill: open-ocean donor for
  surface-opened columns, own shelf-base hold for sub-ice retreated columns). Build with
  `couplings/fesom/remap_restart/build.sh` (**delete stale `*.mod *.o` first**); the binary
  `remap_restart` is what `couple_in` invokes.
- `coupling_ice2fesom_interactive_mesh.functions` — the ice2fesom `couple_in` glue
  (`build_submesh`, `remap_fesom_restart`, `regenerate_oasis_weights`, `remap_oasis_restart`).
- `fesom_mesh_to_cdo.py` — native pyfesom2 griddes (replaces the meshtools container).

## Environments
- Analysis / unstructured-FESOM plotting: `/home/a/a270092/.conda/envs/pyfesom2_env/bin/python`
  (numpy, scipy, netCDF4, matplotlib, pyfesom2).
- `ocp-tool` / grib+cdo / OASIS weight-gen:
  `/work/ab0246/a270092/software/miniforge3/envs/ocp-tool2/bin/python` (has `cdo`, `grib_ls`,
  `grib_get`, pyfesom2). ocp-tool dir: `/work/ab0246/a270092/software/ocp-tool/`.
- Remap build env (same compiler/netcdf as FESOM):
  `source /work/ab0246/a270092/model_codes/awiesm3-develop-is/fesom-2.7/env/levante.dkrz.de/shell`.
- Model codes: `/work/ab0246/a270092/model_codes/awiesm3-develop-is/`
  (`fesom-2.7`, `oifs-48r1`, `oasis`, `pism-github1.2.1`).

## Standalone crash-analysis suite
`scripts/crash_analysis/` (mirror of `/work/ab0246/a270092/software/fesom2_crash_analysis/`):
`run_analysis.py` + per-field plotters, driven by `config.yaml` (point it at a submesh +
`*.oce.blowup.nc`). Produced the `figures/evolution/{tracers_,verify_,*_overview}.png`.

## Reproducing the initstate figures
```
cd scripts/plotting
# edit the hard-coded C=.../movcav8/couple and R=.../restart_remapped_v4 paths at the top
/home/a/a270092/.conda/envs/pyfesom2_env/bin/python plot_layers.py   # planviews + layer-scale sections
/home/a/a270092/.conda/envs/pyfesom2_env/bin/python diag2.py         # dateline-wrapping-triangle artifact check
```
The plotters use `matplotlib.collections.PolyCollection` on the real mesh triangles and
**filter dateline-wrapping triangles** (`lon_range < 180`) — omitting that filter produces the
spurious horizontal stripes shown/explained in `figures/initstate/diag_artifact.png`.
