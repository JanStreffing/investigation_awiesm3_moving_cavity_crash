# AWI-ESM3 ↔ PISM moving-cavity crash investigation

Investigation log, plots, and scripts for two coastal-margin blockers hit on the first
mesh-change leg of an end-to-end **AWI-ESM3 (FESOM2 + OIFS-48r1 + LPJ-GUESS) ↔ Antarctic-PISM
moving-cavity** coupled run (driven by `esm_tools`). On each mesh-change leg the coupling
machinery regenerates the FESOM sub-mesh, the OASIS grids/masks/weights, and the OIFS `ICMGG`
so all components agree on the new coastline.

## The two blockers

- **Blocker I — OIFS step-0 coastline crash — RESOLVED.** A floating-point NaN struck exactly
  the coastline cells whose land–sea type flipped: `ocp-tool` half-converted a flipped cell
  (only mask + soil type), OIFS re-ingested the inconsistency from the regenerated `ICMGG` and
  NaN'd in moist physics. Fixed by a masked NN rebuild of each flipped cell + `ICMGG` re-ingest.
  (A recurring `voskin` floating-overflow variant still shows up state-dependently — see the
  report / `DATA.md`.)

- **Blocker II — FESOM cavity-margin vertical-CFL blowup — LIVE, now root-caused.** ~5–10 model
  days in, FESOM blows up with a vertical-CFL runaway ringing the cavity margin. **Root cause: a
  geostrophically-unbalanced remapped state** — the mesh change plants fresh cavity fronts
  against salty open water at the moving margin but *zeroes* the balancing current, so the
  model's `t=0` adjustment grows `w` at the fronts. Every restart-*value* fix (freezing cap;
  coherent T/S refill to zero density inversions) was tested and does **not** prevent it; 60 s
  clears the *fast* vertical-CFL mode but a *slow* barotropic `eta_n` mode then grows over
  ~10 days. **Proposed fix: initialise the changed-cell velocity in thermal-wind balance instead
  of zero.**

## Repository layout

```
report/            movcav_lsm_investigation.{tex,pdf}  — the full investigation log (read first)
figures/
  initstate/       step-1 (t=0 remapped state) plots — the root-cause evidence (movcav8 v4)
  evolution/       hourly-movie key frames + earlier crash-analysis overviews (movcav4)
  movies/          blowup_TSw_evolution.mp4, blowup_section_evolution.mp4 (movcav4)
scripts/
  plotting/        the initstate PolyCollection plotters + the dateline-artifact diagnostics
  crash_analysis/  the standalone FESOM2 blowup-file analysis suite (run_analysis.py + plots)
  remap_and_coupling/  the actual fix code: F90 restart remap, pyfesom2 griddes, couple glue
  flux_decomp.py   OASIS heat-flux decomposition (solar vs non-solar)
DATA.md            experiments, data paths, key files, jobs, analysis env — start here to reproduce
```

## Key figures
- `figures/initstate/planview_seed_weddell.png` + `section_lat_seed_weddell.png` — the
  **initiation seed**: fresh cavity column (S≈31) at the ice front, **zero velocity**, `w` dipole.
- `figures/initstate/planview_A_90E_L38.png`, `planview_D_60s_etacrash.png` — two other
  high-CFLz sites, both **benign at t=0** (peak-CFLz is a *consequence*, not the seed).
- `figures/initstate/diag_artifact.png` — proof the "horizontal stripes" in early plots were
  **dateline-wrapping triangles** (a plotting artifact), not data.
- `figures/movies/*.mp4` — the blowup growing over days (the evolution view, reconciled in the
  report with the step-1 view).

## Build the report
```
cd report && pdflatex movcav_lsm_investigation.tex    # figure paths are relative to report/? no — see note
```
Note: the `.tex` uses `\includegraphics{movcav4_crash_plots/...}` paths from the original
working tree; the packaged `.pdf` is the built version. To rebuild against this repo, point the
graphics path at `../figures` (the `initstate/`, `evolution/` subdirs match) or keep the shipped PDF.

Analysis Python (unstructured FESOM plotting): `/home/a/a270092/.conda/envs/pyfesom2_env/bin/python`.
