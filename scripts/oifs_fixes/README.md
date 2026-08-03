# The OIFS `voskin:282` fix (Bug IX, 2026-08-03)

Three commits that root-cause and fix the `voskin_mod.F90:282` single-precision floating
overflow — the "mesh-increment crash family", face (i) in the report. Written up in
`report/moving_cavity_investigation.tex`, §"Bug IX".

Source tree: `/work/ab0246/a270234/model_codes/awiesm3-develop-is/oifs-48r1`,
branch `movcav-landice+co2-concdriven`.

| patch | commit | file | what |
|---|---|---|---|
| `01-reresf_part2-setdefault-before-sugridg.patch` | `cfe38c9` | `ifs-source/arpifs/control/reresf_part2.F90` | apply `SETDEFAULT` to every block **before** `SUGRIDG` on the warm-start re-ingest, mirroring the cold-start order — stops uninitialised surface slots being serialised into the flipped cells |
| `02-surface_fields_mix-setdefault-nreqin-le-0.patch` | `8abc199` | `ifs-source/arpifs/module/surface_fields_mix.F90` | `SETDEFAULT` selects `NREQIN <= 0`, not `== -1` — seeds the nine `SD_VD` diagnostic slots it used to skip. Shared `GPOPER` code; wants review before upstream |
| `03-voskin-mask-to-open-water.patch` | `14e0fda` | `ifs-source/surf/module/voskin_mod.F90`, `surfpp_ctl_mod.F90` | optional `PFRWATER` argument; skip points with no water fraction in the four working loops — stops an ocean scheme being evaluated over land at all |

The two mechanisms are independent: 01+02 stop the garbage being created, 03 stops it being
evaluated. Either alone would have prevented the crash.

Apply with `git am` in the `oifs-48r1` tree, in order.

## Not the whole story on disk

The patches are the record; the live tree is the source of truth. Rebuild:
`libarpifs.SP.so` and `libsurf.SP.so` (both rebuilt 2026-08-03 12:27/12:30).

## Debug code left in the tree

Full sweep of the local delta against `origin/main` (17 files, 711 added lines).

**Removed 2026-08-03** (40 lines, uncommitted in the OIFS tree — the libraries are now stale
and need a rebuild):
- `VOSKIN-DIAG` probes ×2 (`voskin_mod.F90:~225`, `~249`) — from `45a7bea`; built to separate
  large `PUSTR/PVSTR` from a collapsed `PTSKM1M`. Neither was the cause and they never fired.
  With the water mask in place they could only have reported an unphysical *ocean* column.
  `ZFRWA` stays: `VOSKIN-ZLA` prints it too.
- `updtim.F90:~1059` — an orphaned six-line `--- TEMPORARY DEBUG (KSTEP=0)` comment block. The
  code it described was already gone; only the comment was left, reading as if a probe were
  still there. Pure litter.

**Keep:**
- `LLREINGEST_DBG` (`reresf_part2.F90:124`) — **not debug** despite the name; it gates the
  re-ingest feature itself. `.FALSE.` disables the moving coastline. Worth renaming.
- `LLREINGEST_VERBOSE` (`:132`) — off. Per-cell dump + bit-pattern scan of the serialised
  column. `PARAMETER`, so it compiles away entirely (verified: no `REINGEST` strings in the
  built `libarpifs.SP.so`). This is the acceptance test for the fix.
- `VOSKIN-ZLA` probe (`voskin_mod.F90:~358`) — the tripwire on the actual failure mode. Covers
  the one route the mask does not close: broken Stokes at a genuine ocean point (e.g. WAM).
  If it goes, `ZFRWA` goes with it.
- `RERESF_PART2/LANDICE` re-ingested-cell-count line, and the `XIOSFPOS` axis-name echo in
  `suxios.F90` — provenance and config logging, not debug.
- The two `ABOR1` guards on `LGPOWMASK` in `GPOPER`'s `PUTALLFLDS_M` branch — assertions on a
  masked write.
- `updtim.F90:~885` `WRITE(*)` duplicating the `NULOUT` reseed count — deliberate: `NULOUT` is
  buffered and lost if a later routine aborts first. Drop once legs stop aborting.

Everything else in the delta (`surfece`, `ecearth`, `callpar`, `srfi/srfis_mod`, the `ece_*`
coupling routines, `cplng2_data_mod`, `yomxios`) is functional — no log-only code.
