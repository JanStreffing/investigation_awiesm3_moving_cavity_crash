# Leg-start crash at the moving cavity front

Working notes — fgdbg02 cycle-7 (1906) step-1 crash. Status: **the SSH-RHS framing is
falsified** (§9). The crash is a vertical-column failure at five specific nodes, seeded by
exact zeros written into `u` (not `v`) at levels 4–5 of three elements during step 1 — a
signature that appears nowhere else in space or in any saved state of the campaign (§12, §13).
See §9–§14; §3, §5 and §8 above them are kept for the record but their conclusions no longer
hold.

## 1. The crash

- Run `fgdbg02` (`which_ale='linfs'`, `step_per_day=72` = 1200 s, fesom omp=1 →
  bit-reproducible, 100-yr target). Crashed at **cycle 7 (year 1906), mstep = 1**.
- Location: **lon −100.38, lat −72.86 — Amundsen Sea, Pine Island / Thwaites ice front**,
  West Antarctica. Several adjacent nodes blow up together (1606/1603…).
- Blow-up detector fires on temperature: `temp → −3.7e243` in one step (full numerical
  explosion, top ~17 levels). `hflux = −90.5 W/m²` (normal), `[FLUXGUARD]` **silent** →
  not flux-driven.
- Log: `run_awiesm3_19060101-19061231` →
  `log/fgdbg02_awiesm3_compute_19060101-19061231_26566252.log` (blow-up dump ~line 19881).

Crash-node dump (node is **open ocean**, `ulevels=1`, `nzmin,nzmax = 1,18`):
```
eta_n      = -1.557      hbar_old = -1.554   (depression already in the restart)
ssh_rhs    =  7,142,987  ssh_rhs_old = -208  (the one exploded field)
d_eta      = -0.0025     (tiny -> not through the SSH increment)
hnode      = 5,5,10,...,50,0,0...   (normal)   CFL_z all < 0.07   Kv/W normal
temp_old   = 2.81,2.73,2.60,1.54,0.28,-0.39,... (physical incoming state)
```

## 2. Not the zstar seam

Earlier the report draft mis-attributed this to the `zstar` cavity-seam collapse of
`sec:linfs-coupled`. **Wrong** — the run is on `which_ale='linfs'` (verified in the staged
`namelist.oce`). `linfs` is active and does not cure this. It is a genuinely different face.

## 3. Mechanism chain

`ssh_rhs` is the ALE external-mode RHS = vertically-integrated **transport divergence**
(`compute_ssh_rhs_ale`, `oce_ale.F90:1818`):
```
ssh_rhs(node) = Σ_edges α·(UV + UV_rhs)·n̂ · helem   + (1-α)·ssh_rhs_old
```
with `UV_rhs = ab1·UV_rhsAB + fresh(−g∇η − ∇p/ρ, Coriolis…)` built in
`oce_ale_vel_rhs.F90:120,155`. `ssh_rhs_old = −208` is normal, so the 7.1e6 is the
**transport-divergence term**. Because `d_eta` stays tiny, the blow-up is the divergence
itself, not the SSH increment. `temp` then explodes via advection by the diverged velocity;
the temp dump is post-explosion, so `ssh_rhs=7.1e6` is the symptom that survives to the dump.

## 4. What we ruled out (all measured, all clean)

| suspect | finding | verdict |
|---|---|---|
| Draft shock magnitude | cycle-7 couple-in is the **gentlest** of the series (max\|Δdraft\| 479 m vs 519–587; fewest nodes >40 m; fewest nlvls changes) | not the trigger |
| Mesh degeneracy | crash box: no new nodes, min-angle 46–56°, uniform areas (1.06–1.64e8) | clean |
| T/S seeding | deepened columns have smooth physical cavity profiles; horiz ΔT to open node ≤ 0.37 °C | clean |
| Baroclinic PGF from density | column-integrated proxy ≈ **0.001 m/s** over 1200 s | negligible |
| Adams-Bashforth history | `urhs_AB` accel at crash region ~1.3e-5, **below** global 99.9th pct; \|urhs_AB\| max ~3e5 in **every** restart incl. healthy chunk-1 | not anomalous |
| Velocity `u/v` | ~0.1 m/s at crash region (global max ~1.9) | physical |
| ssh=0 vs −1.55 "jump" | the ssh=0 nodes are **cavity** nodes (ssh=0 under ice) — normal cavity-front config, ~2400–3355 every cycle | normal |

## 5. The mechanism we DID find: pinched cavity-front peninsula tip

Element upper level = **max** of its node cavity levels (confirmed: `cave==max(node)` for
404 359/404 607 elements). The crash node (open, `ul=1`) is a **peninsula tip**: its
6 surrounding elements, 4 touch the two neighbours that just deepened `ul3→ul7`, so those
4 elements become `cav_elvls=7` (ice at levels 1–6). Only 2 elements stay open at the
surface. Its upper-level control volume, previously bounded by 6 triangles, is now bounded
by **2** — a dead-end water pocket.

Per-level wet-element count at the crash node, by cycle:
```
1901: [3,3,3,6,6,6,6,6]
1902: [2,2,2,2,2,2,2,6]   levels 1-7 trapped ... and SURVIVED
1903: [2,2,6,6,6,6,6,6]
1904: [2,2,6,6,6,6,6,6]
1905: [2,2,6,6,6,6,6,6]
1906: [2,2,2,2,2,2,6,6]   levels 3-6 newly trapped -> CRASH
```
The remapped velocity (near-divergence-free on the OLD 6-connected geometry) now has
nowhere to go through the pinched faces → `ssh_rhs = flux/(tiny pocket volume)` blows up
on step 1.

## 6. But it is ubiquitous, not rare — the open question

The pinch is the **stage**, not the deterministic **trigger**. The trigger event
(an open node gaining newly-trapped upper levels) happens **every cycle, in the dozens**:

| transition | nodes gaining a trapped level | ≥3 levels | ≥4 (as severe as crash) | worst |
|---|---|---|---|---|
| 1901→02 | 148 | 89 | 57 | 10 |
| 1902→03 | 142 | 75 | 53 | 10 |
| 1903→04 | 137 | 60 | 47 | 10 |
| 1904→05 | 123 | 65 | 46 | 10 |
| 1905→06 | 104 | 55 | **36** | 10 |

The crash node's own event (4 trapped levels) is unremarkable — dozens of nodes do the
same or worse every cycle, and **1902 trapped levels 1–7 at this very node and survived.**
Across ~7 years there have been **many hundreds of equally-severe pinch events and exactly
one blew up** (~1-in-hundreds).

**Implication:** a blanket "zero velocity at all pinched peninsula tips" fix would perturb
40–60 healthy cavity-front nodes every cycle and risks the validated melt/current physics.
Too blunt for a 1-in-hundreds event. We need the **extra ingredient** that tips one severe
pinch into a blow-up — most likely the actual transport magnitude in the trapped cells at
that instant.

## 7. Constraints on the fix (from operator)

- **No dt reduction** (leg-start timestep cushion is off the table).
- **No limiting** of the cavity draft change per cycle (would make PISM and AWI-ESM drift
  apart — unphysical).
- The fix belongs in the **restart modification for the AWI-ESM3 leg based on PISM** —
  i.e. `couplings/fesom/remap_restart/mo_remap_fields.F90` (`remap_fesom_restart`).

## 8. Next step

Offline discriminator: compare the crash node against the ~35 *other* cycle-7 nodes that
got an equally-severe pinch but survived. If the crash node has a far larger trapped-cell
transport-divergence proxy (flux through open faces ÷ pocket volume) than every survivor,
a **narrow safety valve clamping only anomalous divergence** in the restart step is
justified — without touching the 99.7 % healthy cases. If not, instrument step 1
(bit-reproducible, immediate crash) to capture the PGF/transport directly.

## 9. 2026-07-31: the SSH-RHS framing is falsified

The 2 GB **`fesom.1906.oce.blowup.nc`** in the crash leg's `work/` — a full 3D dump of
`ssh_rhs`, `d_eta`, `u`, `v`, `u_rhs`, `v_rhs`, `urhs_AB`, `helem`, `hnode`, `w/w_expl/w_impl`,
`cfl_z`, `temp`, `salt`, `Kv`, `N2` at the crash step — was never opened. It settles §3/§5/§8.

~~`ssh_rhs = 7.1e6` is the exploded field.~~ It is an ordinary value for this configuration:

| in the blowup dump | count |
|---|---|
| nodes with \|ssh_rhs\| > 1e6 | 1461 (1458 cavity-front-adjacent) |
| nodes with \|ssh_rhs\| > 1e7 | 506 |
| global max \|ssh_rhs\| | 5.97e7, **8× the crash node** |
| crash node's rank among the 184 pinched open nodes | **143rd of 184** |

and `d_eta` stays in [−1.76, +0.07] m everywhere — the CG solve absorbs the large rhs. There
is no SSH-solver explosion. `ssh_rhs = 7.1e6` with `d_eta = −0.0025` is simply what a
cavity-front node looks like every step of every year.

~~The pinched peninsula tip drives the blow-up.~~ The pinch is real but it is the normal
front geometry, and its ssh_rhs signature is unremarkable at the crash node.

## 10. Why the front carries a huge `UV_rhs` at all (structural, not the trigger)

`ssh_rhs` uses `(UV + UV_rhs)`, and it is `UV_rhs` that is large: at the crash node's four
cavity-covered elements `max|u_rhs|` = 1.11–1.25 m/s, depth-uniform over levels 7–17, while
`u` itself is 0.02–0.03 m/s. Depth-uniform ⇒ a surface-pressure term. The source:

- `eta ≡ 0` at **all 2416 cavity nodes**, exactly, while open ocean south of 60S sits at
  **−1.58 m** (global mean eta −0.47 m).
- This is **FESOM's own convention, not a remap artefact** — the native end-of-year restarts
  for 1900–1905 all have eta identically 0 under the cavity.
- So each of the **1311 mixed cavity/open elements** carries a ~1.59 m eta step (median 1.593,
  max 1.704, all >1.0 m). `compute_vel_rhs` uses `p_eta = g*eta_n(elnodes)` over the three
  nodes ⇒ `-g∇η·dt` ≈ 9.81 × 1.59/12000 × 1200 = **1.56 m/s**. The baroclinic `pgf` cancels
  only ~20 % of it.

Census over element classes (blowup, `max|u_rhs|` over wet levels):

| class | n | p50 | p90 | p99 | max |
|---|---|---|---|---|---|
| front (mixed cavity/open) | 1311 | **1.011** | 1.204 | 1.444 | 1.787 |
| full cavity | 400 (sample) | 0.006 | 0.016 | 0.038 | 0.052 |
| open ocean | 400 (sample) | 0.007 | 0.021 | 0.035 | 0.065 |

The crash elements (1.11–1.25) sit between p50 and p90 of the front population. **They are
typical front elements.** So the §4 "baroclinic PGF ≈ 0.001 m/s, negligible" measurement was
right but irrelevant — the term that matters is the barotropic `-g∇η` across the eta seam, and
it was never measured.

## 11. The actual needle: five columns with an interior `W == 0`

Scanning every node column for a zero of `w` strictly inside `[ulevels, nlevels)`:

```
crash step (blowup):                    5 nodes
end-1900 .. end-1905 restarts:          0 nodes   (~1.25 M columns)
```

and the five are exactly the crash cluster:

```
19342  -100.383 -72.864  ul=1 nl=18   W==0 at level 5          <- log dump, mype 1606
19343  -100.112 -72.711  ul=1 nl=18   W==0 at levels 1-5
19344  -100.564 -72.592  ul=1 nl=19   W==0 at levels 1-5
19345  -100.719 -72.739  ul=1 nl=20   W==0 at level 5
19357  -100.019 -72.556  ul=1 nl=16   W==0 at levels 1-4       <- log dump, mype 1603
```

`w_expl` carries the same zeros, `w_impl` is identically zero (no CFL split), `cfl_z < 0.07`.
Call order in `oce_timestep_ale` is `vert_vel_ale` → `solve_tracers_ale` →
`update_thickness_ale` → `check_blowup`, so **W is computed before the tracer solve in the
same step**: the holes are upstream of the tracer explosion, not a consequence of it.

The tracer signature is a **near-singular vertical tridiagonal solve**, not advection: `temp`
decays geometrically by ~10³ per level away from a single peak level (nz=5 at 19342, nz=4 at
19343/19357, peak ~1e250), and `salt` sits at exactly 3.0 / 45.0 throughout the column.

Falsified along the way: the "orphan level" idea (a node open at the top with every element
starting below) — node 19344 has **all six** elements open (`ul=1`) and still has W==0 at
levels 1–5. `cavity_elvls.out` equals `max(node ulevels)` for 404 359/404 607 elements and
none of the 248 exceptions is in the crash cluster.

## 12. 2026-07-31: `area(nz,n)` is clean — the fault is in `u`

Built `fesom_meshdiag` (a standalone utility already in the tree, `BUILD_MESHDIAG=OFF` in the
production build) with a `dump_areas` routine that gathers `area(nz,n)`, `areasvol(nz,n)`,
the level bookkeeping and `elem_area` to rank 0. Ran it on the 1906 submesh over the existing
`dist_1792` partition (job 26596278, 14 nodes, 17 s).

~~`area(nz,n)` is degenerate at the five columns.~~ **It is not.** Across all 208 109 columns:

```
area     == 0 inside [ulevels, nlevels):  0     < 1 m^2:  0
areasvol == 0 inside [ulevels, nlevels):  0     < 1 m^2:  0
ulevels_nod2D == cavity_nlvls.out for 208109/208109 nodes
ulevels(elem) == cavity_elvls.out for 404607/404607 elements
```

The five W-hole columns have entirely ordinary areas (0.8–2.9e8 m², stepping down with depth
exactly as their neighbours do), and at node 19345 the area is *constant* over levels 1–13
while W(5) is still exactly 0. So the control volume is not the fault.

The dump did fix the element level bounds, though — `nlevels(elem)` is **not** min-over-nodes,
and using the model's own values instead of that guess turned 23 780 spurious hits into a
clean signal:

| interior exact zeros, using the model's `ulevels`/`nlevels` | elements |
|---|---|
| `u`, blowup (post-step-1) | **3** |
| `v`, blowup | 0 |
| `helem`, blowup | 0 |
| `u` and `v`, incoming remapped restart | 39 (identical levels in both components) |

The three are **39028, 39029, 39030** — and they touch exactly the five W-hole nodes,
no others. Their zeros are at **levels 4 and 5, in the zonal component only**:

```
elem 39030 (ul=1 nl=5)   u = 0.06535  0.05930  0.01744  0.       (0)
           restart       u = 0.06534  0.05816  0.01206 -0.01093
                         v = -0.02459 -0.00033 -0.03234 -0.01747   <- v updated normally
elem 39029 (ul=1 nl=14)  u = 0.10995  0.10316  0.05720  0.  0.  0.02065 ...
           restart       u = 0.11033  0.10240  0.05206  0.02921  0.02091  0.01755 ...
elem 39028 (ul=1 nl=18)  u = 0.06315  0.05752  0.00809  0.  0. -0.02747 ...
           restart       u = 0.06592  0.05865  0.00531 -0.01826 -0.02588 -0.02902 ...
```

`helem` at those levels is 10 m, `v` updates normally, and the restart had **no** zeros on
these three elements. So the zeros are *created by the model during step 1*, in one velocity
component, on three consecutive elements.

The 39 restart elements are a separate, benign population: `u` and `v` zero at *identical*
levels, always a prefix 1..4 or 1..5, located in the Ross (~-83°) and Ronne (~-78°) sectors,
touching none of the five crash nodes. That is what a retreating-cavity remap leaves behind
and the model has carried it for years without trouble.

## 13. Is it unique in space and time? Yes, in every state we can observe

`cavity_elvls.out` == model `ulevels(elem)` and `elvls.out` == model `nlevels(elem)`, exactly,
for all 404 607 elements (checked against the meshdiag dump), so the correct bounds are
available for every submesh without rerunning the utility. Two populations, kept apart:

- **BOTH** — `u` and `v` zero at *identical* levels. Retreating-cavity remap residue.
- **U-ONLY** — `u` zero where `v` is not. The step-1 signature.

| state | BOTH | U-ONLY |
|---|---|---|
| 1899 spin-up input (cold start, u=v=0 everywhere) | 405 888 (all) | 0 |
| 1900 native end-of-leg | 0 | 0 |
| 1900→01 incoming remapped | 241 | 0 |
| 1901 native end-of-leg | 0 | 0 |
| 1901→02 incoming | 23 | 0 |
| 1902 native end-of-leg | 0 | 0 |
| 1902→03 incoming | 13 | 0 |
| 1903 native end-of-leg | 0 | 0 |
| 1903→04 incoming | 2 | 0 |
| 1904 native end-of-leg | 0 | 0 |
| 1904→05 incoming | 27 | 0 |
| 1905 native end-of-leg | 0 | 0 |
| 1905→06 incoming (the crash leg's input) | 39 | 0 |
| **crash step 1 (blowup dump)** | **0** | **3** |

Three things fall out:

1. **U-ONLY is zero in all thirteen observable states across seven cycles**, and 3 at the one
   mid-leg instant we can see. In space it is 3 elements of 404 607; in time it is the only
   occurrence anywhere in the campaign's saved states.
2. **The remap's BOTH residue is normally healed.** Every incoming remapped restart carries
   some (2–241), every native end-of-leg restart has none, and at the crash step it is already
   0 — so step 1 fills those zero velocities routinely.
3. **The three U-ONLY elements are not among the 39.** They had no zeros at all in the incoming
   restart. So this is not "the remap residue failed to heal"; it is three elements that were
   fine on input and came out of step 1 with the zonal component zeroed.

Caveat on "time": these are leg-start inputs and end-of-leg outputs, plus the single crash
step. There are no mid-leg snapshots, so this does not prove a U-ONLY zero never appears
transiently inside a healthy year — only that it never survives into any saved state, and that
it is present at the one mid-leg instant that was ever dumped.

## 14. Can the fix live in the remap tool? Not yet justifiable — no input-side handle

§7 puts the fix in `mo_remap_fields.F90`. Two findings bear on that, and both are negative
for a remap-side fix *as currently understood*:

1. **The remap does not write the bad values.** Elements 39028/39029/39030 have no interior
   zeros in `u` or `v` in the incoming restart. The zeros appear only after step 1 (§12).
   The 39 elements the remap *does* leave zeroed are elsewhere and are healed routinely (§13).
2. **The remap does not touch these elements at all.** Comparing the 1905 and 1906 submeshes
   through the stable `map_elem.out` / `map_nod.out` IDs:

```
elem 39028 (id 39720)  1905 ul=1 nl=18  ->  1906 ul=1 nl=18   unchanged
elem 39029 (id 39721)  1905 ul=1 nl=14  ->  1906 ul=1 nl=14   unchanged
elem 39030 (id 39722)  1905 ul=1 nl=5   ->  1906 ul=1 nl=5    unchanged
   and all six distinct nodes (19342-19345, 19357) unchanged in both ulevels and nlevels
```

None is new (129 elements are), none had `nlevels` changed (30 did), none had `ulevels`
changed (1913 did). Geometrically these three elements and their nodes are **identical** to
the previous cycle. What changed is their *neighbourhood* — adjacent elements that did flip
`ulevels` at the couple-in.

So a remap-side fix has no property of these elements to key on. It would have to act on the
neighbourhood of a `ulevels` change, which is the blunt 40–60-nodes-per-cycle intervention
§6 already rejected — unless the instrumented rerun shows the trigger is an input-state
property the remap can precondition. That is the gate.

Note this also further undercuts §5: the elements that acquire the fault are not the pinched
ones. The four elements that flipped `ul3→ul7` at the crash node (39021/39022/39025/39026) are
*not* the three that come out of step 1 with zeroed `u`.

## 15. `[UVZERO]` — the detector is cheap enough to run online

The predicate: **exactly one of `u`,`v` is exactly `0.0` at a level inside
`[ulevels(elem), nlevels(elem)-1]`.** Both components zero together is *legitimate* — that is
what the cavity remap leaves on newly-opened columns and step 1 fills it in normally (§13) —
so only the asymmetric case is flagged.

The wet range must be the **full** one. 39030's zero sits at its bottom wet level, so a
variant that excludes `nlevels-1` catches only 2 of the 3.

| | |
|---|---|
| scan cost, whole global mesh (numpy, vectorised) | 37 ms |
| comparisons per step, globally | 12.1 M |
| per rank at 1792 PEs | ~6 800 |
| false positives over 13 saved states / 7 cycles | **0** |
| hits at the crash step | 3 (39028, 39029, 39030) |

Negligible next to `check_blowup`, which already sweeps every node and element every step.
Implemented as a `[UVZERO]` block in `check_blowup` (`write_step_info.F90`), next to
`[FLUXGUARD]`: sets `found_blowup_loc`, dumps the element's `u`, `v`, `u_rhs`, `v_rhs`,
`helem` columns and its three nodes, then takes the normal blow-up exit — so it fires
*before* the tracer explosion and writes a blowup file with the velocity fault still visible
rather than 1e250 tracers on top of it.

Compile-checked in the uncoupled configuration only. It sits outside the `__oasis` guard so
it compiles in both, but confirming the coupled build means relinking the production
`libfesom.so`/`fesom.x` — not done.

Caveat: this is a **detector, not a fix**. It converts a silent corruption into a labelled,
early abort and would catch a recurrence anywhere in the domain on any future leg. It does
not prevent the zero.

## 16. The values were overwritten, not computed — and what that does NOT yet tell us

`update_vel` computes `UV = UV_old + UV_rhs + Fx`. With `UV_old` from the restart, `UV_rhs`
from the blowup and `Fx` bounded by ~0.0013 m/s from the actual `d_eta` at those nodes:

| element | level | should be | model has |
|---|---|---|---|
| 39028 | 3 | +0.0062 | +0.0081 ✓ |
| 39028 | **4** | **−0.0179** | **0.00000** |
| 39028 | **5** | **−0.0260** | **0.00000** |
| 39029 | 3 | +0.0543 | +0.0572 ✓ |
| 39029 | **4** | **+0.0302** | **0.00000** |
| 39029 | **5** | **+0.0213** | **0.00000** |
| 39030 | 3 | +0.0152 | +0.0174 ✓ |
| 39030 | **4** | **−0.0100** | **0.00000** |

Level 3 is right; levels 4–5 of the same element, same loop iteration, same statement, same
`Fx`, are off by 0.010–0.030 m/s — 8–22× the largest `Fx` can be. The values were overwritten
after being computed correctly.

**Not the atmosphere.** `tx_sur` = −0.013…−0.015, `ty_sur` ≈ −0.034 at the three elements
against a global median |`tx_sur`| of 0.043 — unremarkable. And `stress_surf` enters momentum
only at `nzmin`; levels 4–5 are interior, with no atmospheric path to them.

**Operator's position (2026-07-31), and it is well founded:** a decade of FESOM without this
failure mode, two weeks of moving cavity and it appears. The moving-cavity chain is causally
involved and the fix belongs there. That is compatible with the measurement above — a state
handed to FESOM that makes it index or communicate wrongly produces exactly this.

### Moving-cavity chain outputs checked offline, all clean

| candidate | result |
|---|---|
| element/node `ulevels`,`nlevels` at 39028–39030 | unchanged 1905→1906, none new (§14) |
| `area(nz,n)` / `areasvol(nz,n)` | no zeros, nothing <1 m², whole mesh (§12) |
| `cavity_elvls.out` vs model `ulevels(elem)` | identical, 404 607/404 607 |
| `elvls.out` vs model `nlevels(elem)` | identical, 404 607/404 607 |
| `u`,`v` in the remapped restart at those elements | no zeros (§13) |
| node ownership in `dist_1792` | unique, 0 duplicates, all 6 cycles |
| element "duplicate ownership" | **not a defect** — ~83 000 elements (20 %) are multi-owned in *every* cycle including healthy ones; element sharing is normal in FESOM2 |
| halo comm arrays (`com_info*`), send/recv symmetry | 0 mismatches, nodes and elements, all 6 cycles |

So every *static* output of the moving-cavity chain that can be checked offline is
self-consistent. What cannot be checked offline is the runtime behaviour, and that is the only
place left for the trigger to hide.

### The run that would settle it

A **bounds-checked build** (`-check bounds`, or `-check all`) of the 1906 leg. The leg is
bit-reproducible under `omp=1` and dies on step 1, so it is one short run. If anything writes
outside an array, that build traps it at the exact source line with the array name and the
offending index. If the culprit sits in code that only executes when the cavity geometry
changed, that is the moving-cavity fix site, named directly rather than inferred — and no
FESOM-side repair is needed at all.

Second choice if the bounds build comes back clean: print `UV(1,4:5,39028:39030)` after each
of `compute_vel_rhs`, `viscosity_filter`, `impl_vert_visc_ale`, `update_vel` and
`exchange_elem`, which brackets the write to a single call.

## 17. 2026-07-31 rerun: the write happens in `solve_tracers_ale` — the chain was backwards

Instrumented rerun of the 1906 leg (job 26599924, same script, bit-reproducible). A
`uvzero_probe` — reporting elements where exactly one of `u`,`v` is exactly 0.0 inside the wet
range — was called at seven points of the timestep. Result:

```
after compute_vel_rhs        silent
after viscosity_filter       silent
after impl_vert_visc_ale     silent
inside update_vel pre-exch   silent
inside update_vel post-exch  silent      <- exchange_elem is NOT the writer
after compute_hbar_ale       silent
after vert_vel_ale           silent
after solve_tracers_ale      FIRES: elem 2,3,4 (mype 1603) / 8,9 (mype 1606), nz=4,5, u==0
```

`solve_tracers_ale` has no business writing to `dynamics%uv` at all. The zeros are written
there.

**This inverts §11/§16.** The chain recorded there was: `u` zeros → corrupted divergence →
`W == 0` → singular vertical tracer operator → `temp` ~1e250. But `vert_vel_ale`, which
computes `W`, runs *before* `solve_tracers_ale`, and the probe after it is silent — so `W` was
computed from a clean `UV`. The `W == 0` holes seen in the blowup dump (§11) must therefore
have been written later, during the same tracer solve. ~~The u zeros seed the tracer
explosion.~~ **Both the `u` zeros and the `W` zeros are collateral from whatever goes wrong
inside `solve_tracers_ale`.** The half-zero state remains a valid and uniquely specific marker
(§13) — it is just a symptom, not the seed.

Both detectors fired on the rerun: 86 temperature blow-ups and 9 `[UVZERO]`, all at mstep 1,
confirming the failure reproduces exactly and that the guard catches it.

Still standing from earlier: `area`/`areasvol` clean (§12), the remap does not produce the
state and does not touch those elements (§14), the partition and its comm arrays are
self-consistent (§16).

**Next:** the write is inside one routine, so a bounds-checked build should name the line.
Built via a copy of esm_tools' own `comp-fesom-2.7-main_script.sh` with
`-check bounds -check pointers -check uninit -traceback -g` added through
`CMAKE_Fortran_FLAGS` (`esm_master recomp` itself is broken on this tree by an unrelated
`debm` section error). If the trap lands in code that only executes when the cavity geometry
changed, that is the moving-cavity fix site named directly.

## 18. SOLVED — an unassigned automatic-array element in the GM streamfunction solver

**Root cause.** `fer_solve_Gamma` (`oce_fer_gm.F90`) solves the GM streamfunction over

```
nzmin = ulevels_nod2D_max(n)   ! max over the node's surrounding elements
nzmax = nlevels_nod2D_min(n)   ! min over the node's surrounding elements
```

At node **19343** these are **7 and 5** — inverted. Both sweeps then execute zero times
(`DO nz = nzmin+1, nzmax` = 8..5; `do nz = nzmax-1, nzmin, -1` = 4..7), but this line runs
unconditionally:

```fortran
tr(:,nzmax) = tp(:,nzmax)      ! tp(:,5) was never assigned
```

`tp` is a local automatic array, so an uninitialised stack word is published straight into
`fer_gamma(1,5,19343)` ≈ 4.9e252.

**Measured chain** (probe output, job 26601651):

```
fer_uv(1, nz=4, elem 2) = -1.63611672e+251     fer_uv(2,...) = +1.8e-04   helem = 10.0
fer_uv(1, nz=5, elem 2) = +1.63611672e+251     fer_uv(2,...) = +3.6e-06   helem = 10.0
```

`fer_gamma2vel` forms `fer_uv = Σ(gamma(nz) − gamma(nz+1))/(3·helem)`, so one bad gamma at the
level-5 interface gives that antisymmetric pair, zonal component only, `helem` healthy.
`solve_tracers_ale` then adds `fer_uv` to `UV` (line 206), **advects the tracers with it** —
`temp` → 1e250, `salt` pinned at 3/45 — and subtracts it back (line 339). `(x+y)−y` with
`|y| ≫ |x|` returns exactly `0.0`, which is the half-zero `u`. The identical add/subtract pair
is applied to `Wvel`/`Wvel_e` with `fer_Wvel` five lines below: that is the `W == 0` holes of
§11.

**The fix** (at the defect, not a workaround):

```fortran
if (nzmin >= nzmax) then
    tr(:,:) = 0.0_WP
    cycle
end if
```

A node with no interior range carries no GM transport. Verified: the leg previously died at
`mstep 1` every single run; with the fix it reached **day 27 / step 1920** with `UVZERO=0`,
`STOP=0`, no blow-up (job 26601885, cancelled once past the failure point).

**Why it is a moving-cavity failure, and why it took a decade to appear.** The trigger needs a
node touching both a deep-cavity element (large `ulevels`) and a very shallow one (small
`nlevels`) — a cavity front on steep bathymetry. The predicate
`max(ulevels(elems at n)) > min(nlevels(elems at n))` finds:

| mesh | inverted nodes |
|---|---|
| submesh 1901 | 1 |
| submesh 1902 | 2 |
| submesh 1903 | 1 |
| submesh 1904 | 2 |
| submesh 1905 | 0 |
| submesh 1906 | 2 (incl. **19343**, nzmin=7 nzmax=5) |
| base pre-coupling mesh | 2 |

So the configuration recurs but is not always fatal: `tp(:,nzmax)` is whatever the stack holds
— usually benign, occasionally 1e252. That is the "1-in-hundreds" behaviour of §6, it is why
the bounds-checked build passed with identical source (§17), and it is why `use_cavity` + GM is
required to see it at all. `Fer_GM = .true.` and `redi = .true.` in this configuration.

**Caveats.** The base mesh already contains such nodes, so the latent defect is not created by
the coupling workflow — it is exposed by it. Earlier legs ran with this active; wherever an
inverted node existed and `tp` was small-but-nonzero, GM transport there was junk without
crashing, so those years are suspect.

**Two unrelated out-of-bounds reads found on the way**, both from Fortran's lack of guaranteed
short-circuit `.and.`, both fixed:

- `io_xios.F90:530` — `k >= 1 .and. angles(k) > tmp_ang` reads `angles(0)` at the end of every
  insertion sort. Benign, but aborts any bounds-checked build during init.
- `ice_maEVP.F90:861` and `:1255` — `edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1` reads
  `ulevels(0)` at every boundary edge. Inside `if (use_cavity)`: the garbage value decides
  whether sea-ice velocity is zeroed at the cavity–ocean edge. Not the crash (fixing it left
  the failure bit-identical) but a real bug in cavity physics.

## 19. Where this leaves it

The chain now reads: something in step 1's momentum path writes exact zeros into `u` (not `v`)
at levels 4–5 of elements 39028–39030 → the divergence those elements feed is wrong → `W`
comes out exactly 0 at the five nodes they span → the vertical tracer operator on those five
columns is singular → `temp` ~1e250 and `salt` pinned at 3/45 → `check_blowup`.

Two readings, not yet separated:

1. **A dynamical zero.** Some branch in `compute_vel_rhs` / `viscosity_filter` /
   `impl_vert_visc_ale` / `update_vel` legitimately assigns 0 to `u` on these levels. Then the
   question is which branch and why only the zonal component.
2. **A scribble.** The affected indices are contiguous (elements 39028–39030, nodes
   19342–19345 plus 19357), the values are exact zeros rather than small numbers, only one
   component of one array is hit, and the build uses `-init=zero`. That is the shape of an
   out-of-bounds write. Contiguity alone is not evidence — this submesh is numbered
   geographically, so consecutive indices are also spatial neighbours — but it is worth
   excluding.

Cheapest next step for either reading: the leg is bit-reproducible (`omp=1`), so a rerun with
a print of `UV(1,4:5,39028:39030)` after each of `compute_vel_rhs`, `viscosity_filter`,
`impl_vert_visc_ale` and `update_vel` pins the exact call that writes the zero. That needs a
coupled rerun of the 1906 leg rather than a standalone utility.

## 20. Superseded: the earlier "next step"

`vert_vel_ale` under `linfs` has no explicit zeroing: `Wvel(nz,n)` is the bottom-up cumulative
divergence divided by `area(nz,n)`, and an *exact* 0.0 cannot come out of that sum by chance.
What is needed next is `area(nz,n)` / `areasvol(nz,n)` at these five columns — not in the
blowup file, so it needs either a one-line diagnostic in the instrumented build (cheap: the
leg is bit-reproducible and crashes on step 1) or a reconstruction from the mesh. If
`area(nz,n)` is degenerate at exactly those levels, the same degeneracy explains both the
W hole and the singular tracer tridiagonal — one fault, both symptoms.

The §8 "compare the crash node against equally-pinched survivors" discriminator is no longer
the right test: the discriminator already exists and is exact (5 nodes out of 208 109, zero
false positives across six healthy states). The question is what *makes* those five columns
degenerate, not which pinched node is worst.

## Probes (`notes/ssh_rhs_probes/`)

- `probe_crash.py` — incoming restart fields at crash node (ssh/u/v/urhs_AB/w).
- `pgf.py` — T/S profiles + baroclinic PGF estimate across the crash element.
- `local.py` — mesh quality + new-node check in the crash box.
- `draftchange.py` — per-cycle draft change (via `map_nod.out` stable IDs).
- `pinch.py` / `prevalence.py` — per-level connectivity and trigger-event counts.
- `blowup_uvrhs.py` — opens the blowup dump; per-element `u/v/u_rhs/v_rhs/urhs_AB/helem`
  columns at the crash node, ssh_rhs vs ssh_rhs_old vs d_eta. (§9, §10)
- `eta_seam.py` — eta under cavity vs open ocean; the 1.59 m front step and its implied
  `-g∇η·dt`; native vs remapped restarts. (§10)
- `front_census.py` — `max|u_rhs|` and `|∫u_rhs h|` over all 1311 front elements against
  full-cavity and open controls; ssh_rhs population; pinched-node ranking. (§9, §10)
- `discriminator.py` — scans every column for interior zeros in `hnode` and `w`. (§11)
- `orphan_levels.py` — element-vs-node upper-level conventions; falsifies the orphan-level
  idea. (§11)
- `wdiv.py` — FVM reconstruction of the divergence column at the crash nodes (geometry is
  approximate — use for structure, not magnitude).
- `read_areas.py` — reads `mesh_areas.bin` from the instrumented `fesom_meshdiag`; scans every
  column for degenerate `area`/`areasvol`. (§12)
- `uv_holes.py` — interior exact zeros in `u`/`v`/`helem`, blowup vs incoming restart, using
  the model's own element level bounds. (§12)

Instrumented build (this session, kept):
`/work/ab0246/a270092/tmp/fesom_meshdiag_build/bin/fesom_meshdiag` — separate build dir, the
production `fesom-2.7/build/` and `bin/fesom.x` are untouched. Run dir
`/work/ab0246/a270092/tmp/meshdiag_run/` (namelists copied from the crash leg, mesh symlinked
including `dist_1792`, `run_meshdiag.sh` resubmittable). Two source edits in
`model_codes/awiesm3-develop-is/fesom-2.7/src/`: `fesom_meshdiag.F90` gained `dump_areas`
(called right after `mesh_setup`, before `ocean_setup` — the latter segfaults in the
climatology init and is not needed); `write_step_info.F90` has the `[FLUXGUARD]` block wrapped
in `#if defined(__oasis)` because `ice%atmcoupl` exists only in the coupled build and the tree
otherwise cannot compile any uncoupled utility. Neither changes the coupled build.

Key data: submeshes `couple/submesh_19{01..06}-12-31T00:00:00/` (with
`map_nod.out`, `cavity_nlvls.out`, `cavity_elvls.out`, `nlvls.out`), restart
`run_awiesm3_19060101-19061231/work/fesom.1905.oce.restart/`, and — the one that matters —
`run_awiesm3_19060101-19061231/work/fesom.1906.oce.blowup.nc` (2 GB, full 3D crash-step state).
Note the submesh whose node count matches a restart is the mesh that restart is written on;
match by `N`, not by the directory's year stamp.
