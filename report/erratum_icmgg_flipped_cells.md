# Erratum: coastal land points in the AWI-ESM3 TCO319-DARS2 initial atmosphere

The OpenIFS initial file used by this configuration, `ICMGGabnsINIT_DARS2`, was
prepared with ocp-tool, which adjusts the OpenIFS land-sea mask so that it agrees
with the FESOM coastline. In the version used to build that file, every grid cell
the adjustment turned from ocean into land was assigned soil type 6, organic, and
kept the rest of its surface column at the ocean values it held before the flip.
Those cells therefore enter the simulation as land with organic soil and with no
water in any of the four soil layers. This affects 1047 cells, 0.89 percent of the
117804 land points. They are coastal and form a line one cell wide where the
OpenIFS and FESOM coastlines differed. Interior land is unaffected, and so is the
LPJ-GUESS soil type file `slt_TCO319_DARS2.nc`, only in the sense that it is
extracted from the same initial file and therefore carries the identical cells.

Soil type sets the soil hydraulic properties in HTESSEL, and organic soil holds
considerably more water than the mineral types it replaced, so a cell that begins
empty takes correspondingly longer to reach equilibrium water content. The affected
cells should be treated as unspun for soil moisture through the early part of the
simulation, and results that depend on coastal land surface fluxes or on the
near-coastal soil moisture field should be read with that in mind. Quantities
integrated over land are affected in proportion to the 0.89 percent cell count.
Two further fields were carried over from the ocean state in the same cells. Sea
ice fraction has no effect, because OpenIFS assigns open water and sea ice tile
fractions only where the point is not land, so these points receive no such tile.
Lake cover is carried over as well, and where FLake is active it is used directly
as the lake tile fraction on a land point.

The cause was corrected in ocp-tool commit `863c94a` of 3 July 2026, which rebuilds
the whole surface column of a flipped cell from the nearest stable neighbour of its
new type. Initial files regenerated with the corrected version, together with the
LPJ-GUESS soil type files derived from them, are published under a `_v2` suffix.
The files described here are kept under their existing names so that this
simulation remains reproducible.
