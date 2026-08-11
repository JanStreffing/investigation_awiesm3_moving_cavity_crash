#!/usr/bin/env python
"""Repair the flipped cells of an existing modified ICMGG, in place, into a _v2 copy.

Regenerating a modified ICMGG from scratch does not reproduce the land-sea mask of
files built months ago against an older mesh, and a changed coastline is a far larger
change than the defect being fixed. So instead of regenerating, take the published
file as it stands, find the cells that ocp-tool flipped when it fitted the mask, and
rebuild only those from their nearest stable neighbour with ocp-tool's own routine.
The mask, the field set and every other cell come through untouched.
"""
import sys
import numpy as np
import eccodes as ec

sys.path.insert(0, "/work/ab0246/a270092/software/ocp-tool")
from ocp_tool.lsm import fill_flipped_from_nearest_neighbour


class _Grid:
    """Minimal stand-in for GaussianGrid: the fill only reads the cell centres."""
    def __init__(self, lats, lons):
        self.center_lats = [lats]
        self.center_lons = [lons]


def read_fields(path):
    names, values, latlon = [], [], None
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            sn = ec.codes_get(gid, "shortName")
            names.append(sn)
            values.append(ec.codes_get_values(gid))
            if sn == "lsm" and latlon is None:
                latlon = (ec.codes_get_array(gid, "latitudes"),
                          ec.codes_get_array(gid, "longitudes"))
            ec.codes_release(gid)
    return names, values, latlon


def read_lsm(path):
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                return None
            if ec.codes_get(gid, "shortName") == "lsm":
                v = ec.codes_get_values(gid)
                ec.codes_release(gid)
                return v
            ec.codes_release(gid)


def write_fields(src_path, values, out_path):
    with open(src_path, "rb") as fin, open(out_path, "wb") as fout:
        i = 0
        while True:
            gid = ec.codes_grib_new_from_file(fin)
            if gid is None:
                break
            if len(values[i]) == ec.codes_get_size(gid, "values"):
                ec.codes_set_values(gid, values[i])
            ec.codes_write(gid, fout)
            ec.codes_release(gid)
            i += 1
    return i


def repair(prod_path, base_path, out_path):
    names, values, latlon = read_fields(prod_path)
    base_lsm = read_lsm(base_path)
    lsm_id = names.index("lsm")
    prod_lsm = values[lsm_id]

    to_land = np.where((base_lsm < 0.5) & (prod_lsm >= 0.5))[0].tolist()
    to_ocean = np.where((base_lsm >= 0.5) & (prod_lsm < 0.5))[0].tolist()
    print(f"  flipped to land {len(to_land)}, to ocean {len(to_ocean)}")

    grid = _Grid(np.asarray(latlon[0], dtype=float), np.asarray(latlon[1], dtype=float))
    fill_flipped_from_nearest_neighbour(values, lsm_id, grid, to_land, to_ocean,
                                        verbose=True)

    n = write_fields(prod_path, values, out_path)
    print(f"  wrote {n} messages to {out_path}")
    return to_land, to_ocean


if __name__ == "__main__":
    prod, base, out = sys.argv[1], sys.argv[2], sys.argv[3]
    print(f"repair {prod}")
    repair(prod, base, out)
