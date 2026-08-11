#!/usr/bin/env python
"""Build the _v2 ICMGG and LPJ-GUESS slt set from the published files.

For each published modified ICMGG: repair the cells ocp-tool flipped when it fitted
the land-sea mask, add the CO2 emission fields where the file predates them so the
carbon-cycle setup has what it needs, then extract the soil type file the same way
ocp-tool does. The land-sea mask is never touched, so each _v2 file is a drop-in for
the one beside it.
"""
import os
import subprocess
import sys

import eccodes as ec

sys.path.insert(0, "/work/ab0246/a270092/software/ocp-tool")
sys.path.insert(0, "/work/ab0246/a270092/tmp")
from repair_icmgg_v2 import repair
from ocp_tool.field_interpolation import interpolate_2d_fields_to_icmgg

OCP = "/work/ab0246/a270092/software/ocp-tool"
CO2_EMIS = f"{OCP}/input/openifs_input_default/cams_co2_emissions.grib"
OUT = "/work/ab0246/a270092/tmp/v2out"

# (pool dir, base ICMGG, modified ICMGG, slt grid tag, ocean tag)
CASES = [
    ("oifs-48r1/TCO95L91",   "ICMGGab45INIT", "ICMGGab45INIT_CORE3",     "TCO95",  "CORE3"),
    ("oifs-48r1/TCO95L91",   "ICMGGab45INIT", "ICMGGab45INIT_CORE2",     "TCO95",  "CORE2"),
    ("oifs-48r1/TCO95L91",   "ICMGGab45INIT", "ICMGGab45INIT_CORE2ice",  "TCO95",  "CORE2ice"),
    ("oifs-48r1/TCO95L91",   "ICMGGab45INIT", "ICMGGab45INIT_DARS2",     "TCO95",  "DARS2"),
    ("oifs-48r1/TCO319L137", "ICMGGabnsINIT", "ICMGGabnsINIT_DARS2",     "TCO319", "DARS2"),
    ("oifs-48r1/TCO319L137", "ICMGGabnsINIT", "ICMGGabnsINIT_DARS2cav",  "TCO319", "DARS2cav"),
    ("oifs-48r1/TCO319L137", "ICMGGabnsINIT", "ICMGGabnsINIT_CORE3",     "TCO319", "CORE3"),
    ("oifs-48r1/TCO319L137", "ICMGGabnsINIT", "ICMGGabnsINIT_CORE2",     "TCO319", "CORE2"),
    ("oifs-48r1/TCO319L137", "ICMGGabnsINIT", "ICMGGabnsINIT_DART",      "TCO319", "DART"),
    ("oifs-48r1/TL255L91",   "ICMGGabl7INIT", "ICMGGabl7INIT_CORE2",     "TL255",  "CORE2"),
]
POOL = "/work/ab0246/a270092/input"


def count_co2(path):
    n = 0
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            if ec.codes_get(gid, "shortName").lower().startswith("co2"):
                n += 1
            ec.codes_release(gid)
    return n


def make_slt(icmgg, out_nc):
    tmp = out_nc + ".grb"
    for cmd in (f"grib_copy -w shortName=slt {icmgg} {tmp}",
                f"cdo -s -f nc copy {tmp} {out_nc}",
                f"ncrename -v slt,var43 {out_nc}"):
        rc = subprocess.call(cmd, shell=True)
        if rc != 0 and "ncrename" not in cmd:
            raise RuntimeError(f"failed: {cmd}")
    os.path.exists(tmp) and os.remove(tmp)


def main():
    os.makedirs(OUT, exist_ok=True)
    for subdir, base, mod, grid, ocean in CASES:
        src = f"{POOL}/{subdir}/{mod}"
        basef = f"{POOL}/{subdir}/{base}"
        if not (os.path.exists(src) and os.path.exists(basef)):
            print(f"SKIP {mod}: missing input")
            continue
        dst = f"{OUT}/{mod}_v2"
        print(f"\n=== {mod} ===")
        repair(src, basef, dst)

        if count_co2(src) == 0:
            print("  no CO2 emission fields in the source, adding them for -cc")
            interpolate_2d_fields_to_icmgg(CO2_EMIS, dst, output_file=dst,
                                           variable_name="lsm",
                                           field_type="co2_emissions", verbose=False)
            print(f"  CO2 fields now: {count_co2(dst)}")
        else:
            print(f"  CO2 emission fields already present ({count_co2(src)}), left alone")

        slt = f"{OUT}/slt_{grid}_{ocean}_v2.nc"
        make_slt(dst, slt)
        print(f"  slt -> {slt}")


if __name__ == "__main__":
    main()
