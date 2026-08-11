#!/usr/bin/env python
"""Give every _v2 ICMGG the same complete set of CO2 emission fields.

The published files carry anything from zero to four of them depending on which
ocp-tool wrote them, and the carbon-cycle setup wants the full set. Strip whatever
CO2 messages a file has and add the four from cams_co2_emissions.grib, so the answer
does not depend on the file's history. Nothing else in the file is touched.
"""
import os
import shutil
import sys

import eccodes as ec

sys.path.insert(0, "/work/ab0246/a270092/software/ocp-tool")
from ocp_tool.field_interpolation import interpolate_2d_fields_to_icmgg

CO2_EMIS = "/work/ab0246/a270092/software/ocp-tool/input/openifs_input_default/cams_co2_emissions.grib"
V = "/work/ab0246/a270092/tmp/v2out"


def co2_names(path):
    out = []
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None:
                break
            sn = ec.codes_get(gid, "shortName")
            if sn.lower().startswith(("co2", "fco2")):
                out.append(sn)
            ec.codes_release(gid)
    return out


def strip_co2(path):
    tmp = path + ".stripped"
    kept = dropped = 0
    with open(path, "rb") as fin, open(tmp, "wb") as fout:
        while True:
            gid = ec.codes_grib_new_from_file(fin)
            if gid is None:
                break
            sn = ec.codes_get(gid, "shortName")
            if sn.lower().startswith(("co2", "fco2")):
                dropped += 1
            else:
                ec.codes_write(gid, fout)
                kept += 1
            ec.codes_release(gid)
    os.replace(tmp, path)
    return kept, dropped


def main():
    files = sorted(f for f in os.listdir(V) if f.startswith("ICMGG") and f.endswith("_v2"))
    for name in files:
        p = os.path.join(V, name)
        before = co2_names(p)
        kept, dropped = strip_co2(p)
        interpolate_2d_fields_to_icmgg(CO2_EMIS, p, output_file=p,
                                       variable_name="lsm",
                                       field_type="co2_emissions", verbose=False)
        after = co2_names(p)
        print(f"{name:28s} co2 {len(before)} -> {len(after)}  ({kept} other messages kept)  {sorted(after)}")


if __name__ == "__main__":
    main()
