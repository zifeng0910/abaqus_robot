from __future__ import print_function

import csv
import os
import sys

from odbAccess import openOdb


PROBES = ["PROBE_CENTER", "PROBE_LOW_QUARTER", "PROBE_HIGH_QUARTER", "PROBE_LOW_END_NEAR", "PROBE_HIGH_END_NEAR"]


def mean_scalar(values):
    vals = []
    for value in values:
        data = value.data
        if hasattr(data, "__len__"):
            vals.append(float(data[0]))
        else:
            vals.append(float(data))
    return sum(vals) / len(vals) if vals else float("nan")


def mean_pressure(values):
    vals = []
    for value in values:
        data = value.data
        if len(data) >= 3:
            vals.append(-(float(data[0]) + float(data[1]) + float(data[2])) / 3.0)
    return sum(vals) / len(vals) if vals else float("nan")


def extract(odb_path, csv_path):
    odb = openOdb(odb_path, readOnly=True)
    step = odb.steps[list(odb.steps.keys())[0]]
    inst = odb.rootAssembly.instances["FLUID_EULERIAN-1"]
    rows = []
    for frame in step.frames:
        t = float(frame.frameValue)
        fo = frame.fieldOutputs
        evf_name = next((key for key in fo.keys() if key.startswith("EVF_") and not key.endswith("VOID")), None)
        stress_name = next((key for key in fo.keys() if key.startswith("S_") and "ASSEMBLY_FLUID_EULERIAN-1_WATER" in key), None)
        evf = fo[evf_name] if evf_name else None
        stress = fo[stress_name] if stress_name else None
        vel = fo["V"] if "V" in fo else None
        row = {"time_s": t}
        for probe in PROBES:
            eset = inst.elementSets[probe]
            nset = inst.nodeSets[probe + "_NODES"]
            row[probe + "_p_proxy"] = mean_pressure(stress.getSubset(region=eset).values) if stress else float("nan")
            row[probe + "_evf"] = mean_scalar(evf.getSubset(region=eset).values) if evf else float("nan")
            row[probe + "_v_axial"] = mean_scalar(vel.getSubset(region=nset).values) if vel else float("nan")
        rows.append(row)
    odb.close()
    fields = ["time_s"] + [p + suffix for p in PROBES for suffix in ("_p_proxy", "_evf", "_v_axial")]
    with open(csv_path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: abaqus python extract_f100_cel_boundary_harness.py input.odb output.csv")
    extract(sys.argv[1], sys.argv[2])
