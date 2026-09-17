"""Abaqus Python: export dense RP U/UR histories to one CSV."""
from __future__ import print_function

import csv
import sys
from odbAccess import openOdb


job = sys.argv[1]
odb = openOdb(job + ".odb", readOnly=True)
step = list(odb.steps.values())[-1]
wanted = ["U1", "U2", "U3", "UR1", "UR2", "UR3"]
series = {}
for region in step.historyRegions.values():
    for name in wanted:
        if name in region.historyOutputs:
            candidate = list(region.historyOutputs[name].data)
            if name not in series or len(candidate) > len(series[name]):
                series[name] = candidate
if sorted(series) != sorted(wanted):
    raise RuntimeError("missing RP histories: {}".format(sorted(set(wanted)-set(series))))
times = [row[0] for row in series["U1"]]
with open(job + "_rp_history.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["time_s"] + wanted)
    for index, time_value in enumerate(times):
        writer.writerow([time_value] + [series[name][index][1] for name in wanted])
odb.close()
print("EXTRACTED", job, len(times))
