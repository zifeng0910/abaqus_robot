"""Extract only robot-surface COPEN for the wall-only General Contact case."""
from __future__ import print_function

import sys
from pathlib import Path

import numpy as np
from odbAccess import openOdb

job, destination = sys.argv[1:3]
out = Path(destination)
out.mkdir(parents=True, exist_ok=True)
odb = openOdb(job + ".odb", readOnly=True)
try:
    step = list(odb.steps.values())[-1]
    inst = odb.rootAssembly.instances["ROBOT_SOLID-1"]
    nodes = list(inst.nodes)
    labels = np.asarray([n.label for n in nodes], dtype=np.int32)
    index = {int(label): i for i, label in enumerate(labels)}
    key = next(k for k in step.frames[-1].fieldOutputs.keys()
               if k.startswith("COPEN") and "ROBOT_SOLID-1_ROBOT_SOLID_SURF" in k)
    times = []
    opening = []
    for frame in step.frames:
        times.append(frame.frameValue)
        values = np.full(len(labels), np.nan, dtype=np.float32)
        if key in frame.fieldOutputs:
            for value in frame.fieldOutputs[key].getSubset(region=inst).values:
                if value.nodeLabel in index:
                    values[index[value.nodeLabel]] = float(value.data)
        opening.append(values)
    np.savez_compressed(out / "robot_wall_copen_private.npz",
                        time=np.asarray(times), node_labels=labels,
                        node_coordinates_mm=np.asarray([n.coordinates for n in nodes], dtype=np.float32),
                        COPEN=np.asarray(opening), source_name=key)
    print("EXTRACTED", len(times), "frames", len(labels), "robot nodes", key)
finally:
    odb.close()
