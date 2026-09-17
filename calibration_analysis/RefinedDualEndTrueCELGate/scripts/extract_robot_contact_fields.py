from __future__ import print_function

import os
import sys

import numpy as np
from odbAccess import openOdb


job, outdir = sys.argv[1:]
os.makedirs(outdir, exist_ok=True)
odb = openOdb(job + ".odb", readOnly=True)
try:
    step = list(odb.steps.values())[-1]
    inst = odb.rootAssembly.instances["ROBOT_SOLID-1"]
    labels = np.asarray([node.label for node in inst.nodes], dtype=np.int32)
    index = dict((int(label), i) for i, label in enumerate(labels))
    names = {}
    for prefix in ("CNORMF", "CSHEARF", "COPEN"):
        names[prefix] = next((key for key in step.frames[-1].fieldOutputs.keys()
                              if key.strip().startswith(prefix) and "General_Contact_Domain" in key), None)
    data = {"time": [], "node_labels": labels,
            "node_coordinates_mm": np.asarray([node.coordinates for node in inst.nodes], dtype=np.float32)}
    for prefix in ("CNORMF", "CSHEARF"):
        data[prefix] = []
    data["COPEN"] = []
    for frame in step.frames:
        data["time"].append(frame.frameValue)
        for prefix in ("CNORMF", "CSHEARF"):
            values = np.zeros((len(labels), 3), dtype=np.float32)
            if names[prefix] in frame.fieldOutputs:
                for value in frame.fieldOutputs[names[prefix]].getSubset(region=inst).values:
                    if value.nodeLabel in index:
                        values[index[value.nodeLabel], :] = value.data[:3]
            data[prefix].append(values)
        values = np.full(len(labels), np.nan, dtype=np.float32)
        if names["COPEN"] in frame.fieldOutputs:
            for value in frame.fieldOutputs[names["COPEN"]].getSubset(region=inst).values:
                if value.nodeLabel in index:
                    values[index[value.nodeLabel]] = float(value.data)
        data["COPEN"].append(values)
    np.savez_compressed(os.path.join(outdir, "robot_contact_fields_private.npz"),
                        time=np.asarray(data["time"]), node_labels=labels,
                        node_coordinates_mm=data["node_coordinates_mm"],
                        CNORMF=np.asarray(data["CNORMF"]), CSHEARF=np.asarray(data["CSHEARF"]),
                        COPEN=np.asarray(data["COPEN"]),
                        source_names=np.asarray([names["CNORMF"], names["CSHEARF"], names["COPEN"]]))
    print(names)
finally:
    odb.close()
