"""Extract dense RP histories and robot General Contact fields from one ODB."""
from __future__ import print_function

import json
import os
import sys

import numpy as np
from odbAccess import openOdb

job, destination = sys.argv[1:3]
os.makedirs(destination, exist_ok=True)
odb = openOdb(job + ".odb", readOnly=True)
try:
    step = list(odb.steps.values())[-1]
    rp, contact_history, energy, inventory = {}, {}, {}, []
    wanted = [prefix + str(axis) for prefix in ("U", "UR", "V", "VR", "A", "AR") for axis in (1, 2, 3)]
    for region_name, region in step.historyRegions.items():
        for variable, history in region.historyOutputs.items():
            array = np.asarray(history.data, dtype=float)
            inventory.append((region_name, variable, len(array)))
            if variable in wanted:
                current = rp.get(variable)
                if current is None or len(array) > len(current):
                    rp[variable] = array
            if variable.startswith(("CFN", "CFS", "CFT")):
                contact_history[region_name + "|" + variable] = array
            if variable.startswith(("ALL", "ETOTAL")):
                current = energy.get(variable)
                if current is None or len(array) > len(current):
                    energy[variable] = array
    np.savez_compressed(os.path.join(destination, "rp_history_private.npz"), **rp)
    np.savez_compressed(os.path.join(destination, "contact_history_private.npz"), **contact_history)
    np.savez_compressed(os.path.join(destination, "energy_history_private.npz"), **energy)

    inst = odb.rootAssembly.instances["ROBOT_SOLID-1"]
    labels = np.asarray([node.label for node in inst.nodes], dtype=np.int32)
    label_index = dict((int(label), i) for i, label in enumerate(labels))
    source_names = {}
    for prefix in ("CNORMF", "CSHEARF", "COPEN"):
        source_names[prefix] = next((key for key in step.frames[-1].fieldOutputs.keys()
                                     if key.strip().startswith(prefix) and "General_Contact_Domain" in key), None)
    field = {"time": [], "CNORMF": [], "CSHEARF": [], "COPEN": []}
    for frame in step.frames:
        field["time"].append(frame.frameValue)
        for prefix in ("CNORMF", "CSHEARF"):
            values = np.zeros((len(labels), 3), dtype=np.float32)
            name = source_names[prefix]
            if name and name in frame.fieldOutputs:
                for value in frame.fieldOutputs[name].getSubset(region=inst).values:
                    index = label_index.get(value.nodeLabel)
                    if index is not None:
                        values[index, :] = value.data[:3]
            field[prefix].append(values)
        values = np.full(len(labels), np.nan, dtype=np.float32)
        name = source_names["COPEN"]
        if name and name in frame.fieldOutputs:
            for value in frame.fieldOutputs[name].getSubset(region=inst).values:
                index = label_index.get(value.nodeLabel)
                if index is not None:
                    values[index] = float(value.data)
        field["COPEN"].append(values)
    np.savez_compressed(
        os.path.join(destination, "robot_contact_fields_private.npz"),
        time=np.asarray(field["time"]), node_labels=labels,
        node_coordinates_mm=np.asarray([node.coordinates for node in inst.nodes], dtype=np.float32),
        CNORMF=np.asarray(field["CNORMF"]), CSHEARF=np.asarray(field["CSHEARF"]),
        COPEN=np.asarray(field["COPEN"]),
        source_names=np.asarray([source_names["CNORMF"], source_names["CSHEARF"], source_names["COPEN"]]),
    )
    with open(os.path.join(destination, "odb_inventory.json"), "w") as handle:
        json.dump({"history": inventory, "frames": len(step.frames), "contact_field_names": source_names}, handle, indent=2)
    print("EXTRACTED", job, "frames", len(step.frames), "RP", sorted(rp), "contact fields", source_names)
finally:
    odb.close()
