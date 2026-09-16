"""Extract true-CEL field frames from an Abaqus ODB into a private NPZ."""
from __future__ import print_function

import json
import os
import sys

import numpy as np
from odbAccess import openOdb


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: abaqus python extract_truecel_odb.py JOB OUTDIR")
    job, outdir = sys.argv[1:]
    os.makedirs(outdir, exist_ok=True)
    odb = openOdb(job + ".odb", readOnly=True)
    try:
        step = list(odb.steps.values())[-1]
        instance = odb.rootAssembly.instances["FLUID_EULERIAN-1"]
        node_labels = np.asarray([node.label for node in instance.nodes], dtype=np.int32)
        node_coords = np.asarray([node.coordinates for node in instance.nodes], dtype=float)
        element_labels = np.asarray([element.label for element in instance.elements], dtype=np.int32)
        node_index = dict((int(label), index) for index, label in enumerate(node_labels))
        element_index = dict((int(label), index) for index, label in enumerate(element_labels))
        times, velocity_frames, evf_frames = [], [], []
        field_keys = sorted(set(key for frame in step.frames for key in frame.fieldOutputs.keys()))
        evf_key = next((key for key in field_keys if key == "EVF" or key.startswith("EVF_")), None)
        if "V" not in field_keys or evf_key is None:
            raise RuntimeError("required CEL fields V/EVF missing: {}".format(field_keys))
        for frame in step.frames:
            velocity = np.full((len(node_labels), 3), np.nan, dtype=np.float32)
            for value in frame.fieldOutputs["V"].getSubset(region=instance).values:
                if value.nodeLabel in node_index:
                    velocity[node_index[value.nodeLabel], :] = value.data[:3]
            evf = np.full(len(element_labels), np.nan, dtype=np.float32)
            for value in frame.fieldOutputs[evf_key].getSubset(region=instance).values:
                if value.elementLabel in element_index:
                    data = value.data
                    try:
                        scalar = float(data)
                    except TypeError:
                        scalar = float(np.asarray(data).reshape(-1)[0])
                    evf[element_index[value.elementLabel]] = scalar
            times.append(frame.frameValue)
            velocity_frames.append(velocity)
            evf_frames.append(evf)
        np.savez_compressed(
            os.path.join(outdir, "truecel_field_private.npz"),
            time=np.asarray(times, dtype=float),
            fluid_node_labels=node_labels,
            fluid_node_coordinates_mm=node_coords,
            fluid_velocity_mm_s=np.asarray(velocity_frames, dtype=np.float32),
            fluid_element_labels=element_labels,
            fluid_evf=np.asarray(evf_frames, dtype=np.float32),
        )
        inventory = {
            "job": job, "step": step.name, "frame_count": len(times),
            "fluid_node_count": len(node_labels), "fluid_element_count": len(element_labels),
            "field_keys": field_keys, "evf_field": evf_key, "has_V": True, "has_EVF": True,
        }
        with open(os.path.join(outdir, "truecel_odb_inventory.json"), "w") as handle:
            json.dump(inventory, handle, indent=2)
            handle.write("\n")
        print(json.dumps(inventory, indent=2))
    finally:
        odb.close()


if __name__ == "__main__":
    main()
