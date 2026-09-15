"""Abaqus-Python extraction of independent rigid-node versus RP-UR rotations."""
from odbAccess import openOdb
from pathlib import Path
import csv
import math
import sys

import numpy as np


def rodrigues(vector):
    angle = float(np.linalg.norm(vector))
    if angle < 1.0e-14:
        return np.eye(3)
    axis = vector / angle
    cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * np.dot(cross, cross)


def kabsch(reference, current):
    x = reference - reference.mean(axis=0)
    y = current - current.mean(axis=0)
    u, _, vt = np.linalg.svd(np.dot(x.T, y))
    rotation = np.dot(vt.T, u.T)
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1.0
        rotation = np.dot(vt.T, u.T)
    return rotation


def angle_error(a, b):
    delta = np.dot(a.T, b)
    return math.degrees(math.acos(np.clip((np.trace(delta) - 1.0) / 2.0, -1.0, 1.0)))


odb_path, output_path = sys.argv[1:3]
odb = openOdb(odb_path, readOnly=True)
step = list(odb.steps.values())[-1]
instance = odb.rootAssembly.instances["ROBOT_SOLID-1"]
initial = {node.label: np.asarray(node.coordinates, dtype=float) for node in instance.nodes}

ur = {}
for region in step.historyRegions.values():
    for name in ("UR1", "UR2", "UR3"):
        if name in region.historyOutputs and len(region.historyOutputs[name].data) > len(ur.get(name, [])):
            ur[name] = np.asarray(region.historyOutputs[name].data, dtype=float)

rows = []
for frame in step.frames:
    if "U" not in frame.fieldOutputs:
        continue
    values = [value for value in frame.fieldOutputs["U"].values
              if value.instance is not None and value.instance.name == "ROBOT_SOLID-1"]
    if len(values) < 4:
        continue
    labels = [value.nodeLabel for value in values]
    x0 = np.asarray([initial[label] for label in labels])
    x1 = np.asarray([initial[label] + np.asarray(value.data, dtype=float) for label, value in zip(labels, values)])
    reference = kabsch(x0, x1)
    time_s = float(frame.frameValue)
    vector = np.array([np.interp(time_s, ur[name][:, 0], ur[name][:, 1]) for name in ("UR1", "UR2", "UR3")])
    server = rodrigues(vector)
    fitted = np.dot((x0 - x0.mean(axis=0)), reference.T) + x1.mean(axis=0)
    rms = float(np.sqrt(np.mean(np.sum((fitted - x1) ** 2, axis=1))))
    rows.append([time_s, *vector, math.degrees(np.linalg.norm(vector)), len(labels), rms, angle_error(server, reference)])
odb.close()

with open(output_path, "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["time_s", "UR1_rad", "UR2_rad", "UR3_rad", "UR_norm_deg", "rigid_node_count", "kabsch_rms_mm", "R_server_vs_R_nodes_error_deg"])
    writer.writerows(rows)
print("ORIENTATION_AUDIT", len(rows), "frames", output_path)
