"""Build the L2300 deck from the validated L1800 physics template and frozen mesh."""

from pathlib import Path
import csv
import hashlib
import json
import re

import numpy as np
import pandas as pd

from audit_solver_mesh import exterior, tetra_properties
from audit_frozen_cad_static import A0, RP, rotation_x_to_axis


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
L1800 = REPO / "calibration_analysis/ReducedHydro_FreeCAD_L1800_validation"
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083"
JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083"
STEP = REPO / "cad/freecad_parametric_robot/variants/L2300_D0815_wallwobble/Robot_L2300_D0815_WallWobble.step"
EXPECTED_SHA = "461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2"
DENSITY = 7.80906654321e-9
CAD_COM = np.array([1.2060186937156343, 0.0, 0.0])
MOMENT = 0.0010400401426980412


def mesh_block(labels, xyz, element_labels, elements, faces):
    side_sets = {i: [] for i in range(1, 5)}
    for element_index, side, _ in faces:
        side_sets[side].append(int(element_labels[element_index]))
    lines = [
        "** ROBOT SOURCE: frozen authoritative FreeCAD L2300 STEP; no coordinate scaling",
        "*Part, name=Robot_SOLID", "*Node",
    ]
    lines += [f"{int(label)}, {p[0]:.12g}, {p[1]:.12g}, {p[2]:.12g}" for label, p in zip(labels, xyz)]
    lines.append("*Element, type=C3D4")
    lines += [f"{int(label)}, " + ", ".join(str(int(labels[i])) for i in tet)
              for label, tet in zip(element_labels, elements)]
    for side in range(1, 5):
        lines.append(f"*Elset, elset=ROBOT_EXTERIOR_S{side}")
        values = side_sets[side]
        lines += [", ".join(str(value) for value in values[i:i + 16]) for i in range(0, len(values), 16)]
    lines += [
        "*Nset, nset=ROBOT_CEL_BODY, generate", f"{labels.min()}, {labels.max()}, 1",
        "*Elset, elset=ROBOT_CEL_BODY, generate", f"{element_labels.min()}, {element_labels.max()}, 1",
        "*Elset, elset=ROBOT_SOLID_ALL, generate", f"{element_labels.min()}, {element_labels.max()}, 1",
        "*Surface, type=ELEMENT, name=ROBOT_SOLID_SURF",
        "ROBOT_EXTERIOR_S1, S1", "ROBOT_EXTERIOR_S2, S2",
        "ROBOT_EXTERIOR_S3, S3", "ROBOT_EXTERIOR_S4, S4",
        "** Section: SEC_ROBOT_RIGID",
        "*Solid Section, elset=ROBOT_CEL_BODY, material=MAT_ROBOT_RIGID", ",", "*End Part",
    ]
    return "\n".join(lines)


def main():
    if hashlib.sha256(STEP.read_bytes()).hexdigest() != EXPECTED_SHA:
        raise RuntimeError("Frozen STEP SHA256 mismatch")
    acceptance = json.loads((HERE / "L2300_threshold_aware_acceptance.json").read_text())
    if not acceptance["accepted"]:
        raise RuntimeError("Threshold-aware mesh acceptance is not passed")
    nodes = pd.read_csv(HERE / "L2300_solver_mesh_nodes_local.csv")
    elems = pd.read_csv(HERE / "L2300_solver_mesh_elements.csv")
    labels = nodes.node.to_numpy(int)
    local = nodes[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    lookup = {label: i for i, label in enumerate(labels)}
    element_labels = elems.element.to_numpy(int)
    elements = np.asarray([[lookup[int(v)] for v in row]
                           for row in elems[["n1", "n2", "n3", "n4"]].to_numpy()], int)
    if len(labels) != 5828 or len(elements) != 29141:
        raise RuntimeError("Frozen mesh identity mismatch")
    rotation = rotation_x_to_axis(A0)
    xyz = RP + (local - CAD_COM) @ rotation.T
    volume, mass, com, inertia, signed = tetra_properties(xyz, elements)
    faces = exterior(elements)
    axis_dot = float(np.dot(rotation @ np.array([1.0, 0.0, 0.0]), A0))
    if axis_dot <= .999999 or np.linalg.norm(com - RP) >= .005 or np.any(signed == 0):
        raise RuntimeError("Placement or mesh-volume gate failed")

    template = (L1800 / f"{OLD_JOB}.inp").read_text()
    start = template.index("** ROBOT SOURCE:")
    end = template.index("*End Part", start) + len("*End Part")
    deck = template[:start] + mesh_block(labels, xyz, element_labels, elements, faces) + template[end:]
    deck = deck.replace(OLD_JOB, JOB)
    deck = re.sub(r"(\*Instance, name=Robot_SOLID-1, part=Robot_SOLID\n)[^*]+?(\*End Instance)",
                  r"\1\2", deck, count=1)
    deck = re.sub(r"(\*Elset, elset=ROBOT_SOLID_CEL_ALL, instance=Robot_SOLID-1, generate\n)\d+,\s*\d+,\s*1",
                  rf"\g<1>{element_labels.min()}, {element_labels.max()}, 1", deck, count=1)
    header = ("** FROZEN L2300 WALL-SUPPORTED SCREENING: 29141 C3D4; threshold-aware 5 um acceptance\n"
              "** CAD L=2.300000 mm D=0.815000 mm; no mesh scaling or refinement\n")
    deck = header + deck
    local_inp = HERE / f"{JOB}.inp"
    root_inp = ROOT / f"{JOB}.inp"
    local_inp.write_text(deck)
    root_inp.write_text(deck)

    surface_indices = np.unique(np.asarray([face for _, _, face in faces]).ravel())
    pd.DataFrame({"node": labels[surface_indices], "x_mm": xyz[surface_indices, 0],
                  "y_mm": xyz[surface_indices, 1], "z_mm": xyz[surface_indices, 2]}).to_csv(
                      HERE / "L2300_robot_surface_nodes_global.csv", index=False)
    with (HERE / "L2300_robot_surface_triangles.csv").open("w", newline="") as handle:
        writer = csv.writer(handle); writer.writerow(["element", "side", "n1", "n2", "n3"])
        for element_index, side, face in faces:
            writer.writerow([int(element_labels[element_index]), side] + [int(labels[i]) for i in face])
    identity = {
        "job": JOB, "source_STEP": str(STEP.relative_to(REPO)).replace("\\", "/"),
        "source_STEP_SHA256": EXPECTED_SHA, "scaleFromFile": False,
        "mesh_refinement_performed": False, "nodes": len(labels), "elements": len(elements),
        "element_type": "C3D4", "exterior_triangles": len(faces),
        "RP_mm": RP.tolist(), "mesh_COM_mm": com.tolist(),
        "RP_mesh_COM_error_um": float(np.linalg.norm(com - RP) * 1e3),
        "head_tail_dot_a0": axis_dot, "mesh_volume_mm3": float(volume),
        "mesh_mass_mg": float(mass * 1e9),
        "mesh_principal_inertia_tonne_mm2": np.sort(np.linalg.eigvalsh(inertia)).tolist(),
        "magnetic_moment_Am2": MOMENT,
        "duplicate_mass_keyword_count": len(re.findall(r"(?im)^\*mass\b", deck)),
        "general_contact_count": len(re.findall(r"(?im)^\*contact\s*$", deck)),
        "step_time_s": .008333, "direct_dt_s": 1e-7,
        "field_interval_s": 2.5e-5, "history_frequency": 1,
        "threshold_aware_acceptance": acceptance["classification"],
        "ready_for_datacheck": True,
    }
    if identity["duplicate_mass_keyword_count"] != 0 or identity["general_contact_count"] != 1:
        raise RuntimeError("Deck mass/contact identity failed")
    (HERE / "L2300_deck_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
