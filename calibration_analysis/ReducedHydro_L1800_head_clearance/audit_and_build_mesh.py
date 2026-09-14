"""Audit an imported rigid C3D4 mesh and build a frozen-physics input deck."""
from pathlib import Path
import argparse
import csv
import json
import math
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OLD = REPO / "calibration_analysis" / "ReducedHydro_geometry_L1800_validation"
FREECAD = REPO / "calibration_analysis" / "ReducedHydro_FreeCAD_L1800_validation"
sys.path.insert(0, str(FREECAD))
from build_and_audit_freecad import block_for_mesh, exterior, rotation_x_to_axis, tetra_properties

RP = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
R_BODY = 0.4075


def exact_normal(point, head_radius, head_extent):
    x, y, z = point
    radial = np.array([0.0, y, z]); rr = np.linalg.norm(radial)
    if x >= 1.8 - 1e-7:
        return np.array([1.0, 0.0, 0.0]), "tail"
    if x >= head_extent - 1e-8:
        return radial / max(rr, 1e-30), "body"
    radial_center = R_BODY - head_radius
    normal = np.array([x - head_extent,
                       (rr - radial_center) * y / max(rr, 1e-30),
                       (rr - radial_center) * z / max(rr, 1e-30)])
    return normal / np.linalg.norm(normal), "head"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix")
    parser.add_argument("geometry_json")
    parser.add_argument("job")
    args = parser.parse_args()
    raw = json.loads((HERE / f"{args.prefix}_raw.json").read_text())
    geometry = json.loads(Path(args.geometry_json).read_text())
    head_radius = geometry.get("R_head_design_mm", geometry.get("R_head_mm"))
    head_extent = geometry["head_axial_extent_mm"]
    props_name = Path(args.geometry_json).name.replace("_geometry.json", "_mass_properties.json")
    cad_mass = json.loads((Path(args.geometry_json).parent / props_name).read_text())
    cad_com = np.asarray(cad_mass["center_of_mass"]["value"], float)
    cad_inertia = np.sort(np.asarray(cad_mass["principal_mass_moments"]["value_tonne_mm2"], float))
    cad_mass_tonne = float(cad_mass["mass"]["tonne"])

    nodes = pd.read_csv(HERE / f"{args.prefix}_nodes_local.csv")
    elements_frame = pd.read_csv(HERE / f"{args.prefix}_elements.csv")
    labels = nodes.node.to_numpy(int)
    local = nodes[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    lookup = {label: i for i, label in enumerate(labels)}
    element_labels = elements_frame.element.to_numpy(int)
    elements = np.asarray([[lookup[int(v)] for v in row]
                           for row in elements_frame[["n1", "n2", "n3", "n4"]].to_numpy()], int)
    rotation = rotation_x_to_axis()
    xyz = RP + (local - cad_com) @ rotation.T
    volume, mass, com, inertia, signed = tetra_properties(xyz, elements)
    principal = np.sort(np.linalg.eigvalsh(inertia))
    faces = exterior(elements)

    normal_rows = []
    for element_index, side, face in faces:
        triangle = local[list(face)]
        center = triangle.mean(axis=0)
        normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
        normal /= np.linalg.norm(normal)
        exact, region = exact_normal(center, head_radius, head_extent)
        angle = math.degrees(math.acos(np.clip(abs(np.dot(normal, exact)), -1.0, 1.0)))
        normal_rows.append({"element": int(element_labels[element_index]), "side": side,
                            "region": region, "angle_error_deg": angle})
    normals = pd.DataFrame(normal_rows)
    normals.to_csv(HERE / f"{args.prefix}_surface_normal_audit.csv", index=False)
    region_stats = {region: {"triangle_count": int(len(group)),
                             "median_deg": float(group.angle_error_deg.median()),
                             "P95_deg": float(group.angle_error_deg.quantile(.95)),
                             "max_deg": float(group.angle_error_deg.max())}
                    for region, group in normals.groupby("region")}

    surface_index = np.unique(np.asarray([face for _, _, face in faces]).ravel())
    axial = local[surface_index, 0]
    radial = np.linalg.norm(local[surface_index, 1:], axis=1)
    length = float(axial.max() - axial.min()); diameter = float(2.0 * radial.max())
    mass_error = abs(mass - cad_mass_tonne) / cad_mass_tonne
    inertia_error = np.abs(principal - cad_inertia) / cad_inertia
    com_error = float(np.linalg.norm(com - RP))
    p95 = float(normals.angle_error_deg.quantile(.95))
    head_p95 = region_stats["head"]["P95_deg"]
    gates = {key: bool(value) for key, value in {
        "element_count_preferred": int(len(elements)) < 20000,
        "surface_normal_hard": p95 <= 2.5,
        "surface_normal_acceptable": p95 <= 2.0,
        "head_surface_normal_hard": head_p95 <= 2.5,
        "dimensions": abs(length - 1.8) <= .005 and abs(diameter - .815) <= .005,
        "mass": mass_error < .005,
        "COM": com_error < .005,
        "inertia": float(inertia_error.max()) < .02,
        "head_tail_polarity": bool(np.dot(rotation @ np.array([1.0, 0.0, 0.0]),
                                           np.array([.9647382600216, -.118874237214, .2348382536499])) > .999999),
    }.items()}
    # Element count is preferred rather than a hard geometry gate.
    required = ["surface_normal_hard", "head_surface_normal_hard", "dimensions",
                "mass", "COM", "inertia", "head_tail_polarity"]

    block = block_for_mesh(labels, xyz, element_labels, elements, faces)
    deck = (OLD / f"{OLD_JOB}.inp").read_text()
    start = deck.index("** ROBOT GEOMETRY CANDIDATE:")
    end = deck.index("*End Part", start) + len("*End Part")
    deck = deck[:start] + block + deck[end:]
    deck = deck.replace(OLD_JOB, args.job)
    import re
    deck = re.sub(r"(\*Instance, name=Robot_SOLID-1, part=Robot_SOLID\n)[^*]+?(\*End Instance)",
                  r"\1\2", deck, count=1)
    deck = re.sub(r"(\*Elset, elset=ROBOT_SOLID_CEL_ALL, instance=Robot_SOLID-1, generate\n)\d+,\s*\d+,\s*1",
                  rf"\g<1>{int(element_labels.min())}, {int(element_labels.max())}, 1", deck, count=1)
    inp = HERE / f"{args.job}.inp"; inp.write_text(deck)

    pd.DataFrame({"node": labels[surface_index], "x_mm": xyz[surface_index, 0],
                  "y_mm": xyz[surface_index, 1], "z_mm": xyz[surface_index, 2]}).to_csv(
                      HERE / f"{args.prefix}_surface_nodes_global.csv", index=False)
    with (HERE / f"{args.prefix}_surface_triangles.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["element", "side", "n1", "n2", "n3"])
        for element_index, side, face in faces:
            writer.writerow([int(element_labels[element_index]), side] + [int(labels[i]) for i in face])

    result = {
        "prefix": args.prefix, "job": args.job, "source_step": raw["source_step"],
        "mesh_size_mm": raw["mesh_size_mm"], "preprocessing_wall_clock_s": raw["preprocessing_wall_clock_s"],
        "node_count": len(labels), "element_count": len(elements),
        "exterior_node_count": len(surface_index), "exterior_triangle_count": len(faces),
        "length_mm": length, "diameter_mm": diameter, "volume_mm3": volume,
        "mass_mg": mass * 1e9, "mass_relative_error": mass_error,
        "mesh_COM_mm": com.tolist(), "COM_error_um": com_error * 1e3,
        "principal_inertia_tonne_mm2": principal.tolist(),
        "inertia_relative_error": inertia_error.tolist(),
        "surface_normal_median_deg": float(normals.angle_error_deg.median()),
        "surface_normal_P95_deg": p95, "surface_normal_max_deg": float(normals.angle_error_deg.max()),
        "surface_normal_by_region": region_stats, "input_size_bytes": inp.stat().st_size,
        "gates": gates, "all_required_geometry_gates_pass": all(gates[key] for key in required),
    }
    (HERE / f"{args.prefix}_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
