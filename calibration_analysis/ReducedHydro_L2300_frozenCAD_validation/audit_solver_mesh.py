"""Audit the global-0.060 rigid mesh against frozen CAD and exact-wall poses."""

from pathlib import Path
import json
import math
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CAD = REPO / "cad" / "freecad_parametric_robot" / "variants" / "L2300_D0815_wallwobble"
sys.path.insert(0, str(HERE))
from audit_frozen_cad_static import (
    A0, AZIMUTHS, RP, TILTS, Wall, certified_near_wall_query,
    local_frame, rotation_x_to_axis, sector_metrics,
)

DENSITY = 7.80906654321e-9
L_TOTAL = 2.300
R_BODY = .4075
HEAD_EXTENT = .26149477


def tetra_properties(nodes, elements):
    vertices = nodes[elements]
    signed = np.linalg.det(vertices[:, 1:]-vertices[:, :1])/6.0
    volume_each = np.abs(signed)
    centers = vertices.mean(axis=1)
    volume = volume_each.sum()
    mass = volume*DENSITY
    com = np.average(centers, axis=0, weights=volume_each)
    second = np.zeros((3, 3))
    for points, element_volume in zip(vertices, volume_each):
        total = points.sum(axis=0)
        second += element_volume/20.0*(np.outer(total, total)+points.T@points)
    second_com = DENSITY*second-mass*np.outer(com, com)
    inertia = np.eye(3)*np.trace(second_com)-second_com
    return volume, mass, com, inertia, signed


def exterior(elements):
    definitions = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
    found = {}
    for element_index, tet in enumerate(elements):
        for side, local_ids in enumerate(definitions, 1):
            face = tuple(int(tet[i]) for i in local_ids)
            key = tuple(sorted(face))
            found[key] = None if key in found else (element_index, side, face)
    return [value for value in found.values() if value is not None]


def triangle_samples(triangles, subdivisions=4):
    weights = np.asarray([(i/subdivisions, j/subdivisions, 1.0-(i+j)/subdivisions)
                          for i in range(subdivisions+1)
                          for j in range(subdivisions+1-i)], float)
    return np.einsum("wk,tkd->twd", weights, triangles).reshape(-1, 3), len(weights)


def sample_regions(triangles, samples_per_triangle):
    centroid_x = triangles[:, :, 0].mean(axis=1)
    face_region = np.where(centroid_x >= L_TOTAL-1e-7, "TAIL",
                           np.where(centroid_x <= HEAD_EXTENT+1e-8, "HEAD", "BODY"))
    return np.repeat(face_region, samples_per_triangle)


def exact_normals(points, profile):
    x = points[:, 0]
    radius = np.linalg.norm(points[:, 1:], axis=1)
    normal = np.zeros_like(points)
    region = np.where(x >= L_TOTAL-1e-7, "TAIL",
                      np.where(x <= HEAD_EXTENT+1e-8, "HEAD", "BODY"))
    tail = region == "TAIL"
    body = region == "BODY"
    head = region == "HEAD"
    normal[tail, 0] = 1.0
    normal[body, 1:] = points[body, 1:]/np.maximum(radius[body, None], 1e-30)
    hp = profile[profile.region == "HEAD"]
    xp = hp.x_mm.to_numpy(float)
    rp = hp.radius_mm.to_numpy(float)
    dx = np.gradient(xp)
    dr = np.gradient(rp)
    indices = np.clip(np.searchsorted(xp, x[head]), 1, len(xp)-1)
    choose_left = abs(x[head]-xp[indices-1]) < abs(x[head]-xp[indices])
    indices = indices-choose_left.astype(int)
    radial_unit = points[head, 1:]/np.maximum(radius[head, None], 1e-30)
    normal[head, 0] = -dr[indices]
    normal[head, 1:] = dx[indices, None]*radial_unit
    normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1e-30)
    return normal, region


def mesh_pose(samples, region, wall, frame, cad_com, tilt, azimuth):
    radial = math.cos(math.radians(azimuth))*frame["e1"]+math.sin(math.radians(azimuth))*frame["e2"]
    axis = math.cos(math.radians(tilt))*frame["forward_tangent"]+math.sin(math.radians(tilt))*radial
    rotation = rotation_x_to_axis(axis)
    points = RP+(samples-cad_com)@rotation.T
    gaps, wall_indices, queried = certified_near_wall_query(wall, points, region)
    row = {
        "tilt_deg": tilt, "azimuth_deg": azimuth,
        "mesh_min_gap_um": float(gaps.min()*1e3),
        "mesh_HEAD_min_gap_um": float(gaps[region == "HEAD"].min()*1e3),
        "mesh_BODY_min_gap_um": float(gaps[region == "BODY"].min()*1e3),
        "mesh_TAIL_min_gap_um": float(gaps[region == "TAIL"].min()*1e3),
        "mesh_contact": bool(gaps.min() <= 0.0),
        "exactly_queried_mesh_sample_count": queried,
    }
    clusters, separation, opposing = sector_metrics(
        gaps, wall_indices, wall, frame["e1"], frame["e2"], 20.0)
    row["mesh_sector_clusters_20um"] = clusters
    row["mesh_sector_separation_20um_deg"] = separation
    row["mesh_opposing_bridge_20um"] = opposing
    end_support = min(row["mesh_HEAD_min_gap_um"], row["mesh_TAIL_min_gap_um"]) <= 20.0
    row["mesh_WALL_SUPPORT"] = bool(end_support and not opposing)
    return row


def main():
    raw = json.loads((HERE/"L2300_solver_mesh_raw.json").read_text())
    nodes_df = pd.read_csv(HERE/"L2300_solver_mesh_nodes_local.csv")
    elements_df = pd.read_csv(HERE/"L2300_solver_mesh_elements.csv")
    labels = nodes_df.node.to_numpy(int)
    local = nodes_df[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    lookup = {label: i for i, label in enumerate(labels)}
    element_labels = elements_df.element.to_numpy(int)
    elements = np.asarray([[lookup[int(node)] for node in row]
        for row in elements_df[["n1", "n2", "n3", "n4"]].to_numpy()], int)
    faces = exterior(elements)
    face_nodes = np.asarray([face for _, _, face in faces], int)
    triangles = local[face_nodes]
    samples, samples_per_triangle = triangle_samples(triangles)
    sample_region = sample_regions(triangles, samples_per_triangle)

    props = json.loads((CAD/"Robot_L2300_D0815_WallWobble_mass_properties.json").read_text())
    cad_com = np.asarray(props["center_of_mass"]["value"], float)
    rotation = rotation_x_to_axis(A0)
    global_nodes = RP+(local-cad_com)@rotation.T
    volume, mass, com, inertia, signed = tetra_properties(global_nodes, elements)
    principal = np.sort(np.linalg.eigvalsh(inertia))
    target = np.sort(np.asarray(props["principal_mass_moments"]["value_tonne_mm2"], float))
    mass_target = props["mass"]["tonne"]
    mass_error = abs(mass-mass_target)/mass_target
    com_error = np.linalg.norm(com-RP)
    inertia_error = np.abs(principal-target)/target

    profile = pd.read_csv(CAD/"Robot_L2300_D0815_WallWobble_profile.csv")
    exact, face_region = exact_normals(triangles.mean(axis=1), profile)
    mesh_normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    mesh_normals /= np.maximum(np.linalg.norm(mesh_normals, axis=1)[:, None], 1e-30)
    angles = np.degrees(np.arccos(np.clip(np.abs(np.sum(mesh_normals*exact, axis=1)), -1.0, 1.0)))
    normal_table = pd.DataFrame({
        "element": [element_labels[index] for index, _, _ in faces],
        "side": [side for _, side, _ in faces], "region": face_region,
        "angle_error_deg": angles,
    })
    normal_table.to_csv(HERE/"L2300_solver_surface_normal_diagnostic.csv", index=False)

    wall = Wall()
    initial_points = RP+(samples-cad_com)@rotation.T
    initial_gaps, _, _ = certified_near_wall_query(wall, initial_points, sample_region)
    initial_gap_um = float(initial_gaps.min()*1e3)
    frame = local_frame()
    mesh_rows = []
    for tilt in TILTS:
        for azimuth in AZIMUTHS:
            mesh_rows.append(mesh_pose(samples, sample_region, wall, frame, cad_com, tilt, azimuth))
        print("MESH_REGRESSION_TILT_COMPLETE=%s" % tilt, flush=True)
    mesh_map = pd.DataFrame(mesh_rows)
    cad_map = pd.read_csv(HERE/"L2300_static_clearance_map.csv")
    regression = cad_map.merge(mesh_map, on=["tilt_deg", "azimuth_deg"], validate="one_to_one")
    for region in ("", "HEAD_", "BODY_", "TAIL_"):
        cad_column = (region+"min_gap_um") if region else "min_gap_um"
        mesh_column = "mesh_"+cad_column
        regression[(region+"gap_error_um") if region else "gap_error_um"] = regression[mesh_column]-regression[cad_column]
    regression["CAD_contact"] = regression.min_gap_um <= 0.0
    regression["contact_classification_match"] = regression.CAD_contact == regression.mesh_contact
    regression["wall_support_classification_match"] = regression.WALL_SUPPORT == regression.mesh_WALL_SUPPORT
    regression["opposing_classification_match"] = regression.opposing_bridge_20um == regression.mesh_opposing_bridge_20um
    regression.to_csv(HERE/"L2300_CAD_vs_solver_surface_regression.csv", index=False)

    max_gap_error = float(regression.gap_error_um.abs().max())
    max_region_error = float(regression[["HEAD_gap_error_um", "BODY_gap_error_um", "TAIL_gap_error_um"]].abs().to_numpy().max())
    gates = {
        "element_count_20000_to_50000_preferred": bool(20000 <= len(elements) <= 50000),
        "element_count_le_80000_hard": bool(len(elements) <= 80000),
        "mass_error_lt_0p5pct": bool(mass_error < .005),
        "COM_error_lt_5um": bool(com_error < .005),
        "principal_inertia_error_lt_1p5pct": bool(inertia_error.max() < .015),
        "initial_gap_no_penetration": bool(initial_gap_um >= 0.0),
        "head_tail_polarity": bool(np.dot(rotation@np.array([1.0, 0.0, 0.0]), A0) > .999999),
        "overall_gap_error_le_10um": bool(max_gap_error <= 10.0),
        "regional_gap_error_le_10um": bool(max_region_error <= 10.0),
        "contact_classification_all_match": bool(regression.contact_classification_match.all()),
        "wall_support_classification_all_match": bool(regression.wall_support_classification_match.all()),
        "opposing_classification_all_match": bool(regression.opposing_classification_match.all()),
    }
    hard_gates = {key: value for key, value in gates.items()
                  if key != "element_count_20000_to_50000_preferred"}
    summary = {
        "source_STEP_SHA256": raw["source_STEP_SHA256"],
        "mesh_strategy": "global 0.060 mm; no local refinement",
        "node_count": len(local), "element_count": len(elements),
        "exterior_triangle_count": len(faces), "surface_sample_count": len(samples),
        "mesh_volume_mm3": volume, "mesh_mass_mg": mass*1e9,
        "mass_relative_error": mass_error, "mesh_COM_mm": com.tolist(),
        "COM_error_um": com_error*1e3,
        "mesh_principal_inertia_tonne_mm2": principal.tolist(),
        "inertia_relative_error": inertia_error.tolist(),
        "initial_exact_gap_um": initial_gap_um,
        "surface_normal_diagnostic": {
            "overall_P95_deg": float(np.quantile(angles, .95)),
            "HEAD_P95_deg": float(normal_table.loc[normal_table.region == "HEAD", "angle_error_deg"].quantile(.95)),
            "hard_gate": False,
        },
        "CAD_vs_solver": {
            "maximum_overall_gap_error_um": max_gap_error,
            "maximum_regional_gap_error_um": max_region_error,
            "P95_overall_gap_error_um": float(regression.gap_error_um.abs().quantile(.95)),
            "contact_mismatch_count": int((~regression.contact_classification_match).sum()),
            "wall_support_mismatch_count": int((~regression.wall_support_classification_match).sum()),
            "opposing_mismatch_count": int((~regression.opposing_classification_match).sum()),
        },
        "gates": gates,
        "all_hard_gates_pass": bool(all(hard_gates.values())),
        "HEAD_refinement_required": bool(not all(hard_gates.values())),
    }
    (HERE/"L2300_solver_mesh_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    pd.DataFrame([{
        "source_STEP_SHA256": raw["source_STEP_SHA256"], "nodes": len(local),
        "elements": len(elements), "element_type": "C3D4", "global_size_mm": .060,
        "HEAD_local_refinement": False, "initial_gap_um": initial_gap_um,
    }]).to_csv(HERE/"L2300_solver_mesh_identity.csv", index=False)
    pd.DataFrame([{"quantity": "mass_mg", "CAD": mass_target*1e9, "mesh": mass*1e9,
                   "relative_error": mass_error, "gate": mass_error < .005},
                  {"quantity": "COM_error_um", "CAD": 0.0, "mesh": com_error*1e3,
                   "relative_error": np.nan, "gate": com_error < .005}]
                 + [{"quantity": "principal_inertia_%d_tonne_mm2" % (i+1),
                     "CAD": target[i], "mesh": principal[i], "relative_error": inertia_error[i],
                     "gate": inertia_error[i] < .015} for i in range(3)]).to_csv(
        HERE/"L2300_mass_inertia_audit.csv", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
