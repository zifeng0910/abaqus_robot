"""Run the mesh-regression and unique HeadClearance zero-Abaqus replays."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FREECAD = REPO / "calibration_analysis" / "ReducedHydro_FreeCAD_L1800_validation"
OLD = REPO / "calibration_analysis" / "ReducedHydro_geometry_L1800_validation"
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path[:0] = [str(FREECAD), str(AUDIT)]
from audit_stage_a import mesh_properties
from build_and_audit_freecad import rotation_x_to_axis
from cross_geometry_trajectory_replay import (DT, THRESHOLDS, Wall, clusters,
                                               exact_query, geometry, replay,
                                               summarize, trajectory)

AXIS = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
AXIS /= np.linalg.norm(AXIS)
RP = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
CURRENT_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083"
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
MESH_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_Mesh060"


def center_tree():
    path = REPO.parent / "CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv"
    points = pd.read_csv(path)[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    return cKDTree(points)


def analytic_candidate():
    definition = json.loads((HERE / "Robot_parametric_L1p800_D0p815_HeadClearance_geometry.json").read_text())
    props = json.loads((HERE / "Robot_parametric_L1p800_D0p815_HeadClearance_mass_properties.json").read_text())
    rb = definition["R_body_mm"]; rh = definition["R_head_design_mm"]
    xh = definition["head_axial_extent_mm"]
    # Axial spacing <=15 um and 10-degree circumferential spacing keep the
    # maximum cylindrical chord sag below 1.6 um, within the 2 um replay gate.
    x = np.unique(np.r_[np.linspace(0.0, xh, int(np.ceil(xh/.015))+1),
                        np.linspace(xh, 1.8, int(np.ceil((1.8-xh)/.015))+1)])
    radial = np.where(x < xh, rb-rh+np.sqrt(np.maximum(0.0, rh*rh-(x-xh)**2)), rb)
    angles = np.linspace(0.0, 2.0*np.pi, 36, endpoint=False)
    local = np.asarray([[xx, rr*np.cos(a), rr*np.sin(a)] for xx, rr in zip(x, radial) for a in angles])
    local = np.vstack([local, [1.8, 0.0, 0.0]])
    com_local = np.asarray(props["center_of_mass"]["value"], float)
    rotation = rotation_x_to_axis()
    points = RP + (local-com_local) @ rotation.T
    axial = (points-RP) @ AXIS
    labels = np.where(axial <= np.quantile(axial, .10), "HEAD",
                      np.where(axial >= np.quantile(axial, .90), "TAIL", "BODY"))
    return {"name": "HEADCLEARANCE", "ids": np.arange(len(points)), "points": points,
            "labels": labels, "com": RP, "sampling": {"axial_max_spacing_um": 15.0,
            "circumferential_step_deg": 10.0, "maximum_chord_sag_um": 1.551,
            "point_count": len(points)}}


def replay_prefiltered(geom, traj, wall, tree):
    """Equivalent replay with the conservative centerline rejection applied first."""
    rotations = Rotation.from_rotvec(traj["UR"]).as_matrix()
    traj_com = traj["rp"] + traj["U"] + np.einsum(
        "bij,j->bi", rotations, traj["mesh"]["com"] - traj["rp"])
    local = geom["points"] - geom["com"]
    e1_all = traj["axis"][["e1_x", "e1_y", "e1_z"]].to_numpy(float)
    e2_all = traj["axis"][["e2_x", "e2_y", "e2_z"]].to_numpy(float)
    timeline = []
    for offset in range(0, len(traj["time"]), 250):
        stop = min(len(traj["time"]), offset + 250); rot = rotations[offset:stop]
        positions = traj_com[offset:stop, None, :] + np.einsum("bij,nj->bni", rot, local)
        flat = positions.reshape(-1, 3)
        center_distance = tree.query(flat, k=1)[0]
        eligible_flat = np.flatnonzero(center_distance >= .70)
        possible = np.zeros(flat.shape[0], bool)
        if len(eligible_flat):
            wall_distance = wall.tree.query(flat[eligible_flat], k=1)[0]
            possible[eligible_flat[wall_distance <= wall.radius + .020001]] = True
        rows, cols = np.where(possible.reshape(stop-offset, len(local)))
        exact_gap = np.full((stop-offset, len(local)), .020001)
        exact_tri = np.full((stop-offset, len(local)), -1, int)
        for chunk in range(0, len(rows), 5000):
            part = slice(chunk, chunk+5000)
            gap, tri, _ = exact_query(wall, positions[rows[part], cols[part]])
            exact_gap[rows[part], cols[part]] = gap; exact_tri[rows[part], cols[part]] = tri
        for local_i in range(stop-offset):
            global_i = offset + local_i; gaps = exact_gap[local_i]; tris = exact_tri[local_i]
            valid = tris >= 0
            minimum = float(gaps[valid].min()*1e3) if valid.any() else 20.001
            row = {"geometry": geom["name"], "trajectory": traj["name"],
                   "time_s": traj["time"][global_i], "minimum_gap_um": minimum}
            e1, e2 = e1_all[global_i], e2_all[global_i]
            for threshold in THRESHOLDS:
                selected = valid & (gaps*1e3 <= threshold); labels = geom["labels"][selected]
                angles = (np.arctan2(wall.normals[tris[selected]] @ e2,
                                     wall.normals[tris[selected]] @ e1)
                          if selected.any() else np.array([]))
                centers = clusters(angles)
                separation = max((np.degrees(np.arccos(np.clip(np.cos(a-b), -1, 1)))
                                  for q, a in enumerate(centers) for b in centers[q+1:]), default=0.0)
                tag = str(int(threshold)); row.update({
                    f"bridge_{tag}um": int(len(centers) >= 2 and separation >= 120),
                    f"head_count_{tag}um": int(np.sum(labels == "HEAD")),
                    f"tail_count_{tag}um": int(np.sum(labels == "TAIL")),
                    f"body_count_{tag}um": int(np.sum(labels == "BODY")),
                    f"sector_separation_{tag}um_deg": separation})
            timeline.append(row)
        print(f"{geom['name']} geometry + {traj['name']} trajectory: {stop}/{len(traj['time'])}", flush=True)
    return pd.DataFrame(timeline)


def mesh_regression(wall, tree):
    cached_identity = HERE / "mesh060_replay_identity.json"
    if cached_identity.exists() and (HERE / "mesh060_replay_summary.csv").exists():
        identity = json.loads(cached_identity.read_text())
        print("Reusing completed mesh060 replay regression")
        print(pd.read_csv(HERE / "mesh060_replay_summary.csv").to_string(index=False))
        return identity
    dense_meshes, dense_rp, _ = mesh_properties((FREECAD / f"{CURRENT_JOB}.inp").read_text())
    simple_meshes, _, _ = mesh_properties((HERE / f"{MESH_JOB}.inp").read_text())
    dense_mesh = dense_meshes["Robot_SOLID"]; simple_mesh = simple_meshes["Robot_SOLID"]
    dense_ids = np.unique(pd.read_csv(FREECAD / "freecad_robot_surface_triangles_exact.csv")[["n1","n2","n3"]].to_numpy())
    simple_ids = np.unique(pd.read_csv(HERE / "current_mesh060_surface_triangles.csv")[["n1","n2","n3"]].to_numpy())
    dense_geom = geometry("NEW", dense_mesh, dense_ids, AXIS); dense_geom["name"] = "DENSE035"
    simple_geom = geometry("NEW", simple_mesh, simple_ids, AXIS); simple_geom["name"] = "MESH060"
    traj = trajectory("NEW", FREECAD, dense_mesh, dense_rp, FREECAD / "freecad_8p333_true_axis.csv")
    dense = pd.read_csv(FREECAD / "cross_replay_case_new_new.csv").copy()
    dense["geometry"] = "DENSE035"
    simple = replay(simple_geom, traj, wall, tree)
    timeline = pd.concat([dense, simple], ignore_index=True)
    timeline.to_csv(HERE / "mesh060_replay_timeline.csv", index=False)
    result = summarize(timeline); result.to_csv(HERE / "mesh060_replay_summary.csv", index=False)
    a, b = result.iloc[0], result.iloc[1]
    persistent_a = a.longest20_ms >= .5; persistent_b = b.longest20_ms >= .5
    gates = {"min_gap_difference_le_2um": abs(a.penetration_um-b.penetration_um) <= 2.0,
             "longest20_difference_le_0p05ms": abs(a.longest20_ms-b.longest20_ms) <= .05,
             "bridge_classification_consistent": bool(persistent_a == persistent_b)}
    identity = {"comparison": "dense 0.035 versus simplified 0.060 on current exact-CAD trajectory",
                "gates": {k: bool(v) for k,v in gates.items()}, "all_gates_pass": all(gates.values())}
    (HERE / "mesh060_replay_identity.json").write_text(json.dumps(identity, indent=2)+"\n")
    print(result.to_string(index=False)); print(json.dumps(identity, indent=2))
    return identity


def candidate_replay(wall, tree):
    old_meshes, old_rp, _ = mesh_properties((OLD / f"{OLD_JOB}.inp").read_text())
    new_meshes, new_rp, _ = mesh_properties((FREECAD / f"{CURRENT_JOB}.inp").read_text())
    old_mesh = old_meshes["Robot_SOLID"]; new_mesh = new_meshes["Robot_SOLID"]
    old_traj = trajectory("OLD", OLD, old_mesh, old_rp, OLD / "L1800_8p333_true_axis.csv")
    new_traj = trajectory("NEW", FREECAD, new_mesh, new_rp, FREECAD / "freecad_8p333_true_axis.csv")
    geom = analytic_candidate()
    tables = []
    for traj in (old_traj, new_traj):
        table = replay_prefiltered(geom, traj, wall, tree)
        table.to_csv(HERE / f"head_clearance_replay_{traj['name'].lower()}_private.csv", index=False)
        tables.append(table)
    timeline = pd.concat(tables, ignore_index=True)
    timeline.to_csv(HERE / "head_clearance_replay_timeline.csv", index=False)
    result = summarize(timeline); result.to_csv(HERE / "head_clearance_replay_summary.csv", index=False)
    gates = {}
    for _, row in result.iterrows():
        key = row.trajectory.lower()
        gates[f"{key}_longest20_lt_0p5ms"] = row.longest20_ms < .5
        gates[f"{key}_longest5_no_persistent"] = row.longest5_ms < .5
        gates[f"{key}_longest0_no_persistent"] = row.longest0_ms < .5
        gates[f"{key}_min_gap_gt_minus2um"] = row.penetration_um > -2.0
    identity = {"geometry": "exact analytic HeadClearance circular surface",
                "sampling": geom["sampling"], "wall": "SmoothWall114 exact triangles",
                "gates": {k: bool(v) for k,v in gates.items()},
                "all_zero_abaqus_gates_pass": all(gates.values())}
    (HERE / "head_clearance_replay_identity.json").write_text(json.dumps(identity, indent=2)+"\n")
    print(result.to_string(index=False)); print(json.dumps(identity, indent=2))
    return identity


def main():
    wall = Wall(); tree = center_tree()
    mesh = mesh_regression(wall, tree)
    if not mesh["all_gates_pass"]:
        raise RuntimeError("0.060 mesh replay gate failed; allowed 0.050 refinement is required")
    candidate = candidate_replay(wall, tree)
    if not candidate["all_zero_abaqus_gates_pass"]:
        raise RuntimeError("HeadClearance zero-Abaqus bridge gate failed; STOP")


if __name__ == "__main__":
    main()
