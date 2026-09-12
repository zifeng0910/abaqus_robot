"""Replay the prior jam poses with the shorter rigid mesh as an eligibility gate."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
BASE = REPO / "calibration_analysis" / "ReducedHydro_zeta050_8p333_validation"
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
sys.path.insert(0, str(AUDIT / "normal_contact_damping_probe"))
sys.path.insert(0, str(BASE))
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall
from contact_probe_common import dense_rp
from analyze_zeta050_8p333 import angular_clusters


JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"


def main():
    identity = json.loads((HERE / "L1800_input_identity.json").read_text())
    t, data = dense_rp(BASE / "candidate_private")
    mesh = mesh_properties((HERE / f"{JOB}.inp").read_text())[0]["Robot_SOLID"]
    rp = np.asarray(identity["RP_abaqus_mm"], dtype=float)
    rotation = Rotation.from_rotvec(data["UR"])
    faces = pd.read_csv(ROOT / "Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv")
    node_ids = np.unique(faces[["n1", "n2", "n3"]].to_numpy()).astype(int)
    nodes = np.array([mesh["nodes"][int(i)] for i in node_ids]) + mesh["shift"]
    frame = pd.read_csv(BASE / "zeta050_8p333_true_axis.csv")
    wall = Wall()
    sample_times = np.unique(np.r_[np.arange(0, t[-1] + 1e-12, 2.5e-5),
                                   0.0017712, 0.0018637, 0.0027727, 0.008333])
    indices = np.unique(np.clip(np.rint(sample_times / 1e-7).astype(int), 0, len(t) - 1))
    rows = []
    for i in indices:
        positions = rp + data["U"][i] + rotation[i].apply(nodes - rp)
        gaps, triangles, _ = wall.query(positions)
        e1 = frame.loc[i, ["e1_x", "e1_y", "e1_z"]].to_numpy(float)
        e2 = frame.loc[i, ["e2_x", "e2_y", "e2_z"]].to_numpy(float)
        for threshold_um in (20.0, 5.0, 0.0):
            selected = np.where(gaps * 1e3 <= threshold_um)[0]
            normals = wall.normals[triangles[selected]] if len(selected) else np.empty((0, 3))
            angles = np.arctan2(normals @ e2, normals @ e1) if len(selected) else np.array([])
            centers = angular_clusters(angles)
            separation = max((np.degrees(np.arccos(np.clip(np.cos(a - b), -1, 1)))
                              for q, a in enumerate(centers) for b in centers[q + 1:]), default=0.0)
            rows.append({
                "time_s": t[i], "increment": int(i), "threshold_um": threshold_um,
                "min_gap_um": float(gaps.min() * 1e3), "node_count": len(selected),
                "cluster_count": len(centers), "max_sector_separation_deg": separation,
                "opposing_bridge": bool(len(centers) >= 2 and separation >= 120.0),
            })
    result = pd.DataFrame(rows)
    result.to_csv(HERE / "L1800_old_pose_exact_wall_gate.csv", index=False)
    at20 = result[result.threshold_um == 20.0]
    wall_span = 2.044327
    theta = 34.32513868122562
    envelope = (identity["candidate_pca_axial_span_mm"] * np.sin(np.radians(theta)) +
                identity["candidate_pca_transverse_diameter_mm"] * np.cos(np.radians(theta)))
    summary = {
        "sampled_old_pose_count": int(len(indices)),
        "sampling_interval_us": 25.0,
        "includes_prior_bridge_boundaries": True,
        "minimum_replayed_gap_um": float(result.min_gap_um.min()),
        "opposing_bridge_samples_20um": int(at20.opposing_bridge.sum()),
        "opposing_bridge_samples_5um": int(result[result.threshold_um == 5.0].opposing_bridge.sum()),
        "opposing_bridge_samples_0um": int(result[result.threshold_um == 0.0].opposing_bridge.sum()),
        "conservative_wall_span_mm": wall_span,
        "envelope_at_prior_max_tilt_mm": float(envelope),
        "analytic_clearance_at_prior_max_tilt_um": float((wall_span - envelope) * 1e3),
        "eligible_for_one_dynamic_run": bool(not at20.opposing_bridge.any() and envelope < wall_span),
        "scope": "kinematic old-pose geometry replay only; not a prediction of the new dynamic trajectory",
    }
    (HERE / "L1800_geometry_gate_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
