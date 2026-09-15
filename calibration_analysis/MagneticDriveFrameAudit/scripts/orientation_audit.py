"""Summarize actual ODB pose mapping and independent synthetic pose checks."""
from pathlib import Path
import json
import math

import numpy as np
import pandas as pd

from audit_common import AUDIT, build_production_model, rodrigues, rotation_error_deg, unit


DATA = AUDIT / "field_only"


def quaternion_matrix(axis, angle):
    axis = unit(axis)
    half = 0.5 * angle
    w = math.cos(half)
    x, y, z = axis * math.sin(half)
    return np.array([
        [1 - 2 * (y*y + z*z), 2 * (x*y - z*w), 2 * (x*z + y*w)],
        [2 * (x*y + z*w), 1 - 2 * (x*x + z*z), 2 * (y*z - x*w)],
        [2 * (x*z - y*w), 2 * (y*z + x*w), 1 - 2 * (x*x + y*y)]])


def main():
    actual = pd.read_csv(DATA / "orientation_mapping_actual_odb.csv")
    production, _, _, _ = build_production_model()
    axis = unit(np.array([0.37, -0.51, 0.776]))
    rows = []
    for angle_deg in (0, 15, 30, 45, 60, 90):
        angle = math.radians(angle_deg)
        server = production.MagneticCouplingModel.rotation_matrix_from_ur(axis * angle)
        reference = quaternion_matrix(axis, angle)
        rows.append({"equivalent_pose_deg": angle_deg, "rotation_axis_x": axis[0], "rotation_axis_y": axis[1],
                     "rotation_axis_z": axis[2], "R_server_vs_quaternion_error_deg": rotation_error_deg(server, reference)})
    synthetic = pd.DataFrame(rows)
    synthetic.to_csv(DATA / "orientation_mapping_synthetic_poses.csv", index=False)
    summary = {
        "abaqus_UR_interpretation": "total axis-angle rotation vector in radians for the rigid-body RP",
        "server_conversion": "Rodrigues exponential map, not XYZ Euler angles",
        "actual_odb_frame_count": int(len(actual)),
        "actual_odb_max_rotation_deg": float(actual.UR_norm_deg.max()),
        "actual_odb_max_server_vs_rigid_nodes_error_deg": float(actual.R_server_vs_R_nodes_error_deg.max()),
        "actual_odb_max_kabsch_rms_mm": float(actual.kabsch_rms_mm.max()),
        "synthetic_pose_angles_deg": synthetic.equivalent_pose_deg.tolist(),
        "synthetic_max_server_vs_quaternion_error_deg": float(synthetic.R_server_vs_quaternion_error_deg.max()),
        "classification": "POSE_MAPPING_VALID" if actual.R_server_vs_R_nodes_error_deg.max() < 0.1 else "MAGNETIC_POSE_MAPPING_BUG",
    }
    (DATA / "orientation_mapping_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
