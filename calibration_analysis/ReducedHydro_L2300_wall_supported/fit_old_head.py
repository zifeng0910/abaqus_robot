"""Reconstruct the successful old L1800 head envelope and audit one smooth fit."""
from collections import Counter
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from scipy.special import comb


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OLD = REPO / "calibration_analysis" / "ReducedHydro_geometry_L1800_validation"
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
from audit_stage_a import mesh_properties

OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
AXIS = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
AXIS /= np.linalg.norm(AXIS)
R_BODY = 0.4075
# Deterministic degree-5 Bezier fit obtained against the continuous upper
# convex hull below. The final three poles are radially collinear, imposing a
# zero-curvature (G2) transition into the cylinder.
POLES = np.array([
    [0.0, 0.0],
    [0.0, 0.09664525],
    [0.005, 0.18497527],
    [0.17271875, R_BODY],
    [0.25, R_BODY],
    [0.26149477, R_BODY],
])


def exterior_ids(mesh):
    count = Counter()
    for _, a, b, c, d in mesh["elems"]:
        for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            count[tuple(sorted(face))] += 1
    return np.unique([node for face, number in count.items() if number == 1 for node in face])


def bezier(t):
    degree = len(POLES) - 1
    basis = np.column_stack([
        comb(degree, i) * (1.0 - t) ** (degree - i) * t**i
        for i in range(degree + 1)
    ])
    return basis @ POLES


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    mesh = mesh_properties((OLD / f"{OLD_JOB}.inp").read_text())[0]["Robot_SOLID"]
    ids = exterior_ids(mesh)
    points = np.asarray([mesh["nodes"][int(i)] for i in ids]) + mesh["shift"]
    q = points - mesh["com"]
    axial = q @ AXIS
    axial -= axial.min()
    radial = np.linalg.norm(q - np.outer(q @ AXIS, AXIS), axis=1)

    projected = np.column_stack([axial, radial])
    hull = projected[ConvexHull(projected).vertices]
    upper = hull[(hull[:, 0] <= 0.30 + 1e-12) & (hull[:, 1] >= 0.0)]
    upper = upper[np.argsort(upper[:, 0])]
    # Remove lower-return points of the 2-D hull, retaining its monotone head boundary.
    upper = upper[np.r_[True, np.diff(upper[:, 1]) >= -1e-10]]
    upper = upper[upper[:, 1] <= R_BODY + 1e-6]
    if len(upper) < 6:
        raise RuntimeError("Old exterior hull did not provide a resolvable head envelope")

    x = np.linspace(0.0, POLES[-1, 0], 2001)
    old = np.interp(x, upper[:, 0], upper[:, 1])
    curve = bezier(np.linspace(0.0, 1.0, 20001))
    if np.any(np.diff(curve[:, 0]) < -1e-12) or np.any(np.diff(curve[:, 1]) < -1e-12):
        raise RuntimeError("Bezier head is not monotone")
    fitted = np.interp(x, curve[:, 0], curve[:, 1])
    difference_um = (fitted - old) * 1e3
    dx = np.gradient(curve[:, 0]); dr = np.gradient(curve[:, 1])
    slope = dr / np.maximum(dx, 1e-14)
    finite = np.isfinite(slope) & (curve[:, 0] > 1e-5)
    convex = bool(np.max(np.diff(slope[finite])) <= 1e-3)
    summary = {
        "source": "actual exterior nodes of the successful old scaled L1800 C3D4 mesh",
        "envelope_method": "continuous piecewise-linear upper boundary of the axial-radial convex hull",
        "fit_family": "degree-5 Bezier (single smooth BSpline span, G2 shoulder)",
        "degree": 5,
        "poles_mm": POLES.tolist(),
        "head_axial_extent_mm": float(POLES[-1, 0]),
        "comparison_interval_mm": [0.0, float(POLES[-1, 0])],
        "comparison_samples": len(x),
        "old_exterior_node_count": len(ids),
        "old_head_hull_points": upper.tolist(),
        "max_outward_excess_um": float(difference_um.max()),
        "minimum_difference_um": float(difference_um.min()),
        "profile_rmse_um": float(np.sqrt(np.mean(difference_um**2))),
        "monotone_axial": bool(np.all(np.diff(curve[:, 0]) >= -1e-12)),
        "monotone_radial": bool(np.all(np.diff(curve[:, 1]) >= -1e-12)),
        "convex_profile": convex,
        "nose_tangent": "radial (smooth axisymmetric tip)",
        "shoulder_tangent": "axial with zero endpoint curvature (G2 with cylinder)",
    }
    summary["preferred_max_outward_le_5um"] = bool(summary["max_outward_excess_um"] <= 5.0)
    summary["hard_max_outward_le_10um"] = bool(summary["max_outward_excess_um"] <= 10.0)
    summary["preferred_RMS_le_5um"] = bool(summary["profile_rmse_um"] <= 5.0)
    summary["fit_gate_pass"] = bool(
        summary["max_outward_excess_um"] <= 10.0
        and summary["monotone_axial"] and summary["monotone_radial"] and convex
    )
    pd.DataFrame({
        "nose_axial_x_mm": x,
        "old_L1800_convex_envelope_mm": old,
        "L2300_Bezier_head_mm": fitted,
        "new_minus_old_um": difference_um,
    }).to_csv(HERE / "L2300_head_profile_fit.csv", index=False)
    (HERE / "L2300_head_fit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not summary["fit_gate_pass"]:
        raise RuntimeError("Head fit gate failed")


if __name__ == "__main__":
    main()
