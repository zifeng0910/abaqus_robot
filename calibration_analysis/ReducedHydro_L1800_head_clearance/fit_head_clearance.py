"""Fit one smooth circular-arc head to the successful old L1800 envelope."""
from collections import Counter
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OLD = REPO / "calibration_analysis" / "ReducedHydro_geometry_L1800_validation"
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
from audit_stage_a import mesh_properties

JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
AXIS = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
AXIS /= np.linalg.norm(AXIS)
R_BODY = 0.4075
CURRENT_R_HEAD = 0.4142916666666666
FIT_LIMIT = 0.42
BIN_WIDTH = 0.005


def exterior_ids(mesh):
    count = Counter()
    for _, a, b, c, d in mesh["elems"]:
        for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            count[tuple(sorted(face))] += 1
    return np.unique([n for face, number in count.items() if number == 1 for n in face])


def circular_profile(x, radius):
    extent = np.sqrt(2.0 * radius * R_BODY - R_BODY**2)
    return np.where(
        x < extent,
        R_BODY - radius + np.sqrt(np.maximum(0.0, radius**2 - (x - extent) ** 2)),
        R_BODY,
    )


def main():
    meshes, _, _ = mesh_properties((OLD / f"{JOB}.inp").read_text())
    mesh = meshes["Robot_SOLID"]
    ids = exterior_ids(mesh)
    points = np.asarray([mesh["nodes"][int(i)] for i in ids]) + mesh["shift"]
    q = points - mesh["com"]
    axial = q @ AXIS
    radius = np.linalg.norm(q - np.outer(axial, AXIS), axis=1)
    axial -= axial.min()

    edges = np.arange(0.0, FIT_LIMIT + BIN_WIDTH * 1.01, BIN_WIDTH)
    centers = (edges[:-1] + edges[1:]) / 2.0
    envelope = np.full(len(centers), np.nan)
    bins = np.digitize(axial, edges) - 1
    for i in range(len(centers)):
        values = radius[bins == i]
        if len(values):
            envelope[i] = values.max()
    valid = np.isfinite(envelope)
    envelope = np.interp(centers, centers[valid], envelope[valid])

    # A dense deterministic one-dimensional search avoids optimizer dependence.
    candidates = np.linspace(R_BODY / 2.0 + 1e-4, 0.8, 60001)
    rows = []
    for head_radius in candidates:
        fitted = circular_profile(centers, head_radius)
        difference_um = (fitted - envelope) * 1e3
        maximum = float(difference_um.max())
        rmse = float(np.sqrt(np.mean(difference_um**2)))
        if maximum <= 10.0:
            rows.append((rmse, head_radius, maximum, float(difference_um.min())))
    if not rows:
        raise RuntimeError("No one-arc candidate satisfies the 10 um outward-excess target")
    rmse, fitted_radius, maximum, minimum = min(rows)
    extent = float(np.sqrt(2.0 * fitted_radius * R_BODY - R_BODY**2))
    fitted = circular_profile(centers, fitted_radius)
    current = circular_profile(centers, CURRENT_R_HEAD)
    profile = pd.DataFrame({
        "nose_axial_x_mm": centers,
        "old_L1800_envelope_mm": envelope,
        "current_FreeCAD_mm": current,
        "HeadClearance_candidate_mm": fitted,
        "current_minus_old_um": (current - envelope) * 1e3,
        "candidate_minus_old_um": (fitted - envelope) * 1e3,
    })
    profile.to_csv(HERE / "head_clearance_profile_fit.csv", index=False)
    summary = {
        "method": "deterministic 1D grid fit of one G1 circular arc to old L1800 exterior-node envelope",
        "comparison_interval_mm": [0.0, FIT_LIMIT],
        "bin_width_mm": BIN_WIDTH,
        "L_total_mm": 1.8,
        "D_body_mm": 0.815,
        "R_body_mm": R_BODY,
        "current_R_head_mm": CURRENT_R_HEAD,
        "fitted_R_head_design_mm": float(fitted_radius),
        "R_head_over_R_body": float(fitted_radius / R_BODY),
        "head_axial_extent_mm": extent,
        "straight_body_length_mm": 1.8 - extent,
        "max_outward_excess_um": maximum,
        "minimum_difference_um": minimum,
        "profile_rmse_um": rmse,
        "acceptance_max_outward_um": 10.0,
        "fit_gate_pass": bool(maximum <= 10.0),
    }
    (HERE / "head_clearance_fit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
