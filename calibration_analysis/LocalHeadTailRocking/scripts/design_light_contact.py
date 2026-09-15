"""Derive one light-contact rocking amplitude from exact geometry and ROCK10 gain."""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
JOB = "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
CASE = OUT / "case" / JOB
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
OVERTRAVEL_MM = 0.010
MAX_OVERTRAVEL_DEG = 0.5


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def part_nodes(deck, name):
    part = re.search(rf"(?ms)^\*Part, name={re.escape(name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    return np.asarray([[float(value) for value in row.split(",")[1:4]] for row in block.splitlines() if row.strip()])


def wall_planes(wall, center, c, n, b):
    axial = (wall - center).dot(c)
    ring = wall[np.abs(axial - axial.min()) < 1e-7]
    xy = np.column_stack(((ring - center).dot(n), (ring - center).dot(b)))
    xy = xy[np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))]
    normals, radii = [], []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p
        normal2 = unit([edge[1], -edge[0]])
        if np.dot(normal2, p + q) < 0:
            normal2 *= -1
        normals.append(normal2[0] * n + normal2[1] * b)
        radii.append(np.dot(p, normal2))
    radii = np.asarray(radii)
    if np.ptp(radii) > 1e-8:
        raise RuntimeError(f"Wall face radii are inconsistent: range={np.ptp(radii):.3e} mm")
    return np.asarray(normals), float(radii.mean())


def harmonic_fit(time, values, frequency=5.0):
    omega = 2 * np.pi * frequency
    matrix = np.column_stack((np.sin(omega * time), np.cos(omega * time), np.ones_like(time)))
    sin_coef, cos_coef, mean = np.linalg.lstsq(matrix, values, rcond=None)[0]
    amplitude = float(np.hypot(sin_coef, cos_coef))
    phase = math.degrees(math.atan2(cos_coef, sin_coef))
    return {"amplitude": amplitude, "phase_deg": phase, "mean": float(mean)}


def bisection(func, start, direction, target=0.0, initial_step=0.25, limit=30.0):
    inside = float(start)
    if func(inside) <= target:
        raise RuntimeError("Bisection must begin above the target gap")
    outside = inside + direction * initial_step
    while func(outside) > target and abs(outside - start) < limit:
        outside += direction * initial_step
    if func(outside) > target:
        raise RuntimeError(f"Gap target {target} mm not reached within {limit} deg")
    for _ in range(60):
        midpoint = 0.5 * (inside + outside)
        if func(midpoint) > target:
            inside = midpoint
        else:
            outside = midpoint
    return 0.5 * (inside + outside)


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    if identity["status"] != "SOLVED":
        raise RuntimeError(f"ROCK10 status is {identity['status']}, expected SOLVED")
    deck_path = CASE / f"{JOB}.inp"
    deck = deck_path.read_text()
    nodes = part_nodes(deck, "Robot_SOLID")
    wall = part_nodes(deck, "Pipe_WALL_HELPER")
    c, n, b = map(unit, (identity["initial_axis_aba"], identity["n_rock_aba"], identity["b_rock_aba"]))
    normals, radius = wall_planes(wall, RP0, c, n, b)
    axial = (nodes - RP0).dot(c)
    low, high = axial.min(), axial.max()
    region = np.full(len(nodes), "BODY", dtype=object)
    region[axial <= low + 0.25] = "HEAD"
    region[axial >= high - 0.25] = "TAIL"

    def gaps(theta_deg):
        rotation = Rotation.from_rotvec(math.radians(theta_deg) * b)
        points = RP0 + rotation.apply(nodes - RP0)
        relative = points - RP0
        radial = relative - np.outer(relative.dot(c), c)
        face_gap = radius - radial.dot(normals.T)
        node_gap = face_gap.min(axis=1)
        return {
            "overall": float(node_gap.min()),
            "HEAD": float(node_gap[region == "HEAD"].min()),
            "TAIL": float(node_gap[region == "TAIL"].min()),
            "BODY": float(node_gap[region == "BODY"].min()),
        }

    theta_touch_pos = bisection(lambda angle: gaps(angle)["overall"], 0.0, 1.0)
    theta_touch_neg = bisection(lambda angle: gaps(angle)["overall"], 0.0, -1.0)
    theta_10um_pos = bisection(lambda angle: gaps(angle)["overall"], theta_touch_pos, 1.0, -OVERTRAVEL_MM)
    theta_10um_neg = bisection(lambda angle: gaps(angle)["overall"], theta_touch_neg, -1.0, -OVERTRAVEL_MM)
    theta_target_pos = min(theta_10um_pos, theta_touch_pos + MAX_OVERTRAVEL_DEG)
    theta_target_neg = max(theta_10um_neg, theta_touch_neg - MAX_OVERTRAVEL_DEG)
    target_pos_gaps = gaps(theta_target_pos)
    target_neg_gaps = gaps(theta_target_neg)

    pose = pd.read_csv(OUT / "rocking_pose.csv")
    time = pose.time_s.to_numpy(float)
    theta = pose.theta_rock_deg.to_numpy(float)
    alpha = pose.alpha_B_deg.to_numpy(float)
    robot_fit = harmonic_fit(time, theta)
    field_fit = harmonic_fit(time, alpha)
    gain = robot_fit["amplitude"] / field_fit["amplitude"]
    phase_lag = (robot_fit["phase_deg"] - field_fit["phase_deg"] + 180.0) % 360.0 - 180.0
    theta0 = robot_fit["mean"]
    required_pos = (theta_target_pos - theta0) / gain
    required_neg = (theta0 - theta_target_neg) / gain
    selected = max(required_pos, required_neg)
    predicted_pos = theta0 + gain * selected
    predicted_neg = theta0 - gain * selected

    if selected > 20.0:
        classification = "LIGHT_CONTACT_REQUIRES_UNEXPECTEDLY_LARGE_FIELD_ANGLE"
    else:
        classification = "LIGHT_CONTACT_AMPLITUDE_SELECTED"

    touch_pos_gaps = gaps(theta_touch_pos)
    touch_neg_gaps = gaps(theta_touch_neg)
    touch_pos_region = min((key for key in ("HEAD", "TAIL", "BODY")), key=touch_pos_gaps.get)
    touch_neg_region = min((key for key in ("HEAD", "TAIL", "BODY")), key=touch_neg_gaps.get)

    angle = np.linspace(theta_touch_neg - 1.0, theta_touch_pos + 1.0, 1601)
    gap_table = pd.DataFrame({"theta_rock_deg": angle})
    for key in ("HEAD", "TAIL", "BODY", "overall"):
        gap_table[key + "_gap_um"] = [gaps(value)[key] * 1e3 for value in angle]
    gap_table.to_csv(OUT / "ExactGap_vs_RockingAngle.csv", index=False)

    fig, axis_plot = plt.subplots(figsize=(8.2, 5.2))
    colors = {"HEAD": "#c43c35", "TAIL": "#315a8a", "BODY": "#555555", "overall": "#237a57"}
    for key in ("HEAD", "TAIL", "BODY", "overall"):
        axis_plot.plot(angle, gap_table[key + "_gap_um"], label=key, color=colors[key], linewidth=2.2 if key == "overall" else 1.3)
    axis_plot.axhline(0.0, color="black", linewidth=0.8)
    axis_plot.axhline(-10.0, color="#777777", linewidth=0.8, linestyle="--")
    axis_plot.axvline(theta_touch_neg, color="#315a8a", linewidth=0.9, linestyle=":")
    axis_plot.axvline(theta_touch_pos, color="#c43c35", linewidth=0.9, linestyle=":")
    axis_plot.scatter(
        [theta_target_neg, theta_target_pos],
        [target_neg_gaps["overall"] * 1e3, target_pos_gaps["overall"] * 1e3],
        color="#222222", zorder=5, label="light-contact target",
    )
    axis_plot.set(xlabel="Rocking angle (deg)", ylabel="Exact faceted-wall gap (um)", title="Exact first touch and conservative light-contact target")
    axis_plot.grid(alpha=0.2)
    axis_plot.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "ExactGap_vs_RockingAngle.png", dpi=220)
    fig.savefig(OUT / "ExactGap_vs_RockingAngle.pdf")
    plt.close(fig)

    metrics = json.loads((OUT / "rocking_metrics.json").read_text())
    manifest = {
        "classification": classification,
        "geometry_method": "exact robot mesh vertices against every actual faceted-wall half-space",
        "source_input_sha256": sha256(deck_path),
        "solver_wall_face_count": int(len(normals)),
        "solver_wall_inscribed_radius_mm": radius,
        "theta_touch_pos_deg": theta_touch_pos,
        "theta_touch_neg_deg": theta_touch_neg,
        "first_touch_pos_region": touch_pos_region,
        "first_touch_neg_region": touch_neg_region,
        "gaps_at_touch_pos_um": {key: value * 1e3 for key, value in touch_pos_gaps.items()},
        "gaps_at_touch_neg_um": {key: value * 1e3 for key, value in touch_neg_gaps.items()},
        "theta_10um_free_overtravel_pos_deg": theta_10um_pos,
        "theta_10um_free_overtravel_neg_deg": theta_10um_neg,
        "overtravel_cap_deg": MAX_OVERTRAVEL_DEG,
        "theta_target_pos_deg": theta_target_pos,
        "theta_target_neg_deg": theta_target_neg,
        "gaps_at_target_pos_um": {key: value * 1e3 for key, value in target_pos_gaps.items()},
        "gaps_at_target_neg_um": {key: value * 1e3 for key, value in target_neg_gaps.items()},
        "theta0_deg": theta0,
        "A_robot_current_deg": robot_fit["amplitude"],
        "A_field_current_deg": field_fit["amplitude"],
        "G_theta": gain,
        "phase_lag_robot_vs_field_deg": phase_lag,
        "actual_positive_peak_deg": float(theta.max()),
        "actual_negative_peak_deg": float(theta.min()),
        "final_theta_rock_deg": float(theta[-1]),
        "theta_rock_mean_sampled_deg": float(theta.mean()),
        "required_field_amplitude_from_positive_deg": required_pos,
        "required_field_amplitude_from_negative_deg": required_neg,
        "selected_A_B_required_deg": selected,
        "predicted_positive_robot_peak_deg": predicted_pos,
        "predicted_negative_robot_peak_deg": predicted_neg,
        "selection_basis": "minimum symmetric field amplitude satisfying both asymmetric targets; each target uses the smaller of 10 um free overtravel and 0.5 deg angular overtravel",
        "rock10_validation": {
            "classification": metrics["classification"],
            "theta_cross_rms_deg": metrics["theta_cross_rms_deg"],
            "theta_cross_max_abs_deg": metrics["theta_cross_max_abs_deg"],
            "pca_rocking_linearity": metrics["pca_rocking_linearity"],
            "COM_radial_max_um": metrics["COM_radial_max_um"],
            "axial_displacement_mm": metrics["axial_displacement_mm"],
            "HEAD_only_support_fraction": metrics["HEAD_only_support_fraction"],
            "TAIL_only_support_fraction": metrics["TAIL_only_support_fraction"],
            "both_end_support_fraction": metrics["both_end_support_fraction"],
            "no_contact_fraction": metrics["no_contact_fraction"],
            "HEAD_min_gap_um": metrics["HEAD_min_gap_um"],
            "TAIL_min_gap_um": metrics["TAIL_min_gap_um"],
            "BODY_min_gap_um": metrics["BODY_min_gap_um"],
            "energy_ranges": {key: metrics[key] for key in metrics if key.endswith("_min") or key.endswith("_max")},
        },
    }
    (OUT / "Selected_LightContact_Amplitude.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    if selected > 20.0:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
