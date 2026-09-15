"""Analyze dense straight-control histories and compare them with curved cycle 1."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
BASE = REPO / "calibration_analysis" / "ProductionLocalFrameValidation"
CASE = OUT / "case" / "PROD_LOCAL30_G0_STRAIGHT_CTRL"
JOB = "PROD_LOCAL30_G0_STRAIGHT_CTRL"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600215763, -0.11887423721399708, 0.2348382536498942])
COM_X = 1.25606694407031
PERIOD = 1.0 / 30.0
sys.path.insert(0, str(REPO / "calibration_analysis" / "MagneticDriveFrameAudit" / "scripts"))
from analyze_diagnostics import contact_force, dense, intervals, merge_ranges


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / max(np.linalg.norm(vector), 1.0e-30)


def longest_duration(time, flag):
    return max((time[end] - time[start] + np.median(np.diff(time))
                for start, end in intervals(flag)), default=0.0)


def polygon_area(x, y):
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def surface_nodes():
    deck = (CASE / f"{JOB}.inp").read_text()
    part = re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    initial = np.asarray([[float(value) for value in line.split(",")[1:4]]
                          for line in block.splitlines() if line.strip()])
    axial = np.dot(initial - RP0, A0)
    region = np.full(len(initial), "BODY", dtype=object)
    region[axial <= axial.min() + 0.2615] = "HEAD"
    region[axial >= axial.max() - 0.02] = "TAIL"
    return initial, region


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    time, history = dense(CASE)
    dt = float(np.median(np.diff(time)))
    com_dense = RP0 + history["U"]
    rotation_dense = Rotation.from_rotvec(history["UR"])
    force = contact_force(CASE, time)
    force_magnitude = np.linalg.norm(force, axis=1)
    active = force_magnitude > 1.0e-8
    events = merge_ranges(intervals(active), int(round(5.0e-6 / dt)))

    center0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    tangent = unit(identity["straight_tangent_aba"])
    e1 = unit(identity["initial_e1_aba"])
    e2 = unit(identity["initial_e2_aba"])
    radius = float(identity["wall_radius_mm"])
    sample = np.unique(np.r_[np.arange(0, len(time), 100), len(time) - 1])
    ts = time[sample]
    com = com_dense[sample]
    rotations = rotation_dense[sample]
    axis = rotations.apply(np.broadcast_to(A0, (len(sample), 3)))
    s = np.dot(com - center0, tangent)
    centers = center0 + s[:, None] * tangent
    q1 = axis.dot(e1); q2 = axis.dot(e2); axial = axis.dot(tangent)
    tilt = np.degrees(np.arccos(np.clip(axial, -1.0, 1.0)))
    robot_phase = np.unwrap(np.arctan2(q2, q1))

    telemetry = pd.read_csv(CASE / f"{JOB}_telemetry.csv")
    field = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[name])
                             for name in ("Bx_aba_T", "By_aba_T", "Bz_aba_T")])
    b1 = field.dot(e1); b2 = field.dot(e2); bt = field.dot(tangent)
    field_phase = np.unwrap(np.arctan2(b2, b1))
    phase_lag = np.angle(np.exp(1j * (robot_phase - field_phase)))
    phase_valid = np.sin(np.radians(tilt)) > 0.1

    initial_nodes, region = surface_nodes()
    gap_rows = []
    for local_index, dense_index in enumerate(sample):
        points = com[local_index] + rotations[local_index].apply(initial_nodes - RP0)
        axial_point = np.dot(points - center0, tangent)
        projected = center0 + axial_point[:, None] * tangent
        radial = points - projected
        radial_distance = np.linalg.norm(radial, axis=1)
        gaps = radius - radial_distance
        nearest = int(np.argmin(gaps))
        sector = math.degrees(math.atan2(np.dot(radial[nearest], e2), np.dot(radial[nearest], e1))) % 360.0
        gap_rows.append({
            "time_s": ts[local_index], "sample_index": int(dense_index),
            "min_analytic_gap_um": float(gaps.min() * 1.0e3),
            "HEAD_gap_um": float(gaps[region == "HEAD"].min() * 1.0e3),
            "TAIL_gap_um": float(gaps[region == "TAIL"].min() * 1.0e3),
            "nearest_wall_sector_deg": sector,
            "minimum_end_margin_mm": float(identity["wall_half_length_mm"] - np.max(np.abs(axial_point))),
        })
    gap = pd.DataFrame(gap_rows)
    gap.to_csv(OUT / "straight_exact_gap_timeseries.csv", index=False)
    both_support = (gap.HEAD_gap_um <= 20.0) & (gap.TAIL_gap_um <= 20.0)

    event_rows = []
    separation_steps = int(round(20.0e-6 / dt))
    sector_bins = []
    for number, (start, end) in enumerate(events, 1):
        mid = (start + end) // 2
        direction = -np.mean(force[start:end + 1], axis=0)
        sector = math.degrees(math.atan2(np.dot(direction, e2), np.dot(direction, e1))) % 360.0
        sector_bins.append(int(round(sector / 30.0)) % 12)
        separated = bool(end + separation_steps < len(active) and
                         not active[end + 1:end + 1 + separation_steps].any())
        event_rows.append({
            "event": number, "start_s": float(time[start]), "end_s": float(time[end]),
            "duration_us": float((end - start + 1) * dt * 1.0e6),
            "peak_force_N": float(force_magnitude[start:end + 1].max()),
            "sector_deg": sector, "sector_bin_30deg": sector_bins[-1],
            "true_separation_after": int(separated),
        })
    event_table = pd.DataFrame(event_rows)
    event_table.to_csv(OUT / "straight_contact_events.csv", index=False)
    sector_switches = int(np.sum(np.asarray(sector_bins[1:]) != np.asarray(sector_bins[:-1]))) if len(sector_bins) > 1 else 0
    half_window = max(1, int(round(0.5 * np.median(np.diff(sample)))))
    contact_sample = np.asarray([active[max(0, i-half_window):min(len(active), i+half_window+1)].any() for i in sample])
    merged_active = np.zeros_like(active)
    for start, end in events:
        merged_active[start:end + 1] = True

    energy_archive = np.load(CASE / "private" / "energy_history_private.npz")
    energy_ranges = {}
    for name in energy_archive.files:
        values = energy_archive[name]
        energy_ranges[name] = float(np.ptp(values[:, 1]))
    valid_lag = phase_lag[phase_valid]
    mean_lag = math.degrees(math.atan2(np.mean(np.sin(valid_lag)), np.mean(np.cos(valid_lag))))
    metrics = {
        "case_id": JOB,
        "classification": "ONE_CYCLE_STRAIGHT_CONTROL",
        "local_B_winding": float((field_phase[-1] - field_phase[0]) / (2.0 * np.pi)),
        "robot_axis_winding": float((robot_phase[-1] - robot_phase[0]) / (2.0 * np.pi)),
        "mean_phase_lag_deg": mean_lag,
        "tilt_min_deg": float(tilt.min()), "tilt_mean_deg": float(tilt.mean()), "tilt_max_deg": float(tilt.max()),
        "orbit_area": polygon_area(q1, q2),
        "delta_s_mm": float(s[-1] - s[0]),
        "contact_event_count": len(events),
        "raw_contact_time_fraction": float(active.mean()),
        "merged_contact_time_fraction": float(merged_active.mean()),
        "true_separation_count": int(event_table.true_separation_after.sum()) if len(event_table) else 0,
        "contact_sector_order": "-".join(map(str, sector_bins)),
        "contact_sector_switches": sector_switches,
        "unique_contact_sectors": int(len(set(sector_bins))),
        "longest_contact_us": float(event_table.duration_us.max()) if len(event_table) else 0.0,
        "peak_contact_force_N": float(force_magnitude.max()),
        "both_end_support_fraction": float(both_support.mean()),
        "longest_20um_bridge_ms": float(longest_duration(ts, both_support.to_numpy()) * 1.0e3),
        "minimum_analytic_gap_um": float(gap.min_analytic_gap_um.min()),
        "minimum_end_margin_mm": float(gap.minimum_end_margin_mm.min()),
        "ALLKE_range": energy_ranges.get("ALLKE"), "ETOTAL_range": energy_ranges.get("ETOTAL"),
        "sample_count": len(sample), "dense_increment_count": len(time),
    }
    (OUT / "straight_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pd.DataFrame([metrics]).to_csv(OUT / "straight_metrics.csv", index=False)

    head = com - COM_X * axis
    tail = com + (2.4 - COM_X) * axis
    pose = pd.DataFrame({
        "time_s": ts, "cycle_fraction": ts / PERIOD, "s_mm": s, "delta_s_mm": s - s[0],
        "com_x": com[:, 0], "com_y": com[:, 1], "com_z": com[:, 2],
        "center_x": centers[:, 0], "center_y": centers[:, 1], "center_z": centers[:, 2],
        "head_x": head[:, 0], "head_y": head[:, 1], "head_z": head[:, 2],
        "tail_x": tail[:, 0], "tail_y": tail[:, 1], "tail_z": tail[:, 2],
        "axis_x": axis[:, 0], "axis_y": axis[:, 1], "axis_z": axis[:, 2],
        "t_x": tangent[0], "t_y": tangent[1], "t_z": tangent[2],
        "e1_x": e1[0], "e1_y": e1[1], "e1_z": e1[2],
        "e2_x": e2[0], "e2_y": e2[1], "e2_z": e2[2],
        "q1": q1, "q2": q2, "robot_phase_rad": robot_phase, "B_phase_rad": field_phase,
        "phase_lag_rad": phase_lag, "tilt_deg": tilt,
        "Bx_T": field[:, 0], "By_T": field[:, 1], "Bz_T": field[:, 2],
        "contact_active": contact_sample.astype(int),
    })
    pose.to_csv(OUT / "straight_pose.csv", index=False)

    curved = json.loads((BASE / "three_cycle_summary.json").read_text())["cycles"][0]
    comparison = {"curved_cycle_1": curved, "straight_cycle_1": metrics}
    (OUT / "curved_vs_straight_metrics.json").write_text(json.dumps(comparison, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
