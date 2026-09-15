"""Full-cycle diagnostics centered on robot-axis winding, not displacement."""
from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from audit_common import AUDIT, REPO, RP0, A0, TransportedTubeFrame, unit


HIDDEN = REPO / "calibration_analysis/ReducedHydro_hidden_impact_audit"
VALIDATION = REPO / "calibration_analysis/ReducedHydro_L2300_frozenCAD_validation"
import sys
sys.path[:0] = [str(VALIDATION), str(HIDDEN)]
from audit_frozen_cad_static import Wall, exact_query


CASES = ("FRAME_CURRENT", "FRAME_LOCAL30", "FRAME_LOCAL40", "FRAME_LOCAL30_RAMP", "FRAME_LOCAL_BEST_CONTACT")


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    return max((archive[key] for key in keys), key=len)


def dense(folder):
    archive = np.load(folder / "private/rp_history_private.npz")
    base = get_array(archive, "V1"); time = base[:, 0]
    data = {prefix: np.column_stack([np.interp(time, get_array(archive, prefix + str(i))[:, 0], get_array(archive, prefix + str(i))[:, 1])
                                     for i in (1, 2, 3)]) for prefix in ("U", "UR", "V", "VR")}
    return time, data


def contact_force(folder, time):
    archive = np.load(folder / "private/contact_history_private.npz")
    def vector(prefix):
        columns = []
        for i in (1, 2, 3):
            keys = [key for key in archive.files if "|" + prefix + str(i) + " on surface " in key and "ASSEMBLY_ROBOT" in key]
            array = max((archive[key] for key in keys), key=len)
            columns.append(np.interp(time, array[:, 0], array[:, 1]))
        return np.column_stack(columns)
    return vector("CFN") + vector("CFS")


def intervals(flag):
    edge = np.diff(np.r_[False, flag, False].astype(int))
    return list(zip(np.where(edge == 1)[0], np.where(edge == -1)[0] - 1))


def merge_ranges(ranges, max_gap):
    merged = []
    for start, end in ranges:
        if merged and start - merged[-1][1] - 1 <= max_gap:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def rotation_x_to_axis(axis):
    x = np.array([1.0, 0.0, 0.0]); v = np.cross(x, axis); c = float(np.dot(x, axis))
    k = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + k + k.dot(k) / (1.0 + c)


def longest_boolean_duration(flag, sample_dt):
    return max(((end - start + 1) * sample_dt for start, end in intervals(flag)), default=0.0)


def coarse_proximity(item, com, rotations, tube, indices):
    mesh = REPO / "calibration_analysis/CoarseMotionModeScreen/cases/_meshes/Robot_SCREENING_L2400_D0815_nodes.csv"
    props = json.loads((REPO / "cad/freecad_parametric_robot/screening/coarse_motion_mode/Robot_SCREENING_L2400_D0815_geometry.json").read_text())
    nodes = pd.read_csv(mesh); local = nodes[["x_mm", "y_mm", "z_mm"]].to_numpy(float); x = local[:, 0]
    mask = (x <= .26149478) | (x >= 2.4 - 1e-7) | (np.abs(np.linalg.norm(local[:, 1:], axis=1) - .4075) < 2e-4)
    local = local[mask][::2]; x = x[mask][::2]
    region = np.where(x <= .26149478, "HEAD", np.where(x >= 2.4 - 1e-7, "TAIL", "BODY"))
    initial = RP0 + (local - np.asarray(props["center_of_mass_mm"])).dot(rotation_x_to_axis(A0).T)
    wall = Wall(); rows = []
    for i in indices:
        _, _, c, e1, e2 = tube.at(com[i])
        points = com[i] + rotations[i].apply(initial - RP0)
        gaps, triangle, _ = exact_query(wall, points)
        angles = np.degrees(np.arctan2(wall.normals[triangle].dot(e2), wall.normals[triangle].dot(e1))) % 360.0
        nearest = int(np.argmin(gaps))
        rows.append({"sample_index": int(i), "min_gap_um": float(gaps.min() * 1e3),
                     "HEAD_gap_um": float(gaps[region == "HEAD"].min() * 1e3),
                     "TAIL_gap_um": float(gaps[region == "TAIL"].min() * 1e3),
                     "nearest_wall_sector_deg": float(angles[nearest])})
    return pd.DataFrame(rows)


def analyze(case_id):
    folder = AUDIT / "cases" / case_id
    item = json.loads((folder / "case_identity.json").read_text())
    if item.get("status") != "SOLVED":
        return None
    time, data = dense(folder); dt = float(np.median(np.diff(time)))
    sample = np.unique(np.r_[np.arange(0, len(time), 100), len(time) - 1])
    tube = TransportedTubeFrame(); com = RP0 + data["U"]
    rotations_all = Rotation.from_rotvec(data["UR"]); rotations = rotations_all[sample]
    axis = rotations.apply(np.broadcast_to(A0, (len(sample), 3)))
    frame_rows = [tube.at(position) for position in com[sample]]
    arc = np.array([row[0] for row in frame_rows]); centers = np.array([row[1] for row in frame_rows])
    c = np.array([row[2] for row in frame_rows]); e1 = np.array([row[3] for row in frame_rows]); e2 = np.array([row[4] for row in frame_rows])
    q1 = np.einsum("ij,ij->i", axis, e1); q2 = np.einsum("ij,ij->i", axis, e2)
    axial = np.einsum("ij,ij->i", axis, c)
    robot_phase = np.unwrap(np.arctan2(q2, q1))
    tilt = np.degrees(np.arccos(np.clip(axial, -1.0, 1.0)))
    telemetry = pd.read_csv(next(folder.glob("*_telemetry.csv")))
    field = np.column_stack([np.interp(time[sample], telemetry.t_s, telemetry[name]) for name in ("Bx_aba_T", "By_aba_T", "Bz_aba_T")])
    b1 = np.einsum("ij,ij->i", field, e1); b2 = np.einsum("ij,ij->i", field, e2); bt = np.einsum("ij,ij->i", field, c)
    field_phase = np.unwrap(np.arctan2(b2, b1)); command_phase = math.radians(248.0) + 2.0 * math.pi * 30.0 * time[sample]
    force = contact_force(folder, time); active = np.linalg.norm(force, axis=1) > 1e-8
    events = merge_ranges(intervals(active), int(round(5e-6 / dt)))
    separation_steps = int(round(20e-6 / dt))
    true_separations = sum(end + separation_steps < len(active) and not active[end+1:end+1+separation_steps].any() for _, end in events)
    proximity_indices = sample[::2]
    proximity = coarse_proximity(item, com, rotations_all, tube, proximity_indices)
    proximity["time_s"] = time[proximity.sample_index.to_numpy(int)]
    proximity.to_csv(folder / "proximity_diagnostic.csv", index=False)
    both = (proximity.HEAD_gap_um <= 20.0) & (proximity.TAIL_gap_um <= 20.0)
    prox_dt = float(np.median(np.diff(proximity.time_s)))
    proximity_sector = np.unwrap(np.radians(proximity.nearest_wall_sector_deg.to_numpy()))
    event_sectors = []
    for start, end in events:
        mid = (start + end) // 2
        wall_direction = -np.mean(force[start:end + 1], axis=0)
        if np.linalg.norm(wall_direction) <= 1e-15:
            continue
        _, _, _, event_e1, event_e2 = tube.at(com[mid])
        event_sectors.append(math.atan2(np.dot(wall_direction, event_e2), np.dot(wall_direction, event_e1)))
    event_sectors = np.unwrap(np.asarray(event_sectors))
    sector_transitions = int(np.sum(np.abs(np.diff(event_sectors)) > math.radians(30.0)))
    half_sample_window = max(1, int(round(0.5 * np.median(np.diff(sample)))))
    contact_sample = np.array([active[max(0, i-half_sample_window):min(len(active), i+half_sample_window+1)].any()
                               for i in sample])
    phase_difference = robot_phase - field_phase
    mean_phase_lag = math.atan2(np.mean(np.sin(phase_difference)), np.mean(np.cos(phase_difference)))
    contact_magnitude = np.linalg.norm(force, axis=1)
    peak_contact_index = int(np.argmax(contact_magnitude))
    first_contact_index = int(np.flatnonzero(active)[0]) if active.any() else None
    metrics = dict(item)
    metrics.update(
        analysis_status="VALID", total_robot_phase_winding=float((robot_phase[-1] - robot_phase[0]) / (2*math.pi)),
        axis_winding_number=float((robot_phase[-1] - robot_phase[0]) / (2*math.pi)),
        command_phase_winding=float((command_phase[-1] - command_phase[0]) / (2*math.pi)),
        local_B_phase_winding=float((field_phase[-1] - field_phase[0]) / (2*math.pi)),
        mean_phase_lag_deg=float(np.degrees(mean_phase_lag)), max_directed_tilt_deg=float(tilt.max()),
        final_directed_tilt_deg=float(tilt[-1]), tilt_range_deg=float(np.ptp(tilt)),
        contact_side_transitions=sector_transitions, true_separations=int(true_separations),
        both_end_support_fraction=float(both.mean()), longest_continuous_contact_us=max(((b-a+1)*dt*1e6 for a,b in events), default=0.0),
        longest_opposing_bridge_ms=float(longest_boolean_duration(both.to_numpy(), prox_dt)*1e3),
        delta_s_mm=float(arc[-1]-arc[0]), contact_event_count=len(events),
        first_contact_time_ms=None if first_contact_index is None else float(time[first_contact_index] * 1e3),
        peak_contact_force_N=float(contact_magnitude[peak_contact_index]),
        peak_contact_time_ms=float(time[peak_contact_index] * 1e3),
        wall_sector_span_deg=float(np.degrees(np.ptp(event_sectors))) if len(event_sectors) else 0.0,
        nearest_wall_sector_span_deg=float(np.degrees(np.ptp(proximity_sector))),
        min_coarse_gap_um=float(proximity.min_gap_um.min()))
    (folder / "diagnostic_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    com_x = 1.25606694407031
    head = com[sample] - com_x * axis; tail = com[sample] + (2.4 - com_x) * axis
    pose = pd.DataFrame({"time_s": time[sample], "cycle_fraction": time[sample]*30.0,
        "rp_x_mm": com[sample,0], "rp_y_mm": com[sample,1], "rp_z_mm": com[sample,2],
        "center_x_mm": centers[:,0], "center_y_mm": centers[:,1], "center_z_mm": centers[:,2],
        "head_x_mm": head[:,0], "head_y_mm": head[:,1], "head_z_mm": head[:,2],
        "tail_x_mm": tail[:,0], "tail_y_mm": tail[:,1], "tail_z_mm": tail[:,2],
        "axis_x": axis[:,0], "axis_y": axis[:,1], "axis_z": axis[:,2],
        "c_x": c[:,0], "c_y": c[:,1], "c_z": c[:,2], "e1_x": e1[:,0], "e1_y": e1[:,1], "e1_z": e1[:,2],
        "e2_x": e2[:,0], "e2_y": e2[:,1], "e2_z": e2[:,2], "q1": q1, "q2": q2,
        "robot_phase_unwrapped_rad": robot_phase, "field_phase_unwrapped_rad": field_phase,
        "command_phase_rad": command_phase, "directed_tilt_deg": tilt,
        "Bx_T": field[:,0], "By_T": field[:,1], "Bz_T": field[:,2], "B_t_T": bt, "B_e1_T": b1, "B_e2_T": b2,
        "contact_active": contact_sample.astype(int), "delta_s_mm": arc-arc[0]})
    pose.to_csv(folder / "pose_diagnostic.csv", index=False)
    return metrics


def main():
    results = []
    for case_id in CASES:
        if (AUDIT / "cases" / case_id / "case_identity.json").exists():
            result = analyze(case_id)
            if result is not None: results.append(result)
    if results:
        pd.DataFrame(results).to_csv(AUDIT / "diagnostic_metrics.csv", index=False)
        print(pd.DataFrame(results)[["case_id", "axis_winding_number", "local_B_phase_winding", "max_directed_tilt_deg", "true_separations", "contact_side_transitions", "delta_s_mm"]].to_string(index=False))


if __name__ == "__main__":
    main()
