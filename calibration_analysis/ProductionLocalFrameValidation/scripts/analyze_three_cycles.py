"""Cycle-resolved analysis of the unique production local-frame run."""
from difflib import SequenceMatcher
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from continuous_segment_reference import ContinuousSegmentReference, unit
from regress_production_local_frame import DXF, TRANSFORM, RP0, VENDORED, build, load_module


HERE = Path(__file__).resolve().parent
VALIDATION = HERE.parent
REPO = VALIDATION.parents[1]
CASE = VALIDATION / "case" / "PROD_LOCAL30_G0_3CYCLE"
PRIVATE = CASE / "private"
JOB = "PROD_LOCAL30_G0_3CYCLE"
PERIOD = 1.0 / 30.0
A0 = unit(np.array([0.9647382600216, -0.1188742372140, 0.2348382536499]))
sys.path[:0] = [str(REPO / "calibration_analysis" / "MagneticDriveFrameAudit" / "scripts")]
from analyze_diagnostics import dense, contact_force, intervals, merge_ranges, coarse_proximity


class TubeAdapter:
    def __init__(self, reference):
        self.reference = reference

    def at(self, position_aba):
        point_mag = self.reference.position_to_mag(position_aba)
        s, center_mag = self.reference.project_mag(point_mag)
        tangent, e1, e2 = self.reference.frame_at(s)
        center = self.reference.origin_aba + self.reference.R.T.dot(center_mag)
        return (s, center, self.reference.vector_to_aba(tangent),
                self.reference.vector_to_aba(e1), self.reference.vector_to_aba(e2))


def longest_duration(time, flag):
    return max((time[end] - time[start] + np.median(np.diff(time))
                for start, end in intervals(flag)), default=0.0)


def polygon_area(x, y):
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def cycle_slice(time, cycle):
    lower, upper = (cycle - 1) * PERIOD, cycle * PERIOD
    return (time >= lower - 1.0e-12) & (time <= upper + 1.0e-12)


def main():
    production = load_module("three_cycle_production", VENDORED)
    model, _, _ = build(production, field_mode="ROBOT_LOCAL_TANGENT")
    reference = ContinuousSegmentReference(DXF, TRANSFORM, model.robot_axis_global,
                                            polarity=model.robot_polarity)
    tube = TubeAdapter(reference)
    time_dense, dense_data = dense(CASE)
    contact = contact_force(CASE, time_dense)
    contact_magnitude = np.linalg.norm(contact, axis=1)
    active_dense = contact_magnitude > 1.0e-8
    dt = float(np.median(np.diff(time_dense)))
    events_dense = merge_ranges(intervals(active_dense), int(round(5.0e-6 / dt)))

    sample_index = np.unique(np.r_[np.arange(0, len(time_dense), 100), len(time_dense) - 1])
    time = time_dense[sample_index]
    displacement = dense_data["U"][sample_index]
    ur = dense_data["UR"][sample_index]
    com = RP0 + displacement
    rotations = Rotation.from_rotvec(ur)
    axis = rotations.apply(np.broadcast_to(A0, (len(time), 3)))
    reference.reset()
    frame = [tube.at(position) for position in com]
    s = np.asarray([row[0] for row in frame])
    centers = np.asarray([row[1] for row in frame])
    tangent = np.asarray([row[2] for row in frame])
    e1 = np.asarray([row[3] for row in frame])
    e2 = np.asarray([row[4] for row in frame])
    q1 = np.einsum("ij,ij->i", axis, e1)
    q2 = np.einsum("ij,ij->i", axis, e2)
    axial = np.einsum("ij,ij->i", axis, tangent)
    tilt = np.degrees(np.arccos(np.clip(axial, -1.0, 1.0)))
    robot_phase = np.unwrap(np.arctan2(q2, q1))

    telemetry = pd.read_csv(CASE / f"{JOB}_telemetry.csv")
    field = np.column_stack([np.interp(time, telemetry.t_s, telemetry[column])
                             for column in ("Bx_aba_T", "By_aba_T", "Bz_aba_T")])
    b1 = np.einsum("ij,ij->i", field, e1)
    b2 = np.einsum("ij,ij->i", field, e2)
    field_phase = np.unwrap(np.arctan2(b2, b1))
    phase_lag = np.angle(np.exp(1j * (robot_phase - field_phase)))
    phase_valid = np.sin(np.radians(tilt)) > 0.1

    proximity_index = np.unique(np.r_[np.arange(0, len(time_dense), 1000), len(time_dense) - 1])
    reference.reset()
    proximity = coarse_proximity({}, RP0 + dense_data["U"],
                                 Rotation.from_rotvec(dense_data["UR"]), tube, proximity_index)
    proximity["time_s"] = time_dense[proximity.sample_index.to_numpy(int)]
    proximity.to_csv(VALIDATION / "exact_gap_timeseries.csv", index=False)
    both_support = (proximity.HEAD_gap_um <= 20.0) & (proximity.TAIL_gap_um <= 20.0)

    energy_archive = np.load(PRIVATE / "energy_history_private.npz")
    energy = {key: energy_archive[key] for key in energy_archive.files}
    event_rows = []
    event_sequences = {1: [], 2: [], 3: []}
    separation_steps = int(round(20.0e-6 / dt))
    reference.reset()
    for number, (start, end) in enumerate(events_dense, 1):
        mid = (start + end) // 2
        direction = -np.mean(contact[start:end + 1], axis=0)
        _, _, _, event_e1, event_e2 = tube.at(RP0 + dense_data["U"][mid])
        sector = math.degrees(math.atan2(np.dot(direction, event_e2), np.dot(direction, event_e1))) % 360.0
        cycle = min(3, int(time_dense[mid] / PERIOD) + 1)
        separated = bool(end + separation_steps < len(active_dense) and
                         not active_dense[end + 1:end + 1 + separation_steps].any())
        event_sequences[cycle].append(int(round(sector / 30.0)) % 12)
        event_rows.append({"event": number, "cycle": cycle, "start_s": time_dense[start],
                           "end_s": time_dense[end], "duration_us": (end - start + 1) * dt * 1e6,
                           "peak_force_N": float(contact_magnitude[start:end + 1].max()),
                           "sector_deg": sector, "true_separation_after": int(separated)})
    events = pd.DataFrame(event_rows)
    events.to_csv(VALIDATION / "contact_events.csv", index=False)

    cycle_rows = []
    for cycle in (1, 2, 3):
        mask = cycle_slice(time, cycle)
        local_events = events[events.cycle == cycle] if len(events) else events
        prox_mask = cycle_slice(proximity.time_s.to_numpy(), cycle)
        q1c, q2c = q1[mask], q2[mask]
        centroid = np.array([q1c.mean(), q2c.mean()])
        radial = np.sqrt((q1c - centroid[0]) ** 2 + (q2c - centroid[1]) ** 2)
        lag_mask = mask & phase_valid
        mean_lag = math.degrees(math.atan2(np.mean(np.sin(phase_lag[lag_mask])),
                                           np.mean(np.cos(phase_lag[lag_mask]))))
        energy_ranges = {}
        for name, values in energy.items():
            emask = cycle_slice(values[:, 0], cycle)
            if emask.any():
                energy_ranges[name] = float(np.ptp(values[emask, 1]))
        cycle_rows.append({
            "cycle": cycle,
            "local_B_winding": float((field_phase[mask][-1] - field_phase[mask][0]) / (2 * math.pi)),
            "robot_axis_winding": float((robot_phase[mask][-1] - robot_phase[mask][0]) / (2 * math.pi)),
            "mean_phase_lag_deg": mean_lag, "phase_valid_fraction": float(phase_valid[mask].mean()),
            "tilt_min_deg": float(tilt[mask].min()), "tilt_mean_deg": float(tilt[mask].mean()),
            "tilt_max_deg": float(tilt[mask].max()), "orbit_area": polygon_area(q1c, q2c),
            "orbit_centroid_e1": centroid[0], "orbit_centroid_e2": centroid[1],
            "radial_RMS": float(np.sqrt(np.mean(radial ** 2))),
            "contact_event_count": int(len(local_events)),
            "true_separation_count": int(local_events.true_separation_after.sum()) if len(local_events) else 0,
            "contact_sector_order": "-".join(map(str, event_sequences[cycle])),
            "both_end_support_fraction": float(both_support[prox_mask].mean()),
            "longest_contact_us": float(local_events.duration_us.max()) if len(local_events) else 0.0,
            "longest_bridge_ms": longest_duration(proximity.time_s.to_numpy()[prox_mask], both_support.to_numpy()[prox_mask]) * 1e3,
            "delta_s_mm": float(s[mask][-1] - s[mask][0]),
            "ALLKE_range": energy_ranges.get("ALLKE", float("nan")),
            "ETOTAL_range": energy_ranges.get("ETOTAL", float("nan")),
        })
    cycles = pd.DataFrame(cycle_rows)
    cycles.to_csv(VALIDATION / "three_cycle_metrics.csv", index=False)
    c2, c3 = cycles.iloc[1], cycles.iloc[2]
    sequence_similarity = SequenceMatcher(None, event_sequences[2], event_sequences[3]).ratio()
    comparison = {
        "local_B_winding_difference": float(abs(c2.local_B_winding - c3.local_B_winding)),
        "winding_difference": float(abs(c2.robot_axis_winding - c3.robot_axis_winding)),
        "mean_phase_lag_difference_deg": float(abs(np.degrees(np.angle(np.exp(1j * np.radians(c2.mean_phase_lag_deg - c3.mean_phase_lag_deg)))))),
        "phase_valid_fraction_difference": float(abs(c2.phase_valid_fraction - c3.phase_valid_fraction)),
        "tilt_min_difference_deg": float(abs(c2.tilt_min_deg - c3.tilt_min_deg)),
        "tilt_mean_difference_deg": float(abs(c2.tilt_mean_deg - c3.tilt_mean_deg)),
        "tilt_max_difference_deg": float(abs(c2.tilt_max_deg - c3.tilt_max_deg)),
        "orbit_centroid_distance": float(np.hypot(c2.orbit_centroid_e1 - c3.orbit_centroid_e1,
                                                   c2.orbit_centroid_e2 - c3.orbit_centroid_e2)),
        "orbit_area_difference": float(abs(c2.orbit_area - c3.orbit_area)),
        "radial_RMS_difference": float(abs(c2.radial_RMS - c3.radial_RMS)),
        "tilt_envelope_difference_deg": float(max(abs(c2.tilt_min_deg - c3.tilt_min_deg),
                                                   abs(c2.tilt_max_deg - c3.tilt_max_deg))),
        "contact_event_count_difference": int(abs(c2.contact_event_count - c3.contact_event_count)),
        "true_separation_count_difference": int(abs(c2.true_separation_count - c3.true_separation_count)),
        "contact_sector_sequence_similarity": float(sequence_similarity),
        "both_end_support_fraction_difference": float(abs(c2.both_end_support_fraction - c3.both_end_support_fraction)),
        "longest_contact_difference_us": float(abs(c2.longest_contact_us - c3.longest_contact_us)),
        "longest_bridge_difference_ms": float(abs(c2.longest_bridge_ms - c3.longest_bridge_ms)),
        "delta_s_difference_mm": float(abs(c2.delta_s_mm - c3.delta_s_mm)),
        "ALLKE_range_difference": float(abs(c2.ALLKE_range - c3.ALLKE_range)),
        "ETOTAL_range_difference": float(abs(c2.ETOTAL_range - c3.ETOTAL_range)),
    }
    pd.DataFrame([comparison]).to_csv(VALIDATION / "cycle2_vs_cycle3.csv", index=False)
    pose = pd.DataFrame({"time_s": time, "cycle": np.minimum(3, (time / PERIOD).astype(int) + 1),
        "s_mm": s, "delta_s_mm": s - s[0], "com_x": com[:, 0], "com_y": com[:, 1], "com_z": com[:, 2],
        "center_x": centers[:, 0], "center_y": centers[:, 1], "center_z": centers[:, 2],
        "axis_x": axis[:, 0], "axis_y": axis[:, 1], "axis_z": axis[:, 2],
        "t_x": tangent[:, 0], "t_y": tangent[:, 1], "t_z": tangent[:, 2],
        "e1_x": e1[:, 0], "e1_y": e1[:, 1], "e1_z": e1[:, 2], "e2_x": e2[:, 0], "e2_y": e2[:, 1], "e2_z": e2[:, 2],
        "q1": q1, "q2": q2, "robot_phase_rad": robot_phase, "B_phase_rad": field_phase,
        "phase_lag_rad": phase_lag, "phase_valid": phase_valid.astype(int), "tilt_deg": tilt,
        "Bx_T": field[:, 0], "By_T": field[:, 1], "Bz_T": field[:, 2],
        "contact_active": active_dense[sample_index].astype(int)})
    pose.to_csv(VALIDATION / "three_cycle_pose.csv", index=False)
    gate_results = {
        "full_bounded_rotation": bool(all(0.8 <= abs(row["robot_axis_winding"]) <= 1.2 and
                                           row["tilt_max_deg"] < 90.0 for row in cycle_rows[1:])),
        "orbit_geometry_repeatable": bool(comparison["orbit_centroid_distance"] < 0.05 and
                                            comparison["radial_RMS_difference"] < 0.05 and
                                            comparison["tilt_envelope_difference_deg"] < 5.0 and
                                            comparison["winding_difference"] < 0.1),
        "contact_sequence_repeatable": bool(sequence_similarity >= 0.8 and
                                              comparison["contact_event_count_difference"] <=
                                              max(2, 0.25 * max(c2.contact_event_count, c3.contact_event_count)) and
                                              comparison["longest_contact_difference_us"] <= 50.0),
        "axial_response_repeatable": bool(comparison["delta_s_difference_mm"] <= 0.25),
        "repeated_hit_and_separation": bool(c2.true_separation_count > 0 and c3.true_separation_count > 0),
    }
    stable = all(gate_results.values())
    classification = ("LOCAL_FRAME_STABLE_WALL_WOBBLE" if stable else
                      "LOCAL_FRAME_ROTATIONAL_WOBBLE_NOT_LIMIT_CYCLE" if
                      gate_results["full_bounded_rotation"] else "LOCAL_FRAME_NO_STABLE_WOBBLE")
    summary = {"classification": classification, "cycles": cycle_rows,
               "cycle2_vs_cycle3": comparison, "contact_events": len(events_dense),
               "minimum_exact_gap_um": float(proximity.min_gap_um.min()),
               "energy_variables": sorted(energy), "gate_results": gate_results,
               "gate_thresholds": {"sequence_similarity_min": 0.8,
                   "relative_contact_count_difference_max": 0.25,
                   "longest_contact_difference_us_max": 50.0,
                   "delta_s_difference_mm_max": 0.25,
                   "orbit_centroid_distance_max": 0.05,
                   "radial_RMS_difference_max": 0.05,
                   "tilt_envelope_difference_deg_max": 5.0,
                   "robot_winding_difference_max": 0.1}}
    (VALIDATION / "three_cycle_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(cycles.to_string(index=False))
    print(json.dumps(comparison, indent=2))
    print(classification, json.dumps(gate_results, indent=2))


if __name__ == "__main__":
    main()
