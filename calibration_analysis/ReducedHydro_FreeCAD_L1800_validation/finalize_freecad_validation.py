"""Finalize the exact-CAD validation from existing extraction products only.

This script never opens an ODB and never repeats the whole-trajectory gap scan.
The 10 um bridge state is evaluated only for nodes already proven to be within
20 um in the saved cluster table.
"""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
OLD = HERE.parent / "ReducedHydro_geometry_L1800_validation"
AUDIT = HERE.parent / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall
sys.path.insert(0, str(AUDIT / "normal_contact_damping_probe"))
from contact_probe_common import dense_rp

JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083"
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
DT_COARSE = 1e-6


def interp_at(time, frame, column):
    return float(np.interp(time, frame.time_s, frame[column]))


def true_intervals(times, flags, sample_dt=DT_COARSE):
    times = np.asarray(times, float)
    flags = np.asarray(flags, bool)
    rows = []
    start = None
    previous = None
    for time, flag in zip(times, flags):
        if flag and start is None:
            start = time
        if start is not None and (not flag or (previous is not None and time - previous > 1.5 * sample_dt)):
            end = previous
            rows.append((start, end, end - start + sample_dt))
            start = time if flag else None
        previous = time
    if start is not None:
        rows.append((start, previous, previous - start + sample_dt))
    return rows


def crossing_events(axis, field, phase, bridge, contacts):
    dot = axis.axis_dot_tangent.to_numpy(float)
    time = axis.time_s.to_numpy(float)
    initial_sign = -1 if dot[0] < 0 else 1
    indices = np.where(np.signbit(dot[1:]) != np.signbit(dot[:-1]))[0] + 1
    rows = []
    for number, right in enumerate(indices, 1):
        left = right - 1
        fraction = -dot[left] / (dot[right] - dot[left])
        crossing = time[left] + fraction * (time[right] - time[left])
        before = -1 if dot[left] < 0 else 1
        after = -1 if dot[right] < 0 else 1
        active = bool(np.any((crossing >= contacts.start_s) & (crossing <= contacts.end_s)))
        rows.append({
            "event_number": number,
            "crossing_time_s": crossing,
            "crossing_time_ms": crossing * 1e3,
            "direction": "initial_to_reversed" if after != initial_sign else "reversed_to_initial",
            "dot_before": dot[left],
            "dot_after": dot[right],
            "directed_tilt_deg": 90.0,
            "field_tilt_deg": interp_at(crossing, field, "theta_B_deg"),
            "robot_phase_deg": np.degrees(interp_at(crossing, phase, "phi_robot_rad")),
            "robot_phase_advance_deg": interp_at(crossing, phase, "phi_robot_advance_deg"),
            "B_phase_deg": np.degrees(interp_at(crossing, phase, "phi_B_local_rad")),
            "contact_active": active,
            "bridge_20um": bool(round(interp_at(crossing, bridge, "opposing_bridge"))),
        })
    return pd.DataFrame(rows), initial_sign


def bridge_at_threshold(cluster, threshold_um, mesh, rp, dense, wall, axis, head_ids, tail_ids):
    base20 = cluster[cluster.threshold_um == 20.0].sort_values("increment").copy()
    if threshold_um in (20.0, 5.0, 0.0):
        source = cluster[cluster.threshold_um == threshold_um].sort_values("increment")
        return source[["time_s", "increment", "node_count", "cluster_count",
                       "sector_centers_deg", "max_sector_separation_deg", "head_involved", "tail_involved",
                       "participating_nodes", "opposing_bridge"]].copy()

    rotations = Rotation.from_rotvec(dense["UR"])
    node_map = mesh["nodes"]
    rows = []
    for row in base20.itertuples(index=False):
        ids = [int(value) for value in str(row.participating_nodes).split(";") if value and value != "nan"]
        if not ids:
            rows.append({"time_s": row.time_s, "increment": row.increment, "node_count": 0,
                         "cluster_count": 0, "max_sector_separation_deg": 0.0,
                         "sector_centers_deg": "",
                         "head_involved": False, "tail_involved": False,
                         "participating_nodes": "", "opposing_bridge": False})
            continue
        points0 = np.asarray([node_map[i] for i in ids]) + mesh["shift"]
        i = int(row.increment)
        points = rp + dense["U"][i] + rotations[i].apply(points0 - rp)
        gaps, triangles, _ = wall.query(points)
        keep = gaps * 1e3 <= threshold_um
        kept_ids = np.asarray(ids)[keep]
        if keep.any():
            tangent = axis.loc[i, ["tangent_x", "tangent_y", "tangent_z"]].to_numpy(float)
            e1 = axis.loc[i, ["e1_x", "e1_y", "e1_z"]].to_numpy(float)
            e2 = axis.loc[i, ["e2_x", "e2_y", "e2_z"]].to_numpy(float)
            normals = wall.normals[triangles[keep]]
            angles = np.mod(np.arctan2(normals @ e2, normals @ e1), 2 * np.pi)
            angles = np.sort(angles)
            gaps_a = np.diff(np.r_[angles, angles[0] + 2 * np.pi])
            cut = int(np.argmax(gaps_a))
            ordered = np.r_[angles[cut + 1:], angles[:cut + 1] + 2 * np.pi]
            groups = np.split(ordered, np.where(np.diff(ordered) > np.radians(45.0))[0] + 1)
            centers = [np.mod(np.angle(np.mean(np.exp(1j * group))), 2 * np.pi) for group in groups if len(group)]
            separation = max((np.degrees(np.arccos(np.clip(np.cos(a - b), -1, 1)))
                              for q, a in enumerate(centers) for b in centers[q + 1:]), default=0.0)
        else:
            centers, separation = [], 0.0
        rows.append({
            "time_s": row.time_s, "increment": row.increment, "node_count": len(kept_ids),
            "cluster_count": len(centers), "max_sector_separation_deg": separation,
            "sector_centers_deg": ";".join(f"{np.degrees(value):.3f}" for value in centers),
            "head_involved": bool(set(kept_ids) & head_ids),
            "tail_involved": bool(set(kept_ids) & tail_ids),
            "participating_nodes": ";".join(map(str, kept_ids)),
            "opposing_bridge": len(centers) >= 2 and separation >= 120.0,
        })
    return pd.DataFrame(rows)


def main():
    summary_path = HERE / "freecad_8p333_summary.json"
    summary = json.loads(summary_path.read_text())
    axis = pd.read_csv(HERE / "freecad_8p333_true_axis.csv")
    field = pd.read_csv(HERE / "freecad_8p333_field_orientation.csv")
    phase = pd.read_csv(HERE / "freecad_8p333_local_phase.csv")
    bridge = pd.read_csv(HERE / "freecad_8p333_bridge_timeline.csv")
    contacts = pd.read_csv(HERE / "freecad_8p333_contact_events.csv")
    clusters = pd.read_csv(HERE / "freecad_8p333_contact_clusters.csv")

    events, initial_sign = crossing_events(axis, field, phase, bridge, contacts)
    assert len(events) == 4
    events.to_csv(HERE / "freecad_head_tail_reversal_events.csv", index=False)
    summary["initial_axis_polarity"] = int(initial_sign)
    summary["initial_axis_dot_tangent"] = float(axis.axis_dot_tangent.iloc[0])
    summary["axis_zero_crossing_count"] = int(len(events))
    summary["first_axis_reversal_s"] = float(events.crossing_time_s.iloc[0])
    summary["classification"] = "EXACT_CAD_HEAD_SHAPE_REINTRODUCES_GEOMETRIC_JAM"
    summary["secondary_mechanism"] = "ROTATIONAL_OVERSHOOT_PERSISTS_AND_HEAD_TAIL_POLARITY_REVERSES"
    summary["contact_interpretation"] = "intermittent_microsecond_impacts_not_continuous_contact_jam"

    meshes, rp, _ = mesh_properties((HERE / f"{JOB}.inp").read_text())
    mesh = meshes["Robot_SOLID"]
    _, dense = dense_rp(HERE / "candidate_private")
    wall = Wall()
    head_ids, tail_ids = set(summary["head_node_ids"]), set(summary["tail_node_ids"])
    tables = {}
    sensitivity_rows = []
    for threshold in (20.0, 10.0, 5.0, 0.0):
        table = bridge_at_threshold(clusters, threshold, mesh, rp, dense, wall, axis, head_ids, tail_ids)
        tables[threshold] = table
        intervals = true_intervals(table.time_s, table.opposing_bridge)
        longest = max(intervals, key=lambda item: item[2]) if intervals else (np.nan, np.nan, 0.0)
        sensitivity_rows.append({
            "threshold_um": threshold, "bridge_event_count": len(intervals),
            "first_bridge_s": intervals[0][0] if intervals else np.nan,
            "longest_bridge_start_s": longest[0], "longest_bridge_end_s": longest[1],
            "longest_bridge_ms": longest[2] * 1e3,
            "persistent_over_0p5ms": longest[2] > 0.0005,
            "sampling_identity": "saved 1 us cluster poses; 10 um restricted to saved 20 um candidate nodes",
        })
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(HERE / "freecad_bridge_threshold_sensitivity.csv", index=False)

    longest20 = sensitivity.loc[sensitivity.threshold_um == 20.0].iloc[0]
    mask = ((tables[20.0].time_s >= longest20.longest_bridge_start_s) &
            (tables[20.0].time_s <= longest20.longest_bridge_end_s) &
            tables[20.0].opposing_bridge)
    segment = tables[20.0][mask]
    all_ids = set()
    for value in segment.participating_nodes.fillna(""):
        all_ids.update(int(item) for item in str(value).split(";") if item)
    gap = pd.read_csv(HERE / "freecad_8p333_exact_gap.csv")
    gap_segment = gap[(gap.time_s >= longest20.longest_bridge_start_s) &
                      (gap.time_s <= longest20.longest_bridge_end_s)]
    max_sep_row = segment.loc[segment.max_sector_separation_deg.idxmax()]
    composition = pd.DataFrame([{
        "onset_s": longest20.longest_bridge_start_s,
        "end_s": longest20.longest_bridge_end_s,
        "duration_ms": longest20.longest_bridge_ms,
        "minimum_gap_um": float(gap_segment.gap_um.min()),
        "unique_head_near_wall_nodes": len(all_ids & head_ids),
        "unique_tail_near_wall_nodes": len(all_ids & tail_ids),
        "unique_body_near_wall_nodes": len(all_ids - head_ids - tail_ids),
        "head_involved_fraction": float(segment.head_involved.mean()),
        "tail_involved_fraction": float(segment.tail_involved.mean()),
        "head_tail_simultaneous_fraction": float((segment.head_involved & segment.tail_involved).mean()),
        "wall_sector_A_deg": str(max_sep_row.get("sector_centers_deg", "")).split(";")[0],
        "wall_sector_B_deg": str(max_sep_row.get("sector_centers_deg", "")).split(";")[-1],
        "maximum_sector_separation_deg": float(segment.max_sector_separation_deg.max()),
        "interpretation": "HEAD_and_TAIL_opposing_near_wall_bridge",
    }])
    composition.to_csv(HERE / "freecad_longest_bridge_composition.csv", index=False)

    preflight = json.loads((HERE / "freecad_preflight_identity.json").read_text())
    old = json.loads((OLD / "L1800_8p333_summary.json").read_text())
    comparison = pd.DataFrame([
        {"case": "Old scaled L1800", "contact_events": old["contact_event_count"],
         "longest_contact_us": old["longest_contact_us"], "min_gap_um": old["min_exact_gap_um"],
         "longest_bridge_20um_ms": old["longest_bridge_ms"], "persistent_bridge_20um": old["persistent_bridge"],
         "magnetic_moment_Am2": old["candidate_magnetic_moment_Am2"],
         "Iperp_tonne_mm2": preflight["old_transverse_inertia_tonne_mm2"],
         "max_directed_tilt_deg": np.nan, "final_directed_tilt_deg": np.nan},
        {"case": "Exact FreeCAD L1800", "contact_events": summary["contact_event_count"],
         "longest_contact_us": summary["longest_contact_us"], "min_gap_um": summary["min_exact_gap_um"],
         "longest_bridge_20um_ms": summary["longest_bridge_ms"], "persistent_bridge_20um": summary["persistent_bridge"],
         "magnetic_moment_Am2": preflight["new_moment_Am2"],
         "Iperp_tonne_mm2": preflight["CAD_inertia_principal_tonne_mm2"][1],
         "max_directed_tilt_deg": summary["max_directed_tilt_deg"],
         "final_directed_tilt_deg": summary["final_directed_tilt_deg"]},
        {"case": "New/old ratio", "contact_events": np.nan, "longest_contact_us": np.nan,
         "min_gap_um": np.nan, "longest_bridge_20um_ms": summary["longest_bridge_ms"] / old["longest_bridge_ms"],
         "persistent_bridge_20um": np.nan, "magnetic_moment_Am2": preflight["moment_ratio"],
         "Iperp_tonne_mm2": preflight["new_old_transverse_inertia_ratio"],
         "max_directed_tilt_deg": preflight["new_old_moment_per_Iperp_ratio"],
         "final_directed_tilt_deg": np.nan},
    ])
    comparison.to_csv(HERE / "old_scaled_vs_freecad_summary.csv", index=False)

    torque = pd.read_csv(HERE / "freecad_8p333_magnetic_hydro_torque.csv")
    torque.to_csv(HERE / "freecad_8p333_torque.csv", index=False)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    runtime = pd.DataFrame([summary]).drop(columns=["head_node_ids", "tail_node_ids"])
    runtime.to_csv(HERE / "freecad_8p333_runtime_identity.csv", index=False)
    print(events[["crossing_time_ms", "direction"]].to_string(index=False))
    print(sensitivity.to_string(index=False))
    print(composition.to_string(index=False))


if __name__ == "__main__":
    main()
