"""Analyze prescribed replay reaction torque and geometric support."""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
JOB = "STRAIGHT_PRESCRIBED_WOBBLE_CONTACT_AUDIT"
CASE = HERE / "case" / JOB
BASE = REPO / "calibration_analysis/StraightPipeControl/case/PROD_LOCAL30_G0_STRAIGHT_CTRL"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
A0 /= np.linalg.norm(A0)


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    return max((archive[key] for key in keys), key=len)


def vector(archive, prefix, time=None, surface=None):
    columns = []
    for i in (1, 2, 3):
        if surface:
            candidates = [key for key in archive.files
                          if "|" + prefix + str(i) + " on surface " in key and surface in key]
            array = max((archive[key] for key in candidates), key=len)
        else:
            array = get_array(archive, prefix + str(i))
        if time is None:
            time = array[:, 0]
        columns.append(np.interp(time, array[:, 0], array[:, 1]))
    return time, np.column_stack(columns)


def surface_nodes():
    deck = (CASE / f"{JOB}.inp").read_text()
    part = re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    nodes = np.asarray([[float(value) for value in line.split(",")[1:4]]
                        for line in block.splitlines() if line.strip()])
    axial = (nodes - RP0) @ A0
    region = np.full(len(nodes), "BODY", object)
    region[axial <= axial.min() + 0.2615] = "HEAD"
    region[axial >= axial.max() - 0.02] = "TAIL"
    return nodes, region


def main():
    rp = np.load(CASE / "private/rp_history_private.npz")
    t, u = vector(rp, "U")
    _, ur = vector(rp, "UR", t)
    _, velocity = vector(rp, "V", t)
    _, omega = vector(rp, "VR", t)
    reaction = np.load(CASE / "private/reaction_history_private.npz")
    _, rm = vector(reaction, "RM", t)
    contact = np.load(CASE / "private/contact_history_private.npz")
    _, fn = vector(contact, "CFN", t, "ASSEMBLY_ROBOT")
    _, fs = vector(contact, "CFS", t, "ASSEMBLY_ROBOT")
    fcontact = fn + fs

    identity = json.loads((CASE / "case_identity.json").read_text())
    source_identity = json.loads((BASE / "case_identity.json").read_text())
    tangent = np.asarray(source_identity["straight_tangent_aba"], float)
    e1 = np.asarray(source_identity["initial_e1_aba"], float)
    e2 = np.asarray(source_identity["initial_e2_aba"], float)
    rotations = Rotation.from_rotvec(ur)
    axis = rotations.apply(np.broadcast_to(A0, omega.shape))
    tilt = np.degrees(np.arccos(np.clip(axis @ tangent, -1.0, 1.0)))
    spin = np.einsum("ij,ij->i", omega, axis)
    wobble = np.linalg.norm(omega - spin[:, None] * axis, axis=1)

    initial_nodes, region = surface_nodes()
    sample = np.unique(np.r_[np.arange(0, len(t), 100), len(t) - 1])
    rows = []
    radius = float(source_identity["wall_radius_mm"])
    center0 = np.asarray(source_identity["initial_center_aba_mm"], float)
    for index in sample:
        com = RP0 + u[index]
        points = com + rotations[index].apply(initial_nodes - RP0)
        axial = (points - center0) @ tangent
        projected = center0 + axial[:, None] * tangent
        gaps = radius - np.linalg.norm(points - projected, axis=1)
        rows.append({"time_s": t[index], "sample_index": int(index),
                     "min_gap_um": gaps.min() * 1e3,
                     "HEAD_gap_um": gaps[region == "HEAD"].min() * 1e3,
                     "TAIL_gap_um": gaps[region == "TAIL"].min() * 1e3})
    gap = pd.DataFrame(rows)
    gap["opposing_bridge"] = ((gap.HEAD_gap_um <= 20) & (gap.TAIL_gap_um <= 20)).astype(int)
    gap.to_csv(HERE / "prescribed_exact_gap_timeseries.csv", index=False)
    head_gap = np.interp(t, gap.time_s, gap.HEAD_gap_um)
    tail_gap = np.interp(t, gap.time_s, gap.TAIL_gap_um)
    bridge = (head_gap <= 20) & (tail_gap <= 20)

    baseline_tel = pd.read_csv(BASE / "PROD_LOCAL30_G0_STRAIGHT_CTRL_telemetry.csv")
    tmag = np.column_stack([np.interp(t, baseline_tel.t_s, baseline_tel[name])
                            for name in ("tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm")])
    local = lambda values: np.column_stack((values @ tangent, values @ e1, values @ e2))
    rm_l = local(rm); tmag_l = local(tmag)
    required_transverse = np.linalg.norm(rm_l[:, 1:], axis=1)
    available_transverse = np.linalg.norm(tmag_l[:, 1:], axis=1)
    active = np.linalg.norm(fcontact, axis=1) > 1e-8
    analysis = (t >= 0.0023) & (t <= 0.016667)
    bridge_analysis = analysis & bridge

    output = pd.DataFrame({
        "time_s": t, "RM_t_Nmm": rm_l[:, 0], "RM_e1_Nmm": rm_l[:, 1],
        "RM_e2_Nmm": rm_l[:, 2], "RM_norm_Nmm": np.linalg.norm(rm, axis=1),
        "Tmag_baseline_t_Nmm": tmag_l[:, 0], "Tmag_baseline_e1_Nmm": tmag_l[:, 1],
        "Tmag_baseline_e2_Nmm": tmag_l[:, 2], "Tmag_baseline_norm_Nmm": np.linalg.norm(tmag, axis=1),
        "required_transverse_Nmm": required_transverse,
        "available_transverse_Nmm": available_transverse,
        "omega_spin_rad_s": spin, "omega_wobble_rad_s": wobble, "tilt_deg": tilt,
        "HEAD_gap_um": head_gap, "TAIL_gap_um": tail_gap,
        "contact_active": active.astype(int), "opposing_bridge": bridge.astype(int),
        "com_x_mm": RP0[0] + u[:, 0], "com_y_mm": RP0[1] + u[:, 1], "com_z_mm": RP0[2] + u[:, 2],
        "axis_x": axis[:, 0], "axis_y": axis[:, 1], "axis_z": axis[:, 2],
    })
    source_sample = np.unique(np.r_[np.arange(0, len(output), 100), len(output) - 1])
    output.iloc[source_sample].to_csv(HERE / "prescribed_reaction_vs_magnetic.csv", index=False)

    ratio = required_transverse[bridge_analysis] / np.maximum(available_transverse[bridge_analysis], 1e-12)
    summary = {
        "case_id": JOB, "status": identity["status"],
        "primary_classification": "PRESCRIBED_MOTION_MASKED_GEOMETRIC_BRIDGE",
        "secondary_classification": "GEOMETRIC_INCOMPATIBILITY_DRIVES_PENALTY_REACTION",
        "analysis_window_ms": [2.3, 16.667],
        "prescribed_tilt_mean_deg": float(tilt[analysis].mean()),
        "prescribed_tilt_min_deg": float(tilt[analysis].min()),
        "prescribed_tilt_max_deg": float(tilt[analysis].max()),
        "both_end_support_fraction": float(bridge[analysis].mean()),
        "contact_active_fraction": float(active[analysis].mean()),
        "minimum_gap_um": float(np.minimum(head_gap[analysis], tail_gap[analysis]).min()),
        "required_transverse_RM_peak_Nmm": float(required_transverse[bridge_analysis].max()),
        "required_transverse_RM_RMS_Nmm": float(np.sqrt(np.mean(required_transverse[bridge_analysis] ** 2))),
        "available_magnetic_transverse_RMS_Nmm": float(np.sqrt(np.mean(available_transverse[bridge_analysis] ** 2))),
        "peak_required_to_available_ratio": float(ratio.max()),
        "RMS_required_to_available_ratio": float(
            np.sqrt(np.mean(required_transverse[bridge_analysis] ** 2)) /
            np.sqrt(np.mean(available_transverse[bridge_analysis] ** 2))),
        "bridge_time_required_exceeds_available_fraction": float((ratio > 1).mean()),
        "deep_penetration_warning": True,
        "reaction_ratio_interpretation": "diagnostic_only_not_a_physical_torque_requirement_because_deep_penetration_activates_penalty_stiffness",
        "case_A_supported": False,
        "case_A_reason": "prescribed replay did not produce clean hit-separate wobble",
        "startup_transient_excluded": True,
        "zeta020_magnetic_run_authorized": False,
        "decision": "DO_NOT_CHANGE_ZETA_OR_RUN_GEOMETRY_SWEEP; REDUCE_COMMAND_TRAJECTORY_TO_GEOMETRICALLY_FEASIBLE_TILT_NEXT",
    }
    (HERE / "prescribed_replay_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
