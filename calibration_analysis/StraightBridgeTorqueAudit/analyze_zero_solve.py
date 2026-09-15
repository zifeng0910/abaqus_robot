"""Zero-new-solve torque audit for the completed straight-pipe baseline."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CASE = REPO / "calibration_analysis/StraightPipeControl/case/PROD_LOCAL30_G0_STRAIGHT_CTRL"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
A0 /= np.linalg.norm(A0)


def unit(a):
    return a / max(np.linalg.norm(a), 1e-30)


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    return max((archive[key] for key in keys), key=len)


def vector_history(archive, prefix, time=None, surface=None):
    columns = []
    for component in (1, 2, 3):
        if surface:
            keys = [key for key in archive.files
                    if "|" + prefix + str(component) + " on surface " in key and surface in key]
            array = max((archive[key] for key in keys), key=len)
        else:
            array = get_array(archive, prefix + str(component))
        if time is None:
            time = array[:, 0]
        columns.append(np.interp(time, array[:, 0], array[:, 1]))
    return time, np.column_stack(columns)


def rotation_x_to_axis(axis):
    x = np.array([1.0, 0.0, 0.0])
    v = np.cross(x, axis); c = float(np.dot(x, axis))
    k = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + k + k @ k / (1.0 + c)


def ranges(flag):
    edge = np.diff(np.r_[False, flag, False].astype(int))
    return list(zip(np.where(edge == 1)[0], np.where(edge == -1)[0] - 1))


def merge_ranges(items, max_gap):
    merged = []
    for start, end in items:
        if merged and start - merged[-1][1] - 1 <= max_gap:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def main():
    rp = np.load(CASE / "private/rp_history_private.npz")
    t, u = vector_history(rp, "U")
    _, ur = vector_history(rp, "UR", t)
    _, velocity = vector_history(rp, "V", t)
    _, omega = vector_history(rp, "VR", t)
    contact = np.load(CASE / "private/contact_history_private.npz")
    _, fn = vector_history(contact, "CFN", t, "ASSEMBLY_ROBOT")
    _, fs = vector_history(contact, "CFS", t, "ASSEMBLY_ROBOT")
    fcontact = fn + fs

    identity = json.loads((CASE / "case_identity.json").read_text())
    tangent = unit(np.asarray(identity["straight_tangent_aba"]))
    e1 = unit(np.asarray(identity["initial_e1_aba"]))
    e2 = unit(np.asarray(identity["initial_e2_aba"]))
    rotations = Rotation.from_rotvec(ur)
    matrices = rotations.as_matrix()
    axis = rotations.apply(np.broadcast_to(A0, omega.shape))

    # CAD inertia is scaled to the measured 9.999973 mg solver mass. The
    # axisymmetric body frame is sufficient because the two transverse moments match.
    mass_scale = 9.999973 / 9.311813102590843
    principal_tonne_mm2 = mass_scale * np.array([
        7.596054188364871e-10, 4.456452913183043e-09, 4.456452913183047e-09])
    r0 = rotation_x_to_axis(A0)
    inertia0 = r0 @ np.diag(principal_tonne_mm2 * 1e-9) @ r0.T
    inertia = np.einsum("nij,jk,nlk->nil", matrices, inertia0, matrices)
    angular_momentum = np.einsum("nij,nj->ni", inertia, omega)

    telemetry = pd.read_csv(CASE / "PROD_LOCAL30_G0_STRAIGHT_CTRL_telemetry.csv")
    tmag = np.column_stack([np.interp(t, telemetry.t_s, telemetry[name])
                            for name in ("tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm")])
    spin = np.einsum("ij,ij->i", omega, axis)
    wobble_vector = omega - spin[:, None] * axis
    thydro = -1e-9 * spin[:, None] * axis - 3e-9 * wobble_vector

    dt = float(np.median(np.diff(t)))
    # A 10.1 us Savitzky-Golay derivative suppresses increment-scale numerical
    # chatter. Contact impulses shorter than this remain meaningful by integration,
    # not by the displayed instantaneous peak.
    window = 101
    dhdt = np.column_stack([savgol_filter(angular_momentum[:, i], window, 3,
                                         deriv=1, delta=dt, mode="interp") for i in range(3)])
    tcontact = dhdt * 1e3 - tmag - thydro

    local = lambda x: np.column_stack((x @ tangent, x @ e1, x @ e2))
    tmag_l, thydro_l, tcontact_l = map(local, (tmag, thydro, tcontact))
    gaps = pd.read_csv(REPO / "calibration_analysis/StraightPipeControl/straight_exact_gap_timeseries.csv")
    head_gap = np.interp(t, gaps.time_s, gaps.HEAD_gap_um)
    tail_gap = np.interp(t, gaps.time_s, gaps.TAIL_gap_um)
    bridge = (head_gap <= 20.0) & (tail_gap <= 20.0)
    active = np.linalg.norm(fcontact, axis=1) > 1e-8
    events = merge_ranges(ranges(active), int(round(5e-6 / dt)))

    rows = []
    for number, (start, end) in enumerate(events, 1):
        sl = slice(start, end + 1)
        impulse = np.trapz(fcontact[sl], t[sl], axis=0) if end > start else fcontact[start] * dt
        rows.append({
            "event": number, "start_s": t[start], "end_s": t[end],
            "duration_us": (end - start + 1) * dt * 1e6,
            "peak_force_N": np.linalg.norm(fcontact[sl], axis=1).max(),
            "contact_impulse_Ns": np.linalg.norm(impulse),
            "bridge_fraction": bridge[sl].mean(),
            "restitution_status": "UNRESOLVED_OPPOSING_CONTACT_RESULTANTS" if bridge[sl].any()
                                  else "NO_NODE_LEVEL_CONTACT_POINT_OUTPUT",
        })
    pd.DataFrame(rows).to_csv(HERE / "straight_contact_event_impulses.csv", index=False)

    result = pd.DataFrame({
        "time_s": t, "Tmag_t_Nmm": tmag_l[:, 0], "Tmag_e1_Nmm": tmag_l[:, 1],
        "Tmag_e2_Nmm": tmag_l[:, 2], "Tmag_norm_Nmm": np.linalg.norm(tmag, axis=1),
        "Thydro_t_Nmm": thydro_l[:, 0], "Thydro_e1_Nmm": thydro_l[:, 1],
        "Thydro_e2_Nmm": thydro_l[:, 2], "Thydro_norm_Nmm": np.linalg.norm(thydro, axis=1),
        "Tcontact_inferred_t_Nmm": tcontact_l[:, 0],
        "Tcontact_inferred_e1_Nmm": tcontact_l[:, 1],
        "Tcontact_inferred_e2_Nmm": tcontact_l[:, 2],
        "Tcontact_inferred_norm_Nmm": np.linalg.norm(tcontact, axis=1),
        "omega_spin_rad_s": spin, "omega_wobble_rad_s": np.linalg.norm(wobble_vector, axis=1),
        "HEAD_gap_um": head_gap, "TAIL_gap_um": tail_gap,
        "contact_active": active.astype(int), "opposing_bridge": bridge.astype(int),
    })
    source_sample = np.unique(np.r_[np.arange(0, len(result), 100), len(result) - 1])
    result.iloc[source_sample].to_csv(HERE / "straight_bridge_torque_balance.csv", index=False)

    analysis = t >= 0.0023
    bridged = analysis & bridge
    ratio = np.linalg.norm(tcontact[bridged, 1:], axis=1) / np.maximum(
        np.linalg.norm(tmag[bridged, 1:], axis=1), 1e-12)
    summary = {
        "classification": "PRESCRIBED_REPLAY_REQUIRED_FOR_DECISIVE_TORQUE_COMPARISON",
        "contact_torque_method": "rigid_body_angular_momentum_balance_inversion_10p1us_SG",
        "direct_contact_torque_available": False,
        "reason": "ODB has only whole-surface contact resultants; opposing contact couples are not recoverable from them",
        "analysis_start_ms": 2.3,
        "bridge_fraction": float(bridge[analysis].mean()),
        "Tmag_bridge_rms_Nmm": float(np.sqrt(np.mean(np.linalg.norm(tmag[bridged], axis=1) ** 2))),
        "Tcontact_inferred_bridge_rms_Nmm": float(np.sqrt(np.mean(np.linalg.norm(tcontact[bridged], axis=1) ** 2))),
        "inferred_transverse_contact_to_magnetic_ratio_median": float(np.median(ratio)),
        "inferred_transverse_contact_to_magnetic_ratio_p95": float(np.quantile(ratio, .95)),
        "event_count": len(events),
        "event_restitution_limit": "not identifiable for opposing simultaneous contacts without node-level force/contact-point output",
        "zeta020_authorized": False,
        "zeta_reason": "existing zeta0.50 first-impact e_n=0.621897 and gap reopening do not support damping capture",
        "historical_prescribed_case": "not found in current repository or reachable git history",
    }
    (HERE / "zero_solve_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
