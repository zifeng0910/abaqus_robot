"""Zero-dynamics audit of the production magnetic forcing architecture."""
from pathlib import Path
import json
import math

import numpy as np
import pandas as pd

from audit_common import AUDIT, REPO, RP0, A0, TransportedTubeFrame, build_production_model, local_tangent_field, unit


DATA = AUDIT / "field_only"
L2300 = REPO / "calibration_analysis/CoarseMotionModeScreen/cases/GEO_230/pose_light.csv"


def components(vector, c, e1, e2):
    return np.array([np.dot(vector, c), np.dot(vector, e1), np.dot(vector, e2)])


def cone_fit(vectors, reference_axis=None):
    directions = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    axis = unit(directions.mean(axis=0))
    if reference_axis is not None and np.dot(axis, reference_axis) < 0:
        axis *= -1.0
    angles = np.degrees(np.arccos(np.clip(directions.dot(axis), -1.0, 1.0)))
    return axis, float(angles.mean()), float(np.sqrt(np.mean((angles - angles.mean()) ** 2))), angles


def current_cycle(frame):
    production, model, wrapped, transform = build_production_model()
    times = np.linspace(0.0, 1.0 / 30.0, 1201)
    _, _, c, e1, e2 = frame.at(RP0)
    rows = []
    for time_s in times:
        result = wrapped.evaluate(time_s, RP0, np.zeros(3))
        field = np.asarray(result["B_aba_vec_T"])
        local = components(field, c, e1, e2)
        field_unit = unit(field)
        torque = np.asarray(result["torque_Nmm"])
        rows.append([time_s, result["instantaneous_phase_rad"], result["drive_scale"], *field, *local,
                     math.degrees(math.acos(np.clip(np.dot(field_unit, c), -1.0, 1.0))),
                     math.atan2(local[2], local[1]), *torque])
    columns = ["time_s", "command_phase_rad", "drive_scale", "Bx_T", "By_T", "Bz_T", "B_t_T", "B_e1_T", "B_e2_T",
               "theta_B_t_deg", "psi_B_rad", "Tx_Nmm", "Ty_Nmm", "Tz_Nmm"]
    table = pd.DataFrame(rows, columns=columns)
    table["psi_B_unwrapped_rad"] = np.unwrap(table.psi_B_rad)
    table.to_csv(DATA / "current_production_field_one_cycle.csv", index=False)

    field = table[["Bx_T", "By_T", "Bz_T"]].to_numpy()
    axis, half_angle, residual, angles = cone_fit(field)
    local_axis = components(axis, c, e1, e2)
    start_driver_arc = model._distance_at_time(0.0)
    driver_tangent_mag = model._tangent_at_distance(start_driver_arc)
    driver_normal_mag = model._z_reference_normal(driver_tangent_mag)
    biased_axis_mag = unit(math.cos(math.radians(40.0)) * driver_tangent_mag + math.sin(math.radians(40.0)) * driver_normal_mag)
    biased_axis_aba = transform.vector_to_aba(biased_axis_mag)

    m_body_aba = unit(transform.vector_to_aba(model.robot_polarity * model.robot_axis_global))
    raw_torque = np.cross(model.robot_moment_Am2 * m_body_aba, field) * 1000.0
    torque_local = np.column_stack([raw_torque.dot(c), raw_torque.dot(e1), raw_torque.dot(e2)])
    startup = table.loc[table.time_s <= 0.002].copy()
    startup[["Traw_t_Nmm", "Traw_e1_Nmm", "Traw_e2_Nmm"]] = torque_local[:len(startup)]
    startup.to_csv(DATA / "startup_torque_first_2ms.csv", index=False)

    winding = float((table.psi_B_unwrapped_rad.iloc[-1] - table.psi_B_unwrapped_rad.iloc[0]) / (2.0 * math.pi))
    summary = {
        "production_parameters": {"frequency_Hz": 30.0, "B0_T": 0.010, "cone_half_angle_deg": 30.0,
                                  "cone_axis_bias_deg": 40.0, "phase0_deg": 248.0, "sense": 1,
                                  "gradient_mT": 6.0, "gradient_length_mm": 45.0,
                                  "ramp_time_ms": 1.0, "analytic_follow_robot": True, "adaptive_lead_mm": 0.0},
        "initial_local_frame": {"c_hat_aba": c.tolist(), "e1_aba": e1.tolist(), "e2_aba": e2.tolist(),
                                "c_dot_a0": float(np.dot(c, A0))},
        "initial_driver_arc_mm": start_driver_arc,
        "initial_robot_arc_reported_by_production_model": None,
        "production_basis_selection": "DRIVER_ARC because robot_arc_mm is only populated when adaptive_lead_mm > 0",
        "fitted_cone_axis_aba": axis.tolist(),
        "fitted_cone_axis_local_c_e1_e2": local_axis.tolist(),
        "fitted_half_angle_deg": half_angle,
        "fitted_angular_residual_rms_deg": residual,
        "implemented_initial_biased_axis_aba": biased_axis_aba.tolist(),
        "implemented_axis_vs_fitted_deg": float(np.degrees(np.arccos(np.clip(np.dot(axis, biased_axis_aba), -1.0, 1.0)))),
        "fitted_axis_vs_robot_local_tangent_deg": float(np.degrees(np.arccos(np.clip(np.dot(axis, c), -1.0, 1.0)))),
        "theta_B_t_min_mean_max_deg": [float(table.theta_B_t_deg.min()), float(table.theta_B_t_deg.mean()), float(table.theta_B_t_deg.max())],
        "local_B_azimuth_winding": winding,
        "local_B_azimuth_frequency_Hz": winding * 30.0,
        "B_magnitude_T_min_mean_max": [float(np.linalg.norm(field, axis=1).min()), float(np.linalg.norm(field, axis=1).mean()), float(np.linalg.norm(field, axis=1).max())],
        "B_parallel_T_min_max": [float(table.B_t_T.min()), float(table.B_t_T.max())],
        "B_perp_T_min_max": [float(np.hypot(table.B_e1_T, table.B_e2_T).min()), float(np.hypot(table.B_e1_T, table.B_e2_T).max())],
        "m_body_aba_normalized": m_body_aba.tolist(),
        "m_body_dot_robot_HEAD_to_TAIL_a0": float(np.dot(m_body_aba, A0)),
        "initial_angle_m_B_deg": float(np.degrees(np.arccos(np.clip(np.dot(m_body_aba, unit(field[0])), -1.0, 1.0)))),
        "initial_raw_full_amplitude_torque_Nmm": raw_torque[0].tolist(),
        "initial_raw_full_amplitude_torque_local_t_e1_e2_Nmm": torque_local[0].tolist(),
        "initial_applied_torque_Nmm": table.loc[0, ["Tx_Nmm", "Ty_Nmm", "Tz_Nmm"]].tolist(),
        "initial_drive_scale": float(table.drive_scale.iloc[0]),
    }
    (DATA / "current_production_field_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def local_cycles(frame):
    times = np.linspace(0.0, 1.0 / 30.0, 1201)
    for alpha in (30.0, 40.0):
        rows = []
        for time_s in times:
            field, phase, c, e1, e2 = local_tangent_field(RP0, time_s, alpha, frame)
            local = components(field, c, e1, e2)
            rows.append([time_s, phase, *field, *local,
                         math.degrees(math.acos(np.clip(np.dot(unit(field), c), -1.0, 1.0))),
                         math.atan2(local[2], local[1])])
        table = pd.DataFrame(rows, columns=["time_s", "command_phase_rad", "Bx_T", "By_T", "Bz_T", "B_t_T", "B_e1_T", "B_e2_T", "theta_B_t_deg", "psi_B_rad"])
        table["psi_B_unwrapped_rad"] = np.unwrap(table.psi_B_rad)
        table.to_csv(DATA / ("local_tangent_field_%02ddeg_one_cycle.csv" % alpha), index=False)


def trajectory_audit(frame):
    pose = pd.read_csv(L2300)
    rows = []
    for row in pose.itertuples():
        position = np.array([row.rp_x_mm, row.rp_y_mm, row.rp_z_mm])
        _, _, c, e1, e2 = frame.at(position)
        field = np.array([row.B_x_T, row.B_y_T, row.B_z_T])
        local = components(field, c, e1, e2)
        rows.append([row.time_s, row.cycle_fraction, *position, *field, *local,
                     math.degrees(math.acos(np.clip(np.dot(unit(field), c), -1.0, 1.0))),
                     math.atan2(local[2], local[1]), math.radians(248.0) + 2.0 * math.pi * 30.0 * row.time_s])
    table = pd.DataFrame(rows, columns=["time_s", "cycle_fraction", "rp_x_mm", "rp_y_mm", "rp_z_mm", "Bx_T", "By_T", "Bz_T", "B_t_T", "B_e1_T", "B_e2_T", "theta_B_t_deg", "psi_B_rad", "command_phase_rad"])
    table["psi_B_unwrapped_rad"] = np.unwrap(table.psi_B_rad)
    table["local_field_phase_rate_Hz"] = np.gradient(table.psi_B_unwrapped_rad, table.time_s) / (2.0 * math.pi)
    table.to_csv(DATA / "production_field_along_existing_L2300.csv", index=False)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    frame = TransportedTubeFrame()
    summary = current_cycle(frame)
    local_cycles(frame)
    trajectory_audit(frame)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
