"""Compare the 0.5 ms live-socket and precomputed-table dynamics."""
from __future__ import annotations

import json

import numpy as np

from table_model import OUT, config, read_table, table_loads


def load_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    cfg = config()
    live_dir = OUT / "cases" / "AB_LIVE_SOCKET_0P5MS"
    table_dir = OUT / "cases" / "AB_PRECOMPUTED_TABLE_0P5MS_CPU1"
    live = load_csv(live_dir / "AB_LIVE_SOCKET_0P5MS_rp_history.csv")
    table = load_csv(table_dir / "AB_PRECOMPUTED_TABLE_0P5MS_CPU1_rp_history.csv")
    names = ["U1", "U2", "U3", "UR1", "UR2", "UR3"]
    table_interp = {name: np.interp(live["time_s"], table["time_s"], table[name]) for name in names}
    du = np.column_stack([live[name] - table_interp[name] for name in names[:3]])
    dur = np.column_stack([live[name] - table_interp[name] for name in names[3:]])
    c = np.asarray(cfg["canonical_plus_s_axis_aba"])
    axial = np.abs(du.dot(c))
    rotation_deg = np.linalg.norm(dur, axis=1) * 180.0 / np.pi
    phase_table, s_table = read_table(OUT / cfg["table"]["path"])
    live_load = load_csv(live_dir / "AB_LIVE_SOCKET_0P5MS_telemetry.csv")
    predicted = []
    for row in live_load:
        t = float(row["t_s"])
        # Evaluate the table backend on its own dynamic trajectory, not on the
        # live trajectory. This makes the force/torque comparison a true A/B.
        position = np.asarray(cfg["initial_center_aba_mm"]) + np.array([
            np.interp(t, table["time_s"], table[name]) for name in ("U1", "U2", "U3")])
        ur = np.array([np.interp(t, table["time_s"], table[name]) for name in ("UR1", "UR2", "UR3")])
        _, _, force, torque, _, _ = table_loads(t, position, ur, cfg, phase_table, s_table)
        predicted.append(np.r_[force, torque])
    predicted = np.asarray(predicted)
    direct = np.column_stack([live_load[name] for name in ("fx_aba_N", "fy_aba_N", "fz_aba_N",
                                                           "tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm")])
    force_scale = np.maximum(np.linalg.norm(direct[:,:3],axis=1),1e-12)
    torque_scale = np.maximum(np.linalg.norm(direct[:,3:],axis=1),1e-12)
    report = {
        "duration_s": cfg["dynamic_ab_duration_s"],
        "max_axial_position_difference_mm": float(axial.max()),
        "max_rotation_vector_difference_deg": float(rotation_deg.max()),
        "magnetic_force_max_relative_error_on_live_states": float(np.max(np.linalg.norm(predicted[:,:3]-direct[:,:3],axis=1)/force_scale)),
        "magnetic_torque_max_relative_error_on_live_states": float(np.max(np.linalg.norm(predicted[:,3:]-direct[:,3:],axis=1)/torque_scale)),
        "thresholds": {"axial_position_difference_mm": 0.01, "rocking_difference_deg": 0.2},
    }
    report["passed"] = report["max_axial_position_difference_mm"] < 0.01 and report["max_rotation_vector_difference_deg"] < 0.2
    (OUT / "dynamic_ab_validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit("DYNAMIC_AB_FAILED")


if __name__ == "__main__":
    main()
