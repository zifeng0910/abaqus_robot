"""Numerical and delivery regression tests for the 8.333 ms validation."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    summary = json.loads((HERE / "zeta050_8p333_summary.json").read_text())
    runtime = pd.read_csv(HERE / "zeta050_8p333_runtime_identity.csv").iloc[0]
    phase = pd.read_csv(HERE / "zeta050_8p333_local_phase.csv")
    axis = pd.read_csv(HERE / "zeta050_8p333_true_axis.csv")
    torque = pd.read_csv(HERE / "zeta050_8p333_magnetic_hydro_torque.csv")
    translation = pd.read_csv(HERE / "zeta050_8p333_canonical_translation.csv")
    events = pd.read_csv(HERE / "zeta050_8p333_contact_events.csv")
    comparison = pd.read_csv(HERE / "cel_vs_reducedhydro_zeta050_summary.csv")

    assert np.isclose(runtime.completed_s, 0.008333, atol=2e-12)
    assert int(runtime.increments) == 83330
    assert np.isclose(runtime.dt_s, 1e-7, atol=1e-15)
    assert len(phase) == len(axis) == len(torque) == len(translation) == 83331
    assert np.all(np.diff(phase.time_s) > 0)
    assert np.all(np.isfinite(axis[["axis_x", "axis_y", "axis_z", "tilt_deg"]]))
    assert np.max(np.abs(np.linalg.norm(axis[["axis_x", "axis_y", "axis_z"]], axis=1) - 1)) < 1e-10
    assert (torque.Thydro_dot_omega_W <= 1e-12).all()
    assert int(summary["contact_event_count"]) == len(events)
    assert set(comparison.case) == {"CEL", "RH old contact", "RH zeta0.50"}
    assert np.isclose(summary["canonical_delta_s_mm"], translation.delta_s_mm.iloc[-1])
    for name in (
        "true_axis_phase_8p333", "phase_difference_8p333", "tilt_angle_8p333",
        "gap_and_contact_8p333", "bridge_state_8p333", "magnetic_vs_hydro_torque_8p333",
        "canonical_translation_8p333", "cel_vs_rh_zeta050_phase",
    ):
        assert (HERE / f"{name}.png").stat().st_size > 10_000
        assert (HERE / f"{name}.pdf").stat().st_size > 5_000
    report = (HERE / "Wobble30Hz_ReducedHydro_Zeta050_8p333_report.md").read_text(encoding="utf-8")
    for section in range(1, 21):
        assert f"## {section}." in report
    decisions = [
        "REDUCED_HYDRO_ZETA050_WOBBLE_SURVIVES_8P333MS",
        "REDUCED_HYDRO_ZETA050_PARTIAL_WOBBLE_SURVIVAL",
        "GEOMETRY_CONTACT_STILL_JAMS_WITH_REDUCED_HYDRO",
        "REDUCED_HYDRO_DAMPING_OR_ROTATIONAL_MODEL_LIMITS_WOBBLE",
        "ZETA050_LONG_RUN_IMPLEMENTATION_INVALID",
    ]
    assert sum(f"`{decision}`" in report for decision in decisions) == 1
    print("PASS: zeta050 8.333 ms numerical and delivery checks")


if __name__ == "__main__":
    main()
