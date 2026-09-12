"""Numerical and delivery gates for the L=1.800 mm candidate."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    s = json.loads((HERE / "L1800_8p333_summary.json").read_text())
    identity = json.loads((HERE / "L1800_input_identity.json").read_text())
    phase = pd.read_csv(HERE / "L1800_8p333_local_phase.csv")
    axis = pd.read_csv(HERE / "L1800_8p333_true_axis.csv")
    events = pd.read_csv(HERE / "L1800_8p333_contact_events.csv")
    torque = pd.read_csv(HERE / "L1800_8p333_magnetic_hydro_torque.csv")
    comp = pd.read_csv(HERE / "baseline_vs_L1800_summary.csv")
    assert np.isclose(identity["candidate_pca_axial_span_mm"], 1.8, atol=1e-9)
    assert np.isclose(identity["base_pca_transverse_diameter_mm"], identity["candidate_pca_transverse_diameter_mm"], atol=1e-9)
    assert np.isclose(identity["mass_volume_ratio"], identity["candidate_magnetic_moment_Am2"] / identity["base_magnetic_moment_Am2"], atol=1e-12)
    assert s["completed_s"] == .008333 and s["increments"] == 83330 and s["dt_s"] == 1e-7
    assert len(phase) == len(axis) == len(torque) == 83331
    assert len(events) == s["contact_event_count"] == 18
    assert not s["persistent_bridge"] and s["longest_bridge_ms"] < .5
    assert (torque.Thydro_dot_omega_W <= 1e-12).all()
    assert np.max(np.abs(np.linalg.norm(axis[["axis_x", "axis_y", "axis_z"]], axis=1)-1)) < 1e-10
    assert set(comp.case) == {"Frozen baseline", "L=1.800 mm"}
    for name in ("L1800_phase_and_tilt", "L1800_contact_and_bridge", "L1800_torque_and_energy", "baseline_vs_L1800_dynamics"):
        assert (HERE / f"{name}.png").stat().st_size > 10_000
        assert (HERE / f"{name}.pdf").stat().st_size > 5_000
    report = (HERE / "Wobble30Hz_ReducedHydro_L1800_8p333_report.md").read_text(encoding="utf-8")
    assert "`L1800_PARTIAL_WOBBLE_SURVIVAL_NO_PERSISTENT_BRIDGE`" in report
    assert report.count("## ") == 12
    print("PASS: L1800 8.333 ms numerical and delivery checks")


if __name__ == "__main__":
    main()
