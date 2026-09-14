from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    s = json.loads((HERE / "freecad_8p333_summary.json").read_text())
    crossings = pd.read_csv(HERE / "freecad_head_tail_reversal_events.csv")
    sensitivity = pd.read_csv(HERE / "freecad_bridge_threshold_sensitivity.csv")
    torque = pd.read_csv(HERE / "freecad_8p333_torque.csv")
    assert s["first_axis_reversal_s"] > 0
    assert abs(s["first_axis_reversal_s"] - 0.0028303) < 1e-6
    assert len(crossings) == s["axis_zero_crossing_count"] == 4
    assert s["longest_contact_us"] < 10
    assert abs(s["min_exact_gap_um"] + 0.791552) < 0.002
    value20 = sensitivity.loc[sensitivity.threshold_um == 20, "longest_bridge_ms"].iloc[0]
    assert abs(value20 - 0.547) < 0.002
    assert sensitivity.loc[sensitivity.threshold_um == 10, "longest_bridge_ms"].iloc[0] < 0.5
    assert np.all(sensitivity.loc[sensitivity.threshold_um.isin([5, 0]), "longest_bridge_ms"] == 0)
    assert s["max_directed_tilt_deg"] > 150
    assert s["final_directed_tilt_deg"] > 90
    assert s["Thydro_power_max_W"] <= 1e-15
    assert torque.Thydro_dot_omega_W.max() <= 1e-15
    assert s["contact_event_count"] == 19
    assert s["initial_axis_polarity"] == -1
    assert "1 us" in sensitivity.sampling_identity.iloc[0]
    assert s["classification"] == "EXACT_CAD_HEAD_SHAPE_REINTRODUCES_GEOMETRIC_JAM"
    print("PASS: exact-FreeCAD validation numerical gates")


if __name__ == "__main__":
    main()
