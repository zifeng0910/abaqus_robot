"""Regression gates for the completed unique L2300 dynamic analysis."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    result = json.loads((HERE/"L2300_dynamic_summary.json").read_text())
    assert result["run_completed_successfully"] and result["increments"] == 83330
    assert abs(result["completed_s"]-.008333) < 1e-12 and result["dt_s"] == 1e-7
    assert result["mesh_elements"] == 29141 and not result["mesh_refinement_performed"]
    assert result["classification"] == "L2300_WALL_SUPPORTED_FORWARD_WOBBLE"
    assert not result["axis_crossed_90deg"] and not result["polarity_reversal"]
    assert not result["persistent_jam"] and not result["phase_plateau_over_0p5ms"]
    assert not result["late_phase_reverse"] and result["late_phase_rate_Hz"]*result["late_field_phase_rate_Hz"] > 0
    assert result["strict_support_fraction"] > .80
    assert result["loose_support_fraction"]-result["strict_support_fraction"] < .02
    assert result["delta_s_2ms_to_end_mm"] > 0 and result["late_median_Vt_mm_s"] > 0
    assert result["hydro_power_max_W"] <= 1e-15 and result["hydro_power_positive_fraction"] == 0

    events = pd.read_csv(HERE/"L2300_contact_events.csv")
    assert len(events) == result["contact_event_count"] == 263
    assert abs(events.start_s.iloc[0]-result["first_contact_s"]) < 1e-12
    assert events.duration_us.max() > 400 and events.duration_us.max() < 430
    assert events.iloc[0].separable_rebound and .3 < events.iloc[0].effective_restitution < .7

    exact = pd.read_csv(HERE/"L2300_exact_gap.csv")
    coarse = exact[exact.sampling.str.startswith("coarse")]
    assert np.diff(coarse.time_s).max() <= 2.01e-6
    assert abs(exact.min_gap_um.min()-result["min_exact_gap_um"]) < 1e-10
    assert set(pd.read_csv(HERE/"Long_L1800_L2300_summary.csv").case) == {"Long L~2.9", "L1.8", "L2.3"}

    forbidden = {".odb", ".sim", ".dll", ".tiff"}
    assert not any(path.suffix.lower() in forbidden for path in HERE.iterdir())
    print("PASS: completed-run identity, wall-supported wobble gates, sampling, and delivery exclusions")


if __name__ == "__main__":
    main()
