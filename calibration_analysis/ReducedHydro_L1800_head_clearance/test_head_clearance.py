"""Regression gates for the stopped HeadClearance study."""
from pathlib import Path
import json

import pandas as pd


HERE=Path(__file__).resolve().parent


def main():
    mesh=json.loads((HERE/"current_mesh060_audit.json").read_text())
    replay=json.loads((HERE/"mesh060_replay_identity.json").read_text())
    fit=json.loads((HERE/"head_clearance_fit.json").read_text())
    candidate=json.loads((HERE/"head_clearance_replay_identity.json").read_text())
    summary=pd.read_csv(HERE/"head_clearance_replay_summary.csv").set_index("trajectory")
    assert mesh["all_required_geometry_gates_pass"]
    assert mesh["surface_normal_P95_deg"] <= 1.5
    assert mesh["surface_normal_by_region"]["head"]["P95_deg"] <= 2.0
    assert replay["all_gates_pass"]
    assert fit["fit_gate_pass"] and fit["max_outward_excess_um"] <= 10.0
    assert not candidate["all_zero_abaqus_gates_pass"]
    assert summary.loc["OLD","longest20_ms"] >= .5
    assert summary.loc["NEW","longest20_ms"] < .5
    assert not (HERE/"headclearance_dynamic_run_identity.json").exists()
    print("PASS: mesh accepted, sole HeadClearance candidate correctly stopped before dynamic")


if __name__=="__main__":main()
