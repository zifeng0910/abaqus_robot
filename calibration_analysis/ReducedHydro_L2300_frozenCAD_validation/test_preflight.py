"""Numerical regression tests for the frozen-L2300 preflight stop."""

from pathlib import Path
import hashlib
import json

import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SHA = "461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2"


def main():
    step = REPO / "cad/freecad_parametric_robot/variants/L2300_D0815_wallwobble/Robot_L2300_D0815_WallWobble.step"
    assert hashlib.sha256(step.read_bytes()).hexdigest() == SHA
    static = json.loads((HERE / "L2300_static_geometry_summary.json").read_text())
    mesh = json.loads((HERE / "L2300_solver_mesh_summary.json").read_text())
    decision = json.loads((HERE / "L2300_preflight_decision.json").read_text())
    regression = pd.read_csv(HERE / "L2300_CAD_vs_solver_surface_regression.csv")
    assert static["classification"] == "L2300_STATIC_WALL_SUPPORT_WINDOW_PRESENT"
    assert mesh["node_count"] == 5828 and mesh["element_count"] == 29141
    assert mesh["mass_relative_error"] < .005
    assert mesh["COM_error_um"] < 5
    assert max(mesh["inertia_relative_error"]) < .015
    assert mesh["CAD_vs_solver"]["maximum_regional_gap_error_um"] <= 10
    assert (~regression.contact_classification_match).sum() == 2
    assert (~regression.wall_support_classification_match).sum() == 2
    assert (~regression.opposing_classification_match).sum() == 0
    assert set(map(tuple, regression.loc[~regression.contact_classification_match,
                                          ["tilt_deg", "azimuth_deg"]].to_numpy())) == {(32, 350), (35, 340)}
    assert set(map(tuple, regression.loc[~regression.wall_support_classification_match,
                                          ["tilt_deg", "azimuth_deg"]].to_numpy())) == {(30, 230), (30, 290)}
    assert decision["classification"] == "L2300_IMPLEMENTATION_INVALID"
    assert decision["dynamic_run_count"] == 0 and decision["dynamic_budget_unconsumed"]
    print("PASS: frozen identity, static window, mesh metrics, mismatch poses, and stop decision")


if __name__ == "__main__":
    main()
