"""Create public preflight tables after the mandatory mesh-gate stop."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OLD = REPO / "calibration_analysis" / "ReducedHydro_geometry_L1800_validation"
STEM = "Robot_parametric_L2p300_D0p815_WallWobble"


def main():
    geometry=json.loads((HERE/f"{STEM}_geometry.json").read_text())
    mass=json.loads((HERE/f"{STEM}_mass_properties.json").read_text())
    mesh=json.loads((HERE/"L2300_mesh_gate.json").read_text())
    fit=json.loads((HERE/"L2300_head_fit.json").read_text())
    old=json.loads((OLD/"L1800_input_identity.json").read_text())
    old_volume=old["candidate_mass_mg"]/7.80906654321
    density=old["candidate_magnetic_moment_Am2"]/old_volume
    moment=density*geometry["volume_mm3"]
    pd.DataFrame([{
        "model":STEM,"L_total_mm":geometry["L_total_mm"],"D_body_mm":geometry["D_body_mm"],
        "head_axial_extent_mm":geometry["head_axial_extent_mm"],
        "added_straight_cylinder_mm":geometry["length_change_vs_old_L1800_mm"],
        "head_fit_max_outward_um":fit["max_outward_excess_um"],
        "head_fit_RMS_um":fit["profile_rmse_um"],"fit_hard_gate_pass":fit["fit_gate_pass"],
        "CAD_volume_mm3":geometry["volume_mm3"],"CAD_mass_mg":mass["mass"]["mg"],
    }]).to_csv(HERE/"L2300_geometry_identity.csv",index=False)
    magnetic={"constant_magnetization_rule":True,"reference":"successful old scaled L1800",
        "reference_volume_mm3":old_volume,"reference_moment_Am2":old["candidate_magnetic_moment_Am2"],
        "moment_density_Am2_per_mm3":density,"L2300_CAD_volume_mm3":geometry["volume_mm3"],
        "L2300_magnetic_moment_Am2":moment}
    pd.DataFrame([magnetic]).to_csv(HERE/"L2300_magnetic_moment.csv",index=False)
    trials=pd.DataFrame([
        {"trial":1,"deviation_factor":.20,"min_size_factor":1.0,"elements":28419,"overall_P95_deg":1.397578,"HEAD_P95_deg":3.716038,"mass_error_pct":.411677,"outcome":"FAIL_HEAD_NORMAL"},
        {"trial":2,"deviation_factor":.10,"min_size_factor":.5,"elements":48785,"overall_P95_deg":.923749,"HEAD_P95_deg":3.995883,"mass_error_pct":.195222,"outcome":"FAIL_HEAD_NORMAL"},
        {"trial":3,"deviation_factor":.05,"min_size_factor":.05,"elements":29552,"overall_P95_deg":1.397579,"HEAD_P95_deg":5.130274,"mass_error_pct":.407058,"outcome":"FAIL_HEAD_NORMAL_PRE_G2"},
        {"trial":4,"deviation_factor":.005,"min_size_factor":.05,"elements":mesh["element_count"],"overall_P95_deg":mesh["surface_normal_P95_deg"],"HEAD_P95_deg":mesh["surface_normal_by_region"]["head"]["P95_deg"],"mass_error_pct":mesh["mass_relative_error"]*100,"outcome":"FAIL_HEAD_NORMAL_FINAL"},
    ])
    trials.to_csv(HERE/"L2300_mesh_trial_summary.csv",index=False)
    decision={"classification":"L2300_MESH_GATE_FAILED_STOP","dynamic_job_run":False,
        "static_pose_scan_run":False,"reason":"HEAD surface-normal P95 is not below 2.0 deg",
        "failed_value_deg":mesh["surface_normal_by_region"]["head"]["P95_deg"],
        "hard_limit_deg":2.0,"unique_dynamic_budget_consumed":False,
        "required_next_action":"user geometry/meshing decision; no automatic length or physics change"}
    (HERE/"L2300_preflight_decision.json").write_text(json.dumps(decision,indent=2)+"\n")
    print(json.dumps({"magnetic":magnetic,"decision":decision},indent=2))


if __name__=="__main__":main()
