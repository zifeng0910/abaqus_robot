from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def main():
    fit=json.loads((HERE/"L2300_head_fit.json").read_text())
    mesh=json.loads((HERE/"L2300_mesh_gate.json").read_text())
    stop=json.loads((HERE/"L2300_preflight_decision.json").read_text())
    assert fit["fit_gate_pass"] and fit["max_outward_excess_um"] <= 10
    assert mesh["surface_normal_P95_deg"] < 2
    assert mesh["surface_normal_by_region"]["head"]["P95_deg"] >= 2
    assert not mesh["all_mesh_gates_pass"]
    assert stop["classification"] == "L2300_MESH_GATE_FAILED_STOP"
    assert not stop["dynamic_job_run"] and not stop["unique_dynamic_budget_consumed"]
    for ext in ("step","brep","FCStd"):
        assert (HERE/f"Robot_parametric_L2p300_D0p815_WallWobble.{ext}").exists()
    print("PASS: L2300 CAD accepted and mandatory mesh-gate stop enforced")

if __name__=="__main__":main()
