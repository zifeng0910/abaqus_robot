"""Apply the authorized 5 um threshold-aware solver-surface gate."""

from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
EPS_GAP_UM = 5.0


def contact_status(cad, mesh):
    if (cad <= 0) == (mesh <= 0):
        return "MATCH"
    if -EPS_GAP_UM <= cad <= EPS_GAP_UM or -EPS_GAP_UM <= mesh <= EPS_GAP_UM:
        return "CONTACT_THRESHOLD_AMBIGUOUS"
    return "TRUE_CONTACT_CLASSIFICATION_FAILURE"


def support_status(cad, mesh):
    if (cad <= 20) == (mesh <= 20):
        return "MATCH"
    if 15 <= cad <= 25 and 15 <= mesh <= 25:
        return "WALL_SUPPORT_THRESHOLD_AMBIGUOUS"
    return "TRUE_WALL_SUPPORT_CLASSIFICATION_FAILURE"


def main():
    table = pd.read_csv(HERE / "L2300_CAD_vs_solver_surface_regression.csv")
    table["EPS_GAP_um"] = EPS_GAP_UM
    table["contact_threshold_status"] = [
        contact_status(cad, mesh) for cad, mesh in zip(table.min_gap_um, table.mesh_min_gap_um)
    ]
    cad_end = table[["HEAD_min_gap_um", "TAIL_min_gap_um"]].min(axis=1)
    mesh_end = table[["mesh_HEAD_min_gap_um", "mesh_TAIL_min_gap_um"]].min(axis=1)
    table["CAD_end_gap_um"] = cad_end
    table["mesh_end_gap_um"] = mesh_end
    table["wall_support_threshold_status"] = [
        support_status(cad, mesh) for cad, mesh in zip(cad_end, mesh_end)
    ]
    table["opposing_strict_match"] = table.opposing_classification_match
    table["threshold_aware_pose_pass"] = (
        ~table.contact_threshold_status.str.contains("FAILURE")
        & ~table.wall_support_threshold_status.str.contains("FAILURE")
        & table.opposing_strict_match
    )
    table.to_csv(HERE / "L2300_threshold_aware_mesh_acceptance.csv", index=False)

    mesh = json.loads((HERE / "L2300_solver_mesh_summary.json").read_text())
    true_contact = int(table.contact_threshold_status.str.contains("FAILURE").sum())
    ambiguous_contact = int((table.contact_threshold_status == "CONTACT_THRESHOLD_AMBIGUOUS").sum())
    true_support = int(table.wall_support_threshold_status.str.contains("FAILURE").sum())
    ambiguous_support = int((table.wall_support_threshold_status == "WALL_SUPPORT_THRESHOLD_AMBIGUOUS").sum())
    gates = {
        "preferred_max_gap_error_le_5um": mesh["CAD_vs_solver"]["maximum_overall_gap_error_um"] <= 5,
        "hard_max_gap_error_le_10um": mesh["CAD_vs_solver"]["maximum_overall_gap_error_um"] <= 10,
        "no_non_ambiguous_contact_mismatch": true_contact == 0,
        "no_non_ambiguous_wall_support_mismatch": true_support == 0,
        "opposing_classification_all_match": bool(table.opposing_strict_match.all()),
        "mass_error_lt_0p5pct": mesh["mass_relative_error"] < .005,
        "COM_error_lt_5um": mesh["COM_error_um"] < 5,
        "principal_inertia_error_lt_1p5pct": max(mesh["inertia_relative_error"]) < .015,
        "mesh_frozen_29141_C3D4": mesh["element_count"] == 29141,
    }
    accepted = bool(all(gates.values()))
    summary = {
        "classification": "L2300_SOLVER_SURFACE_ACCEPTED_WITH_5UM_UNCERTAINTY" if accepted else "L2300_SOLVER_SURFACE_REJECTED",
        "EPS_GAP_um": EPS_GAP_UM,
        "pose_count": len(table),
        "maximum_CAD_solver_gap_error_um": mesh["CAD_vs_solver"]["maximum_overall_gap_error_um"],
        "contact_threshold_ambiguous_count": ambiguous_contact,
        "true_contact_failure_count": true_contact,
        "wall_support_threshold_ambiguous_count": ambiguous_support,
        "true_wall_support_failure_count": true_support,
        "opposing_mismatch_count": int((~table.opposing_strict_match).sum()),
        "gates": gates,
        "accepted": accepted,
        "mesh_refinement_performed": False,
    }
    (HERE / "L2300_threshold_aware_acceptance.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not accepted:
        raise RuntimeError(json.dumps(summary, indent=2))

    report_path = HERE / "Wobble30Hz_L2300_WallSupported_FinalValidation.md"
    report = report_path.read_text()
    marker = "## 20. Threshold-aware solver-surface acceptance"
    if marker not in report:
        report += f"""

## 20. Threshold-aware solver-surface acceptance

The prior `L2300_IMPLEMENTATION_INVALID` decision above is retained as an audit-history result of the former exact Boolean gate. The subsequently authorized engineering rule uses `EPS_GAP = 5 um`, consistent with the measured maximum CAD-to-solver gap error of {summary['maximum_CAD_solver_gap_error_um']:.3f} um. The two sub-micron contact sign changes are `CONTACT_THRESHOLD_AMBIGUOUS`; the two 20 um wall-support crossings are `WALL_SUPPORT_THRESHOLD_AMBIGUOUS`. There are zero non-ambiguous contact failures, zero non-ambiguous wall-support failures, and zero opposing-bridge mismatches across all {len(table)} poses.

Current pre-datacheck classification: **L2300_SOLVER_SURFACE_ACCEPTED_WITH_5UM_UNCERTAINTY**. The 29,141-C3D4 mesh is frozen; no refinement is performed.
"""
        report_path.write_text(report)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
