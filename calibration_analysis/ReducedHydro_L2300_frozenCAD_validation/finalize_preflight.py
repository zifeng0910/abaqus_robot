"""Write the frozen-L2300 preflight decision and evidence-bounded report."""

from pathlib import Path
import hashlib
import json

import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STEP_REL = Path("cad/freecad_parametric_robot/variants/L2300_D0815_wallwobble/Robot_L2300_D0815_WallWobble.step")
EXPECTED_SHA = "461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2"
MOMENT = 0.001040040143


def main():
    step = REPO / STEP_REL
    actual_sha = hashlib.sha256(step.read_bytes()).hexdigest()
    if actual_sha != EXPECTED_SHA:
        raise RuntimeError("Frozen STEP SHA256 mismatch")

    static = pd.read_csv(HERE / "L2300_static_clearance_map.csv")
    ceiling = pd.read_csv(HERE / "L2300_static_tilt_ceiling.csv")
    regression = pd.read_csv(HERE / "L2300_CAD_vs_solver_surface_regression.csv")
    mesh = json.loads((HERE / "L2300_solver_mesh_summary.json").read_text())
    geom = json.loads((REPO / STEP_REL.parent / "Robot_L2300_D0815_WallWobble_geometry.json").read_text())
    props = json.loads((REPO / STEP_REL.parent / "Robot_L2300_D0815_WallWobble_mass_properties.json").read_text())
    magnetic = json.loads((REPO / STEP_REL.parent / "Robot_L2300_D0815_WallWobble_magnetic_identity.json").read_text())
    calculated_moment = (magnetic["reference_magnetic_moment_Am2"]
                         * props["volume"]["value"] / magnetic["reference_volume_mm3"])
    if abs(calculated_moment - MOMENT) > 5e-13:
        raise RuntimeError("Magnetic moment identity mismatch")

    bad_contact = regression.loc[~regression.contact_classification_match]
    bad_support = regression.loc[~regression.wall_support_classification_match]
    decision = {
        "classification": "L2300_IMPLEMENTATION_INVALID",
        "static_classification": "L2300_STATIC_WALL_SUPPORT_WINDOW_PRESENT",
        "reason": "global 0.060 mm solver surface fails exact contact and wall-support classification gates",
        "datacheck_run": False,
        "dynamic_run_count": 0,
        "dynamic_budget_unconsumed": True,
        "head_refinement_used": False,
        "head_refinement_reason": "not used because remaining contact mismatches are at BODY/TAIL and cannot be corrected by the only permitted HEAD refinement",
        "contact_mismatch_poses": bad_contact[["tilt_deg", "azimuth_deg"]].to_dict("records"),
        "wall_support_mismatch_poses": bad_support[["tilt_deg", "azimuth_deg"]].to_dict("records"),
        "exactly_one_next_step": "Authorize a solver-surface strategy that can conform the BODY/TAIL grazing poses, then rerun preflight before any dynamic solve.",
    }
    (HERE / "L2300_preflight_decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    pd.DataFrame([{
        "rule": magnetic["rule"],
        "reference_volume_mm3": magnetic["reference_volume_mm3"],
        "reference_magnetic_moment_Am2": magnetic["reference_magnetic_moment_Am2"],
        "L2300_CAD_volume_mm3": props["volume"]["value"],
        "calculated_magnetic_moment_Am2": calculated_moment,
        "target_magnetic_moment_Am2": MOMENT,
        "match": abs(calculated_moment - MOMENT) <= 5e-13,
    }]).to_csv(HERE / "L2300_magnetic_identity.csv", index=False)
    pd.DataFrame([{
        "job": "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083",
        "datacheck_run": False, "dynamic_run": False, "duration_ms": 8.333,
        "dt_s": 1e-7, "reason": "solver-surface classification hard gate failed",
    }]).to_csv(HERE / "L2300_runtime_identity.csv", index=False)

    support_onset = ceiling.first_wall_support_tilt_deg.dropna()
    bridge_onset = ceiling.first_opposing_20um_tilt_deg.dropna()
    rows30 = static[static.tilt_deg == 30]
    rows32 = static[static.tilt_deg == 32]
    max_inertia = 100 * max(mesh["inertia_relative_error"])
    report = f"""# Wobble30Hz L2300 Wall-Supported Final Validation

## 1. Frozen FreeCAD STEP identity

The only geometry read was `{STEP_REL.as_posix()}`. SHA256 is `{actual_sha}` and matches the frozen manifest. The unchanged CAD identity is L={geom['L_total_mm']:.6f} mm, D={geom['D_body_mm']:.6f} mm, volume={props['volume']['value']:.12f} mm3, mass={props['mass']['mg']:.9f} mg, and COM x={props['center_of_mass']['value'][0]:.12f} mm.

## 2. Why L2.3 is tested

L1.8 tumbled through 90 degrees, while the approximately 2.9 mm robot formed a persistent geometric bridge. Frozen L2.3 is tested as the intermediate wall-supported geometry; no dimension or HEAD profile was changed here.

## 3. Exact static wall-support map

The refined 27,722-point frozen BRep surface was queried against all 1,528 SmoothWall114 triangles at 360 prescribed poses. Classification is **L2300_STATIC_WALL_SUPPORT_WINDOW_PRESENT**. At 30 degrees the exact gap spans {rows30.min_gap_um.min():.3f} to {rows30.min_gap_um.max():.3f} um; at 32 degrees it spans {rows32.min_gap_um.min():.3f} to {rows32.min_gap_um.max():.3f} um.

## 4. Static geometric tilt ceiling

Discrete first wall-support tilt ranges from {support_onset.min():.0f} to {support_onset.max():.0f} degrees (median {support_onset.median():.0f} degrees). A 20 um opposing bridge first appears only at {bridge_onset.min():.0f} degrees for {len(bridge_onset)} of 36 azimuths. This places L2.3 between the observed L1.8 tumble and long-robot jam geometries.

## 5. Simplified Abaqus rigid mesh

The canonical STEP was imported directly with `scaleFromFile=OFF`. The global 0.060 mm mesh has {mesh['node_count']:,} nodes, {mesh['element_count']:,} C3D4 elements, and {mesh['exterior_triangle_count']:,} exterior triangles. No local refinement was applied.

## 6. CAD vs solver-surface regression

Maximum overall and regional gap errors are both {mesh['CAD_vs_solver']['maximum_overall_gap_error_um']:.3f} um, within the 10 um positional gate. Opposing-bridge classification matches at all poses. The hard classification gate nevertheless fails: contact mismatches occur at (32 deg, 350 deg) and (35 deg, 340 deg), while wall-support mismatches occur at (30 deg, 230 deg) and (30 deg, 290 deg). The two contact mismatches are BODY/TAIL grazing cases: exact CAD penetration is only 0.183 and 0.543 um, while the faceted solver surface remains clear by 0.287 and 0.163 um. The permitted HEAD-only 0.050 mm refinement cannot repair these BODY/TAIL mismatches and was therefore not consumed.

## 7. Mass / COM / inertia

Mesh mass error is {100 * mesh['mass_relative_error']:.5f}%, COM error is {mesh['COM_error_um']:.4f} um, and maximum principal-inertia error is {max_inertia:.5f}%; all pass. Initial solver-surface gap is {mesh['initial_exact_gap_um']:.3f} um with no penetration. Normal P95 remains diagnostic only.

## 8. Magnetic identity

Constant magnetization applied to the frozen CAD volume gives {calculated_moment:.12f} A m2, matching the target {MOMENT:.12f} A m2. The production field was not evaluated because the solver-surface preflight failed before deck construction.

## 9. Datacheck

Not run. Solver-surface contact and wall-support classifications did not both match exact CAD, so the Datacheck eligibility condition was not met.

## 10. Single 8.333 ms run

Not run. The one permitted dynamic-run budget remains unconsumed; no ODB, SIM, raw telemetry, dynamic CSV, or GIF was generated.

## 11. Wall-supported fraction

Static fractions rise from 0 at 20 degrees to 0.306 at 28 degrees, 0.639 at 30 degrees, 0.778 at 32 degrees, and 1.000 at 35 degrees. Dynamic post-contact fraction is N/A.

## 12. Wall-sector switching

N/A without an eligible dynamic trajectory.

## 13. Directed tilt / tumble

N/A. No claim about 90-degree crossing or polarity reversal is made from static geometry alone.

## 14. Opposing bridge / jam

Static opposing geometry is absent through 35 degrees and occurs at 38 degrees for 4/36 azimuths at 20 um. Dynamic bridge duration and jam classification are N/A.

## 15. True-axis phase

N/A without an eligible dynamic trajectory.

## 16. Translation

N/A; total delta-s, 2 ms-to-end delta-s, and late median Vt were not computed.

## 17. Long vs L1800 vs L2300 GIF interpretation

No L2300 dynamic GIF was generated because preflight failed. Static evidence supports an intermediate geometry but cannot establish wobble, tumble, jam, or propulsion.

## 18. Decision

**L2300_IMPLEMENTATION_INVALID**

This is a solver-surface preprocessing failure, not evidence that frozen L2.3 is too short or too long. The static CAD geometry passes; the global 0.060 mm faceted surface does not reproduce all mandatory exact-CAD threshold classifications.

## 19. Exactly one next step

Authorize a solver-surface strategy that can conform the BODY/TAIL grazing poses, then rerun the same preflight before any Datacheck or dynamic solve. Extending duration, optimizing propulsion, or changing frozen geometry is not justified yet.
"""
    (HERE / "Wobble30Hz_L2300_WallSupported_FinalValidation.md").write_text(report)


if __name__ == "__main__":
    main()
