"""Write the final stopped-at-gate report from public audit products."""
from pathlib import Path
import json

import pandas as pd


HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]


def main():
    long=json.loads((REPO/"calibration_analysis"/"ReducedHydro_zeta050_8p333_validation"/"zeta050_8p333_summary.json").read_text())
    old=json.loads((REPO/"calibration_analysis"/"ReducedHydro_geometry_L1800_validation"/"L1800_8p333_summary.json").read_text())
    current=json.loads((REPO/"calibration_analysis"/"ReducedHydro_FreeCAD_L1800_validation"/"freecad_8p333_summary.json").read_text())
    dense=json.loads((REPO/"calibration_analysis"/"ReducedHydro_FreeCAD_L1800_validation"/"freecad_preflight_identity.json").read_text())
    mesh=json.loads((HERE/"current_mesh060_audit.json").read_text())
    mesh_replay=pd.read_csv(HERE/"mesh060_replay_summary.csv")
    fit=json.loads((HERE/"head_clearance_fit.json").read_text())
    geom=json.loads((HERE/"Robot_parametric_L1p800_D0p815_HeadClearance_geometry.json").read_text())
    mass=json.loads((HERE/"Robot_parametric_L1p800_D0p815_HeadClearance_mass_properties.json").read_text())
    candidate=pd.read_csv(HERE/"head_clearance_replay_summary.csv").set_index("trajectory")
    replay_id=json.loads((HERE/"head_clearance_replay_identity.json").read_text())
    datacheck=json.loads((HERE/"mesh060_datacheck_identity.json").read_text())
    e_factor=dense["abaqus_import"]["element_count"]/mesh["element_count"]
    n_factor=dense["abaqus_import"]["node_count"]/mesh["node_count"]
    inp_factor=3714205/mesh["input_size_bytes"]
    old_contact_improvement=100*(1-old["longest_contact_us"]/long["longest_contact_us"])
    current_contact_improvement=100*(1-current["longest_contact_us"]/long["longest_contact_us"])
    old_row=candidate.loc["OLD"];new_row=candidate.loc["NEW"]
    text=f"""# L1800 head clearance and mesh simplification report

## 1. Was shortening L beneficial?

Yes. Freezing `L_total=1.800 mm` and `D_body=0.815 mm` remains justified. The long baseline has {long['contact_event_count']} contact events, a {long['longest_contact_us']:.1f} us longest solver contact, and a {long['longest_bridge_ms']:.4f} ms persistent opposing bridge. Old scaled L1800 reduces these to {old['contact_event_count']} events, {old['longest_contact_us']:.1f} us, and {old['longest_bridge_ms']:.4f} ms with no persistent bridge. Exact FreeCAD L1800 has {current['contact_event_count']} events, {current['longest_contact_us']:.1f} us, and {current['longest_bridge_ms']:.4f} ms at 20 um, with no 5 or 0 um bridge.

Longest solver-contact duration improves by {old_contact_improvement:.3f}% for old L1800 and {current_contact_improvement:.3f}% for exact FreeCAD L1800 relative to the long baseline.

## 2. Three-way visual evidence

`Long_vs_L1800Scaled_vs_L1800FreeCAD_8p333.gif` synchronizes physical time, camera, center, zoom, wall opacity, B-vector scale, and robot-axis scale. Its sampling is deliberately denser at 0-1.5, 2.5-4.5, and 5.8-7.0 ms. The long case remains opposing-wall constrained while both short cases retain visibly greater clearance; the exact circular head briefly restores the 20 um HEAD+TAIL state late in the record.

`GeometryImprovement_SlowMotion_5p8_to_7p0ms.gif` isolates old scaled and current exact FreeCAD L1800. Red surface nodes and orange wall sectors expose why the fuller circular head creates the late opposing-wall state.

## 3. Why L=1.800 is frozen

Shortening reduces persistent bridge from {long['longest_bridge_ms']:.4f} ms to {old['longest_bridge_ms']:.4f} ms in the successful old geometry. No evidence supports returning to 2.9 mm or changing diameter.

## 4. Exact FreeCAD head causality

The prior four-way replay classified the mechanism as `EXACT_HEAD_GEOMETRY_DOMINATES_BRIDGE_REINTRODUCTION`: the current trajectory on old geometry has zero bridge, whereas current exact geometry on the old trajectory reaches 0.787 ms. The current head's maximum outward excess is 78.2716 um near body-fixed x=-0.943 mm.

## 5. Dense-mesh cost

The dense rigid mesh has {dense['abaqus_import']['node_count']:,} nodes, {dense['abaqus_import']['element_count']:,} C3D4, and {dense['surface_triangle_count']:,} exterior triangles. Its INP is 3,714,205 bytes and its completed 8.333 ms dynamic wall clock was 1,014 s. Because all volume elements move rigidly, internal density is computational cost rather than deformation resolution.

## 6. Simplified mesh design

The single permitted `0.060 mm` candidate has {mesh['node_count']:,} nodes, {mesh['element_count']:,} C3D4, {mesh['exterior_node_count']:,} exterior nodes, and {mesh['exterior_triangle_count']:,} exterior triangles. Preprocessing took {mesh['preprocessing_wall_clock_s']:.3f} s. Element, node, and INP reductions are {e_factor:.3f}x, {n_factor:.3f}x, and {inp_factor:.3f}x. The preferred element target below 20,000 was missed, but it was not a hard gate and geometry accuracy passed.

## 7. Simplified mesh geometry regression

Measured dimensions are {mesh['length_mm']:.9f} x {mesh['diameter_mm']:.9f} mm. Mass error is {100*mesh['mass_relative_error']:.4f}%, COM error {mesh['COM_error_um']:.4f} um, and maximum principal-inertia error {100*max(mesh['inertia_relative_error']):.4f}%. Overall surface-normal P95 is {mesh['surface_normal_P95_deg']:.4f} deg and HEAD P95 is {mesh['surface_normal_by_region']['head']['P95_deg']:.4f} deg.

On the same current trajectory, dense and simplified results are respectively {mesh_replay.iloc[0].penetration_um:.4f}/{mesh_replay.iloc[1].penetration_um:.4f} um minimum gap and {mesh_replay.iloc[0].longest20_ms:.4f}/{mesh_replay.iloc[1].longest20_ms:.4f} ms longest20. Differences are {abs(mesh_replay.iloc[0].penetration_um-mesh_replay.iloc[1].penetration_um):.4f} um and {abs(mesh_replay.iloc[0].longest20_ms-mesh_replay.iloc[1].longest20_ms):.4f} ms, so regression and bridge classification pass without a 0.050 refinement.

Datacheck has zero errors and {datacheck['distorted_C3D4']} distorted C3D4 versus {datacheck['dense_reference_distorted_C3D4']} dense. Abaqus kernel datacheck time is {datacheck['abaqus_kernel_wall_clock_s']} s versus {datacheck['dense_reference_abaqus_kernel_wall_clock_s']} s dense; the measured end-to-end invocation including compilation is {datacheck['end_to_end_wall_clock_s']:.2f} s.

## 8. Head-profile fit

The one-dimensional deterministic fit changes only the circular arc. `R_head_design={fit['fitted_R_head_design_mm']:.8f} mm`, `R_head/R_body={fit['R_head_over_R_body']:.8f}`, head axial extent is {fit['head_axial_extent_mm']:.8f} mm, and straight-body length is {fit['straight_body_length_mm']:.8f} mm. Maximum outward excess falls from 78.2716 to {fit['max_outward_excess_um']:.4f} um.

## 9. Head-clearance candidate

The candidate remains one convex G1 circular arc with unchanged cylinder, flat tail, total length, and diameter. Volume is {geom['volume_mm3']:.9f} mm3 ({100*geom['volume_ratio_vs_current_exact_CAD']:.4f}% of current exact CAD), mass is {mass['mass']['mg']:.6f} mg, COM x is {mass['center_of_mass']['value'][0]:.9f} mm, and the principal inertias are {', '.join(f'{v:.9g}' for v in mass['principal_mass_moments']['value_tonne_mm2'])} tonne mm2. Constant magnetization would scale magnetic moment with this volume if a dynamic run became eligible.

## 10. Offline bridge replay

On the OLD trajectory the candidate gives minimum gap {old_row.penetration_um:.4f} um and longest bridges 20/10/5/0 um of {old_row.longest20_ms:.4f}/{old_row.longest10_ms:.4f}/{old_row.longest5_ms:.4f}/{old_row.longest0_ms:.4f} ms. On the current exact-FreeCAD trajectory the values are {new_row.penetration_um:.4f} um and {new_row.longest20_ms:.4f}/{new_row.longest10_ms:.4f}/{new_row.longest5_ms:.4f}/{new_row.longest0_ms:.4f} ms.

The old-trajectory 20 um duration exceeds 0.5 ms and both trajectories violate the preferred -2 um minimum-gap gate. Classification: `HEAD_CLEARANCE_CANDIDATE_INSUFFICIENT`.

## 11. Dynamic eligibility

`all_zero_abaqus_gates_pass={replay_id['all_zero_abaqus_gates_pass']}`. The unique 8.333 ms dynamic job was **not run**. No candidate mesh, retry, physics change, B-field change, or second radius was attempted after failure.

## 12. One 8.333 ms result, if run

Not applicable. Contact events, solver-contact duration, dynamic longest20, dynamic exact gap, directed tilt, late phase reversal, and dynamic runtime for the HeadClearance candidate are not available because the mandatory offline gate failed.

## 13. Contact / gap / bridge

Offline replay shows that reducing the current-trajectory bridge to 0.088 ms is not sufficient: the same geometry introduces a 0.599 ms OLD-trajectory bridge and severe negative gaps. Therefore the candidate cannot be claimed to eliminate HEAD+TAIL geometric clamping robustly.

## 14. Directed rotation

No candidate dynamic rotation exists. The current exact CAD reference still reaches {current['max_directed_tilt_deg']:.4f} deg, crosses 90 deg four times, and exhibits late polarity/phase reversals. This remains context, not a candidate outcome.

## 15. Dense-vs-simplified performance

The accepted current-geometry mesh cuts elements by {e_factor:.3f}x and INP size by {inp_factor:.3f}x while preserving replay. A dynamic speedup, ODB size, and extraction-time ratio were not measured because running the candidate after gate failure was prohibited. No speedup is inferred from element count alone.

## 16. Visual GIF interpretation

The three-way GIF makes the benefit of shortening visible; the slow-motion GIF localizes the remaining exact-head bridge. These are synchronized comparisons of existing solutions and do not substitute animation appearance for exact-gap metrics.

## 17. Decision

Length shortening is confirmed and the `0.060 mm` mesh strategy is accepted for future rigid-robot work. The only HeadClearance candidate fails zero-Abaqus qualification. Geometry jam is not demonstrated as solved, and rotational overshoot cannot be reassessed for this candidate.

## 18. Exactly one next step

Before any further dynamic run, perform a geometry-only audit of why matching the old radial envelope in nose coordinates creates a severe COM-relative OLD-trajectory interference. A future head candidate requires separate authorization; B/rotational tuning is premature while geometry remains unresolved.
"""
    (HERE/"L1800_head_clearance_and_mesh_simplification_report.md").write_text(text)


if __name__=="__main__":main()
