"""Write the public exact-FreeCAD validation report from finalized CSV files."""
from pathlib import Path
import json

import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    s = json.loads((HERE / "freecad_8p333_summary.json").read_text())
    p = json.loads((HERE / "freecad_preflight_identity.json").read_text())
    r = pd.read_csv(HERE / "freecad_head_tail_reversal_events.csv")
    b = pd.read_csv(HERE / "freecad_bridge_threshold_sensitivity.csv")
    c = pd.read_csv(HERE / "freecad_longest_bridge_composition.csv").iloc[0]
    comparison = pd.read_csv(HERE / "old_scaled_vs_freecad_summary.csv")
    old, new = comparison.iloc[0], comparison.iloc[1]
    crossing_text = ", ".join(f"{value:.4f}" for value in r.crossing_time_ms)
    sensitivity = "\n".join(
        f"- {row.threshold_um:g} um: {row.longest_bridge_ms:.3f} ms; persistent={bool(row.persistent_over_0p5ms)}"
        for row in b.itertuples())
    text = f"""# Wobble 30 Hz exact FreeCAD L1800 D0815 Reduced-Hydrodynamics validation

## 1. Executive decision

Primary classification: `{s['classification']}`.

Secondary mechanism: `{s['secondary_mechanism']}`.

The exact shape restores a 20 um opposing-wall bridge lasting {s['longest_bridge_ms']:.3f} ms, but it does **not** produce continuous solver contact. The robot also crosses transverse orientation four times.

## 2. Scope

This is the final audit of the single 8.333 ms exact-CAD run. No parameter was tuned after observing its response.

## 3. Authoritative geometry

The robot is the direct Abaqus/CAE import of `Robot_parametric_L1p800_D0p815.step`, not a bounding cylinder or analytically substituted nose.

## 4. Mesh identity

The C3D4 mesh contains {p['abaqus_import']['node_count']:,} nodes and {p['abaqus_import']['element_count']:,} elements; {p['surface_triangle_count']:,} exterior triangles define the audited surface.

## 5. Dimensions

The mesh length is {p['mesh_length_mm']:.9f} mm and diameter is {p['mesh_diameter_mm']:.9f} mm.

## 6. Mass and inertia

Mesh mass is {p['mesh_mass_mg']:.6f} mg (relative CAD error {100*p['mass_relative_error']:.4f}%). Maximum principal-inertia error is {100*max(p['inertia_relative_error']):.4f}%.

## 7. Placement and polarity

RP-to-mesh COM error is {p['RP_mesh_COM_error_um']:.4f} um and the imported HEAD-to-TAIL axis has dot product {p['head_to_tail_dot_a0']:.1f} with production `a0`.

## 8. Surface-quality gate

Exterior-normal median/P95/maximum errors are {p['surface_normal_median_deg']:.4f}/{p['surface_normal_P95_deg']:.4f}/{p['surface_normal_max_deg']:.4f} deg.

## 9. Initial clearance

The exact initial gap is +{p['initial_exact_gap_um']:.3f} um; datacheck found no overclosure or node adjustment.

## 10. Magnetic similarity

The exact-CAD moment is {p['new_moment_Am2']:.12g} A m2. New/old moment and transverse-inertia ratios are {p['moment_ratio']:.6f} and {p['new_old_transverse_inertia_ratio']:.6f}; their acceleration-scale ratio is {p['new_old_moment_per_Iperp_ratio']:.6f}. The large rotation is therefore not explained by a sudden magnetic/inertial scale increase.

## 11. Same-pose regression

The magnetic field components and 10 mT magnitude are unchanged at the same pose; force and torque scale only with magnetic moment.

## 12. Datacheck

Datacheck completed with zero errors and 13 audited warnings: 12 standard VUAMP warnings and 45 distorted tetrahedra among 78,024. The rigid-body definition is valid and has no duplicate mass.

## 13. Single dynamic run identity

Job `{s['job']}` completed {s['completed_s']*1e3:.3f} ms in {s['increments']:,} fixed increments of {s['dt_s']:.1e} s. It was the only dynamic run; there was no retry.

## 14. Sampling identity

RP history and contact history were recorded per increment. Exact geometry was scanned over the whole trajectory at 1 us and within solver-active contact windows at 0.1 us. The 10 um sensitivity uses only nodes already identified inside 20 um at saved poses.

## 15. Contact events

There are {s['contact_event_count']} solver-active events. First contact occurs at {s['first_contact_s']*1e3:.4f} ms and the longest event lasts {s['longest_contact_us']:.1f} us.

## 16. Exact gap

Minimum exact gap is {s['min_exact_gap_um']:.6f} um. Penetration is less severe than the old scaled case ({old.min_gap_um:.6f} um).

## 17. Opposing bridge definition

A bridge requires near-wall surface nodes in at least two wall-normal sectors separated by 120 deg or more. It is a geometric screen, distinct from solver contact.

## 18. Threshold sensitivity

{sensitivity}

Persistence above 0.5 ms exists only for the 20 um near-wall criterion.

## 19. Longest bridge

The longest 20 um bridge spans {c.onset_s*1e3:.3f}-{c.end_s*1e3:.3f} ms and covers {c.duration_ms:.3f} ms. Its local minimum gap is {c.minimum_gap_um:.3f} um.

## 20. Bridge composition

HEAD and TAIL participate simultaneously in {100*c.head_tail_simultaneous_fraction:.1f}% of this interval; {int(c.unique_head_near_wall_nodes)} unique HEAD nodes, {int(c.unique_tail_near_wall_nodes)} TAIL nodes, and {int(c.unique_body_near_wall_nodes)} BODY nodes participate. Maximum sector separation is {c.maximum_sector_separation_deg:.2f} deg. This is HEAD+TAIL opposing near-wall clamping, with no BODY participation.

## 21. Directed rotation and polarity

Maximum directed tilt is {s['max_directed_tilt_deg']:.6f} deg and final directed tilt is {s['final_directed_tilt_deg']:.6f} deg. The canonical increasing-s tangent is opposite production `a0` initially (`a dot t`={s['initial_axis_dot_tangent']:.6f}); this initial sign is not a dynamic reversal.

## 22. True reversal events

The four interpolated `a dot t=0` crossings occur at {crossing_text} ms. The first true reversal is {s['first_axis_reversal_s']*1e3:.4f} ms, not t=0. Each event table row records direction, field tilt, robot/B phase, contact, and 20 um bridge state.

## 23. Field and torque context

Field tilt spans {s['theta_B_min_deg']:.3f}-{s['theta_B_max_deg']:.3f} deg. Magnetic torque does not collapse. Maximum hydro power is {s['Thydro_power_max_W']:.3g} W and its positive fraction is {s['Thydro_power_positive_fraction']:.3g}, so the reduced-hydrodynamics term remains dissipative.

## 24. Translation outcome

Canonical displacement is {s['canonical_delta_s_mm']:.6f} mm and final tangent velocity is {s['final_Vt_mm_s']:.6f} mm/s. These are outcomes, not pass/fail geometry gates.

## 25. Comparison, limits, and next audit

Old scaled L1800 had {int(old.contact_events)} contact events, {old.longest_contact_us:.1f} us longest contact, {old.longest_bridge_20um_ms:.4f} ms longest bridge, and no persistent bridge. Exact CAD has {int(new.contact_events)} events, {new.longest_contact_us:.1f} us longest contact, and {new.longest_bridge_20um_ms:.4f} ms longest bridge. The exact shape improves peak penetration while worsening 20 um bridge duration.

This validates the stated classification but does not by itself prove head shape is the sole cause, because geometry and trajectory changed together. The next and only authorized step is a zero-Abaqus OLD/NEW geometry x OLD/NEW trajectory replay against the same SmoothWall114 mesh.
"""
    (HERE / "Wobble30Hz_FreeCAD_L1800_D0815_ReducedHydro_8p333_report.md").write_text(text)


if __name__ == "__main__":
    main()
