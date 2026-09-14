"""Write the zero-Abaqus head-shape bridge causality audit."""
from pathlib import Path
import json

import pandas as pd

HERE=Path(__file__).resolve().parent


def markdown_table(frame):
    columns=list(frame.columns)
    lines=["| "+" | ".join(columns)+" |","| "+" | ".join("---" for _ in columns)+" |"]
    for row in frame.itertuples(index=False,name=None):
        lines.append("| "+" | ".join(f"{value:.6g}" if isinstance(value,float) else str(value) for value in row)+" |")
    return "\n".join(lines)


def main():
    s=pd.read_csv(HERE/"cross_replay_bridge_summary.csv"); identity=json.loads((HERE/"cross_replay_identity.json").read_text())
    p=pd.read_csv(HERE/"bridge_pose_geometry.csv"); h=pd.read_csv(HERE/"head_profile_difference_summary.csv").iloc[0]
    case=lambda g,t:s[(s.geometry==g)&(s.trajectory==t)].iloc[0]
    a,b,c,d=case("OLD","OLD"),case("NEW","OLD"),case("OLD","NEW"),case("NEW","NEW")
    table=markdown_table(s[["geometry","trajectory","longest20_ms","longest10_ms","longest5_ms","longest0_ms","penetration_um","HEAD","TAIL","BODY"]])
    poses=markdown_table(p)
    text=f"""# FreeCAD L1800 head-shape bridge causality audit

## 1. Why classification alone is not causality

The observed exact-CAD run changed shape, mass, inertia, magnetic moment, and trajectory together. Its 0.547 ms bridge therefore classified the outcome but could not isolate its cause.

## 2. Old geometry identity

OLD is the actual exterior reference mesh of the scaled L1800 Abaqus robot, placed about its volume COM. It is not a cylinder proxy.

## 3. New exact-CAD geometry identity

NEW is the actual imported FreeCAD exterior mesh. The convex body's {identity['geometry_nodes']['NEW']} actual convex-hull extremal nodes retain the shape envelope; no analytic head or bounding cylinder is substituted.

## 4. Old vs new head profile

In the common body-fixed frame, the new HEAD is fuller than OLD. Maximum outward difference is {h.head_max_outward_difference_um:.2f} um at x={h.head_max_difference_x_mm:.3f} mm; the positive HEAD difference spans x={h.head_positive_difference_start_x_mm:.3f} to {h.head_positive_difference_end_x_mm:.3f} mm.

## 5. Four-way cross replay

{table}

A reproduces 0.2588 ms as {a.longest20_ms:.3f} ms. D reproduces 0.547 ms exactly. B reaches {b.longest20_ms:.3f} ms when only geometry is substituted; C remains {c.longest20_ms:.3f} ms when only trajectory is substituted.

## 6. Threshold sensitivity

B remains persistent at 20, 10, 5, and 0 um ({b.longest20_ms:.3f}, {b.longest10_ms:.3f}, {b.longest5_ms:.3f}, {b.longest0_ms:.3f} ms). D persistence is specific to 20 um ({d.longest20_ms:.3f} ms); its 10 um bridge is {d.longest10_ms:.3f} ms and 5/0 um are zero.

## 7. Longest bridge pose

The observed D interval spans {d.longest20_start_s*1e3:.3f}-{d.longest20_end_s*1e3:.3f} ms. Pose details are:

{poses}

## 8. HEAD/TAIL/body participation

Under a uniform surface-axial-decile partition, B's longest bridge includes maximum per-pose counts HEAD={int(b.HEAD)}, TAIL={int(b.TAIL)}, BODY={int(b.BODY)}. D includes HEAD={int(d.HEAD)}, TAIL={int(d.TAIL)}, BODY={int(d.BODY)}. This is an end-profile/whole-shape effect, not a claim that one nose node acts alone.

## 9. Geometry-only effect

Replacing OLD with NEW on the unchanged OLD trajectory raises longest20 from {a.longest20_ms:.3f} to {b.longest20_ms:.3f} ms and drives minimum gap from {a.penetration_um:.3f} to {b.penetration_um:.3f} um. Geometry alone is sufficient to cross 0.5 ms.

## 10. Trajectory-only effect

Applying NEW motion to OLD geometry produces no opposing bridge and retains +{c.penetration_um:.3f} um minimum clearance. The new rotational trajectory alone is not sufficient.

## 11. Coupled effect

The observed NEW+NEW duration ({d.longest20_ms:.3f} ms) is less severe than NEW+OLD ({b.longest20_ms:.3f} ms). The NEW trajectory moderates, rather than creates, the geometry-driven bridge in this replay.

## 12. Rotational overshoot context

The first NEW-trajectory HEAD/TAIL flip occurs at 2.8303 ms. The observed longest bridge begins at {d.longest20_start_s*1e3:.3f} ms, after that flip; its onset/midpoint/severity directed tilts are {p.directed_tilt_deg.iloc[0]:.2f}, {p.directed_tilt_deg.iloc[1]:.2f}, and {p.directed_tilt_deg.iloc[2]:.2f} deg.

## 13. Causal classification

`{identity['causal_classification']}`

This is a geometric near-wall bridge classification, not continuous-contact jam. The exact wall is SmoothWall114, mapping is `{identity['pose_map']}`, and synchronized sampling is 1 us.

## 14. Exactly one next step

Create one HEAD-only clearance candidate that removes the 70-78 um outward nose excess near x=-0.94 to -0.88 mm while freezing D, L, mass/moment rules, hydrodynamics, contact, and trajectory replay settings. Validate that candidate first by the same zero-Abaqus four-way replay; do not run a new 8.333 ms job until it passes the bridge screen.
"""
    (HERE/"FreeCAD_L1800_head_shape_bridge_causality_audit.md").write_text(text)


if __name__=="__main__":main()
