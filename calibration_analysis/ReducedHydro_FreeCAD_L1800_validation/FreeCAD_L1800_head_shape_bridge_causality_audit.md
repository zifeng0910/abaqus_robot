# FreeCAD L1800 head-shape bridge causality audit

## 1. Why classification alone is not causality

The observed exact-CAD run changed shape, mass, inertia, magnetic moment, and trajectory together. Its 0.547 ms bridge therefore classified the outcome but could not isolate its cause.

## 2. Old geometry identity

OLD is the actual exterior reference mesh of the scaled L1800 Abaqus robot, placed about its volume COM. It is not a cylinder proxy.

## 3. New exact-CAD geometry identity

NEW is the actual imported FreeCAD exterior mesh. The convex body's 1419 actual convex-hull extremal nodes retain the shape envelope; no analytic head or bounding cylinder is substituted.

## 4. Old vs new head profile

In the common body-fixed frame, the new HEAD is fuller than OLD. Maximum outward difference is 78.27 um at x=-0.943 mm; the positive HEAD difference spans x=-0.943 to -0.557 mm.

## 5. Four-way cross replay

| geometry | trajectory | longest20_ms | longest10_ms | longest5_ms | longest0_ms | penetration_um | HEAD | TAIL | BODY |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OLD | OLD | 0.259 | 0.167 | 0 | 0 | -1.33398 | 3 | 2 | 0 |
| NEW | OLD | 0.787 | 0.628 | 0.645 | 0.556 | -47.6555 | 109 | 32 | 95 |
| OLD | NEW | 0 | 0 | 0 | 0 | 5.87328 | 0 | 0 | 0 |
| NEW | NEW | 0.547 | 0.038 | 0 | 0 | -0.658342 | 49 | 13 | 7 |

A reproduces 0.2588 ms as 0.259 ms. D reproduces 0.547 ms exactly. B reaches 0.787 ms when only geometry is substituted; C remains 0.000 ms when only trajectory is substituted.

## 6. Threshold sensitivity

B remains persistent at 20, 10, 5, and 0 um (0.787, 0.628, 0.645, 0.556 ms). D persistence is specific to 20 um (0.547 ms); its 10 um bridge is 0.038 ms and 5/0 um are zero.

## 7. Longest bridge pose

The observed D interval spans 6.246-6.792 ms. Pose details are:

| pose | time_s | directed_tilt_deg | local_azimuth_deg | head_gap_um | tail_gap_um | body_gap_um | wall_sector_A_deg | wall_sector_B_deg | wall_sector_separation_deg | head_nodes_20um | tail_nodes_20um | body_nodes_20um | relative_to_first_flip |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| onset | 0.006246 | 85.5593 | 1.48872 | 10.4295 | 6.11643 | 5.26216 | 17.4898 | 141.11 | 123.62 | 29 | 6 | 2 | after |
| midpoint | 0.006519 | 79.842 | 4.11125 | 0.793214 | 14.2527 | 14.7885 | 17.4378 | 167.638 | 150.2 | 55 | 3 | 3 | after |
| maximum_severity | 0.006367 | 82.2523 | 3.23622 | 17.516 | 0.0943813 | -0.344254 | 17.4259 | 163.101 | 145.675 | 10 | 17 | 9 | after |

## 8. HEAD/TAIL/body participation

Under a uniform surface-axial-decile partition, B's longest bridge includes maximum per-pose counts HEAD=109, TAIL=32, BODY=95. D includes HEAD=49, TAIL=13, BODY=7. This is an end-profile/whole-shape effect, not a claim that one nose node acts alone.

## 9. Geometry-only effect

Replacing OLD with NEW on the unchanged OLD trajectory raises longest20 from 0.259 to 0.787 ms and drives minimum gap from -1.334 to -47.655 um. Geometry alone is sufficient to cross 0.5 ms.

## 10. Trajectory-only effect

Applying NEW motion to OLD geometry produces no opposing bridge and retains +5.873 um minimum clearance. The new rotational trajectory alone is not sufficient.

## 11. Coupled effect

The observed NEW+NEW duration (0.547 ms) is less severe than NEW+OLD (0.787 ms). The NEW trajectory moderates, rather than creates, the geometry-driven bridge in this replay.

## 12. Rotational overshoot context

The first NEW-trajectory HEAD/TAIL flip occurs at 2.8303 ms. The observed longest bridge begins at 6.246 ms, after that flip; its onset/midpoint/severity directed tilts are 85.56, 79.84, and 82.25 deg.

## 13. Causal classification

`EXACT_HEAD_GEOMETRY_DOMINATES_BRIDGE_REINTRODUCTION`

This is a geometric near-wall bridge classification, not continuous-contact jam. The exact wall is SmoothWall114, mapping is `COM(t)+R(t)*(x_ref-COM_ref)`, and synchronized sampling is 1 us.

## 14. Exactly one next step

Create one HEAD-only clearance candidate that removes the 70-78 um outward nose excess near x=-0.94 to -0.88 mm while freezing D, L, mass/moment rules, hydrodynamics, contact, and trajectory replay settings. Validate that candidate first by the same zero-Abaqus four-way replay; do not run a new 8.333 ms job until it passes the bridge screen.
