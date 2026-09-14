# Wobble 30 Hz exact FreeCAD L1800 D0815 Reduced-Hydrodynamics validation

## 1. Executive decision

Primary classification: `EXACT_CAD_HEAD_SHAPE_REINTRODUCES_GEOMETRIC_JAM`.

Secondary mechanism: `ROTATIONAL_OVERSHOOT_PERSISTS_AND_HEAD_TAIL_POLARITY_REVERSES`.

The exact shape restores a 20 um opposing-wall bridge lasting 0.547 ms, but it does **not** produce continuous solver contact. The robot also crosses transverse orientation four times.

## 2. Scope

This is the final audit of the single 8.333 ms exact-CAD run. No parameter was tuned after observing its response.

## 3. Authoritative geometry

The robot is the direct Abaqus/CAE import of `Robot_parametric_L1p800_D0p815.step`, not a bounding cylinder or analytically substituted nose.

## 4. Mesh identity

The C3D4 mesh contains 15,215 nodes and 78,024 elements; 8,934 exterior triangles define the audited surface.

## 5. Dimensions

The mesh length is 1.799999952 mm and diameter is 0.815000027 mm.

## 6. Mass and inertia

Mesh mass is 6.753664 mg (relative CAD error 0.1517%). Maximum principal-inertia error is 0.2801%.

## 7. Placement and polarity

RP-to-mesh COM error is 0.2252 um and the imported HEAD-to-TAIL axis has dot product 1.0 with production `a0`.

## 8. Surface-quality gate

Exterior-normal median/P95/maximum errors are 0.8224/0.8224/2.4314 deg.

## 9. Initial clearance

The exact initial gap is +268.010 um; datacheck found no overclosure or node adjustment.

## 10. Magnetic similarity

The exact-CAD moment is 0.000790028260812 A m2. New/old moment and transverse-inertia ratios are 1.097811 and 1.111290; their acceleration-scale ratio is 0.987872. The large rotation is therefore not explained by a sudden magnetic/inertial scale increase.

## 11. Same-pose regression

The magnetic field components and 10 mT magnitude are unchanged at the same pose; force and torque scale only with magnetic moment.

## 12. Datacheck

Datacheck completed with zero errors and 13 audited warnings: 12 standard VUAMP warnings and 45 distorted tetrahedra among 78,024. The rigid-body definition is valid and has no duplicate mass.

## 13. Single dynamic run identity

Job `Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083` completed 8.333 ms in 83,330 fixed increments of 1.0e-07 s. It was the only dynamic run; there was no retry.

## 14. Sampling identity

RP history and contact history were recorded per increment. Exact geometry was scanned over the whole trajectory at 1 us and within solver-active contact windows at 0.1 us. The 10 um sensitivity uses only nodes already identified inside 20 um at saved poses.

## 15. Contact events

There are 19 solver-active events. First contact occurs at 0.8900 ms and the longest event lasts 1.9 us.

## 16. Exact gap

Minimum exact gap is -0.791552 um. Penetration is less severe than the old scaled case (-1.357477 um).

## 17. Opposing bridge definition

A bridge requires near-wall surface nodes in at least two wall-normal sectors separated by 120 deg or more. It is a geometric screen, distinct from solver contact.

## 18. Threshold sensitivity

- 20 um: 0.547 ms; persistent=True
- 10 um: 0.038 ms; persistent=False
- 5 um: 0.000 ms; persistent=False
- 0 um: 0.000 ms; persistent=False

Persistence above 0.5 ms exists only for the 20 um near-wall criterion.

## 19. Longest bridge

The longest 20 um bridge spans 6.246-6.792 ms and covers 0.547 ms. Its local minimum gap is -0.349 um.

## 20. Bridge composition

HEAD and TAIL participate simultaneously in 100.0% of this interval; 56 unique HEAD nodes, 28 TAIL nodes, and 0 BODY nodes participate. Maximum sector separation is 150.21 deg. This is HEAD+TAIL opposing near-wall clamping, with no BODY participation.

## 21. Directed rotation and polarity

Maximum directed tilt is 166.657693 deg and final directed tilt is 104.558215 deg. The canonical increasing-s tangent is opposite production `a0` initially (`a dot t`=-0.973009); this initial sign is not a dynamic reversal.

## 22. True reversal events

The four interpolated `a dot t=0` crossings occur at 2.8303, 4.1556, 6.1099, 7.3522 ms. The first true reversal is 2.8303 ms, not t=0. Each event table row records direction, field tilt, robot/B phase, contact, and 20 um bridge state.

## 23. Field and torque context

Field tilt spans 81.727-97.326 deg. Magnetic torque does not collapse. Maximum hydro power is 0 W and its positive fraction is 0, so the reduced-hydrodynamics term remains dissipative.

## 24. Translation outcome

Canonical displacement is -1.267754 mm and final tangent velocity is -178.296035 mm/s. These are outcomes, not pass/fail geometry gates.

## 25. Comparison, limits, and next audit

Old scaled L1800 had 18 contact events, 2.6 us longest contact, 0.2588 ms longest bridge, and no persistent bridge. Exact CAD has 19 events, 1.9 us longest contact, and 0.5470 ms longest bridge. The exact shape improves peak penetration while worsening 20 um bridge duration.

This validates the stated classification but does not by itself prove head shape is the sole cause, because geometry and trajectory changed together. The next and only authorized step is a zero-Abaqus OLD/NEW geometry x OLD/NEW trajectory replay against the same SmoothWall114 mesh.
