# Wobble30Hz L2300 Wall-Supported Final Validation

## 1. Frozen FreeCAD STEP identity

The only geometry read was `cad/freecad_parametric_robot/variants/L2300_D0815_wallwobble/Robot_L2300_D0815_WallWobble.step`. SHA256 is `461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2` and matches the frozen manifest. The unchanged CAD identity is L=2.300000 mm, D=0.815000 mm, volume=1.140268022399 mm3, mass=8.904428864 mg, and COM x=1.206018693716 mm.

## 2. Why L2.3 is tested

L1.8 tumbled through 90 degrees, while the approximately 2.9 mm robot formed a persistent geometric bridge. Frozen L2.3 is tested as the intermediate wall-supported geometry; no dimension or HEAD profile was changed here.

## 3. Exact static wall-support map

The refined 27,722-point frozen BRep surface was queried against all 1,528 SmoothWall114 triangles at 360 prescribed poses. Classification is **L2300_STATIC_WALL_SUPPORT_WINDOW_PRESENT**. At 30 degrees the exact gap spans -67.961 to 73.557 um; at 32 degrees it spans -97.763 to 49.473 um.

## 4. Static geometric tilt ceiling

Discrete first wall-support tilt ranges from 25 to 35 degrees (median 30 degrees). A 20 um opposing bridge first appears only at 38 degrees for 4 of 36 azimuths. This places L2.3 between the observed L1.8 tumble and long-robot jam geometries.

## 5. Simplified Abaqus rigid mesh

The canonical STEP was imported directly with `scaleFromFile=OFF`. The global 0.060 mm mesh has 5,828 nodes, 29,141 C3D4 elements, and 3,846 exterior triangles. No local refinement was applied.

## 6. CAD vs solver-surface regression

Maximum overall and regional gap errors are both 4.025 um, within the 10 um positional gate. Opposing-bridge classification matches at all poses. The hard classification gate nevertheless fails: contact mismatches occur at (32 deg, 350 deg) and (35 deg, 340 deg), while wall-support mismatches occur at (30 deg, 230 deg) and (30 deg, 290 deg). The two contact mismatches are BODY/TAIL grazing cases: exact CAD penetration is only 0.183 and 0.543 um, while the faceted solver surface remains clear by 0.287 and 0.163 um. The permitted HEAD-only 0.050 mm refinement cannot repair these BODY/TAIL mismatches and was therefore not consumed.

## 7. Mass / COM / inertia

Mesh mass error is 0.40778%, COM error is 0.5516 um, and maximum principal-inertia error is 0.77864%; all pass. Initial solver-surface gap is 260.861 um with no penetration. Normal P95 remains diagnostic only.

## 8. Magnetic identity

Constant magnetization applied to the frozen CAD volume gives 0.001040040143 A m2, matching the target 0.001040040143 A m2. The production field was not evaluated because the solver-surface preflight failed before deck construction.

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
