# L1800 head clearance and mesh simplification report

## 1. Was shortening L beneficial?

Yes. Freezing `L_total=1.800 mm` and `D_body=0.815 mm` remains justified. The long baseline has 204 contact events, a 415.2 us longest solver contact, and a 5.5604 ms persistent opposing bridge. Old scaled L1800 reduces these to 18 events, 2.6 us, and 0.2588 ms with no persistent bridge. Exact FreeCAD L1800 has 19 events, 1.9 us, and 0.5470 ms at 20 um, with no 5 or 0 um bridge.

Longest solver-contact duration improves by 99.374% for old L1800 and 99.542% for exact FreeCAD L1800 relative to the long baseline.

## 2. Three-way visual evidence

`Long_vs_L1800Scaled_vs_L1800FreeCAD_8p333.gif` synchronizes physical time, camera, center, zoom, wall opacity, B-vector scale, and robot-axis scale. Its sampling is deliberately denser at 0-1.5, 2.5-4.5, and 5.8-7.0 ms. The long case remains opposing-wall constrained while both short cases retain visibly greater clearance; the exact circular head briefly restores the 20 um HEAD+TAIL state late in the record.

`GeometryImprovement_SlowMotion_5p8_to_7p0ms.gif` isolates old scaled and current exact FreeCAD L1800. Red surface nodes and orange wall sectors expose why the fuller circular head creates the late opposing-wall state.

## 3. Why L=1.800 is frozen

Shortening reduces persistent bridge from 5.5604 ms to 0.2588 ms in the successful old geometry. No evidence supports returning to 2.9 mm or changing diameter.

## 4. Exact FreeCAD head causality

The prior four-way replay classified the mechanism as `EXACT_HEAD_GEOMETRY_DOMINATES_BRIDGE_REINTRODUCTION`: the current trajectory on old geometry has zero bridge, whereas current exact geometry on the old trajectory reaches 0.787 ms. The current head's maximum outward excess is 78.2716 um near body-fixed x=-0.943 mm.

## 5. Dense-mesh cost

The dense rigid mesh has 15,215 nodes, 78,024 C3D4, and 8,934 exterior triangles. Its INP is 3,714,205 bytes and its completed 8.333 ms dynamic wall clock was 1,014 s. Because all volume elements move rigidly, internal density is computational cost rather than deformation resolution.

## 6. Simplified mesh design

The single permitted `0.060 mm` candidate has 4,490 nodes, 22,282 C3D4, 1,526 exterior nodes, and 3,048 exterior triangles. Preprocessing took 2.630 s. Element, node, and INP reductions are 3.502x, 3.389x, and 2.998x. The preferred element target below 20,000 was missed, but it was not a hard gate and geometry accuracy passed.

## 7. Simplified mesh geometry regression

Measured dimensions are 1.799999952 x 0.815000036 mm. Mass error is 0.4315%, COM error 0.6037 um, and maximum principal-inertia error 0.8023%. Overall surface-normal P95 is 1.3976 deg and HEAD P95 is 1.6883 deg.

On the same current trajectory, dense and simplified results are respectively -0.6583/-0.5845 um minimum gap and 0.5470/0.5290 ms longest20. Differences are 0.0739 um and 0.0180 ms, so regression and bridge classification pass without a 0.050 refinement.

Datacheck has zero errors and 4 distorted C3D4 versus 45 dense. Abaqus kernel datacheck time is 1 s versus 2 s dense; the measured end-to-end invocation including compilation is 16.13 s.

## 8. Head-profile fit

The one-dimensional deterministic fit changes only the circular arc. `R_head_design=0.48129821 mm`, `R_head/R_body=1.18109990`, head axial extent is 0.47560676 mm, and straight-body length is 1.32439324 mm. Maximum outward excess falls from 78.2716 to 9.9989 um.

## 9. Head-clearance candidate

The candidate remains one convex G1 circular arc with unchanged cylinder, flat tail, total length, and diameter. Volume is 0.848275917 mm3 (97.9349% of current exact CAD), mass is 6.624243 mg, COM x is 0.983262698 mm, and the principal inertias are 1.76305345e-09, 1.76305345e-09, 5.28854908e-10 tonne mm2. Constant magnetization would scale magnetic moment with this volume if a dynamic run became eligible.

## 10. Offline bridge replay

On the OLD trajectory the candidate gives minimum gap -52.9121 um and longest bridges 20/10/5/0 um of 0.5990/0.3110/0.2450/0.0960 ms. On the current exact-FreeCAD trajectory the values are -5.4171 um and 0.0880/0.0000/0.0000/0.0000 ms.

The old-trajectory 20 um duration exceeds 0.5 ms and both trajectories violate the preferred -2 um minimum-gap gate. Classification: `HEAD_CLEARANCE_CANDIDATE_INSUFFICIENT`.

## 11. Dynamic eligibility

`all_zero_abaqus_gates_pass=False`. The unique 8.333 ms dynamic job was **not run**. No candidate mesh, retry, physics change, B-field change, or second radius was attempted after failure.

## 12. One 8.333 ms result, if run

Not applicable. Contact events, solver-contact duration, dynamic longest20, dynamic exact gap, directed tilt, late phase reversal, and dynamic runtime for the HeadClearance candidate are not available because the mandatory offline gate failed.

## 13. Contact / gap / bridge

Offline replay shows that reducing the current-trajectory bridge to 0.088 ms is not sufficient: the same geometry introduces a 0.599 ms OLD-trajectory bridge and severe negative gaps. Therefore the candidate cannot be claimed to eliminate HEAD+TAIL geometric clamping robustly.

## 14. Directed rotation

No candidate dynamic rotation exists. The current exact CAD reference still reaches 166.6577 deg, crosses 90 deg four times, and exhibits late polarity/phase reversals. This remains context, not a candidate outcome.

## 15. Dense-vs-simplified performance

The accepted current-geometry mesh cuts elements by 3.502x and INP size by 2.998x while preserving replay. A dynamic speedup, ODB size, and extraction-time ratio were not measured because running the candidate after gate failure was prohibited. No speedup is inferred from element count alone.

## 16. Visual GIF interpretation

The three-way GIF makes the benefit of shortening visible; the slow-motion GIF localizes the remaining exact-head bridge. These are synchronized comparisons of existing solutions and do not substitute animation appearance for exact-gap metrics.

## 17. Decision

Length shortening is confirmed and the `0.060 mm` mesh strategy is accepted for future rigid-robot work. The only HeadClearance candidate fails zero-Abaqus qualification. Geometry jam is not demonstrated as solved, and rotational overshoot cannot be reassessed for this candidate.

## 18. Exactly one next step

Before any further dynamic run, perform a geometry-only audit of why matching the old radial envelope in nose coordinates creates a severe COM-relative OLD-trajectory interference. A future head candidate requires separate authorization; B/rotational tuning is premature while geometry remains unresolved.
