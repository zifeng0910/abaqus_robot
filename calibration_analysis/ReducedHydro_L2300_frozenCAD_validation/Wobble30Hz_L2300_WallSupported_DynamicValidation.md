# L2300 wall-supported wobble: dynamic validation

Primary classification: **`L2300_WALL_SUPPORTED_FORWARD_WOBBLE`**.

The unique 8.333 ms calculation shows the intended intermediate length regime. The L2.3 robot remains on its initial-polarity side, repeatedly uses the wall at a bounded 30-41 deg tilt, advances in local true-axis phase without a sustained plateau, and moves forward in canonical increasing-s. This is a quarter-cycle screening result, not evidence of 30 Hz periodic lock.

## 1. Why the old Boolean mesh gate was over-conservative

The old gate treated any sign difference across the artificial 0 um contact threshold or any class difference across the artificial 20 um support threshold as an implementation failure. All four disagreements were smaller than the measured solver-surface representation error and occurred at the threshold itself:

| Tilt | Azimuth | CAD gap | Solver gap | Revised status |
|---:|---:|---:|---:|---|
| 32 deg | 350 deg | -0.182656 um | +0.287417 um | Contact-threshold ambiguous |
| 35 deg | 340 deg | -0.542964 um | +0.163119 um | Contact-threshold ambiguous |
| 30 deg | 230 deg | 19.151227 um | 22.567383 um | Support-threshold ambiguous |
| 30 deg | 290 deg | 18.271328 um | 20.077680 um | Support-threshold ambiguous |

## 2. 5 um solver-surface uncertainty band

`EPS_GAP = 5 um` covers every mismatch. Across 360 audited poses there are two ambiguous contact classifications, two ambiguous wall-support classifications, zero non-ambiguous failures, and zero opposing-bridge mismatches. Opposing classification is therefore 100% consistent.

## 3. Frozen CAD and accepted mesh

The authoritative STEP SHA256 remains `461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2`. L and D remain 2.300000 and 0.815000 mm. The accepted global 0.060 mm rigid mesh is frozen at 5,828 nodes, 29,141 C3D4 elements, and 3,846 exterior triangles. No mesh refinement, scaling, CAD regeneration, or extra `*MASS` was performed.

## 4. Datacheck

Datacheck passed with zero errors, zero initial overclosure, zero node adjustments, a valid rigid body, one valid General Contact definition, and a successfully compiled/linked VUAMP bridge. The 13 `.dat` warnings comprise 12 standard USER-amplitude interface reminders and one distorted-element warning covering six C3D4 elements; all adjusted-node flags are `NO`. Two contact-output scope warnings and the direct-time-control energy-audit warning are nonblocking and are retained in the audit CSVs.

## 5. Single 8.333 ms solve

Exactly one dynamic job ran to 8.333 ms with fixed `dt = 1e-7 s`, 83,330 increments, and 335 field frames. There was no retry, second mesh, second geometry, or second dynamic job. Dense RP and robot-side contact histories contain 83,331 samples.

## 6. Contact events

Abaqus robot-side `CFN + CFS` identifies 263 events. First contact starts at 0.9097 ms, ends at 0.9110 ms, lasts 1.4 us, peaks at 5.727 N, and has a 4.482e-6 N s resultant impulse. Its frozen-CAD minimum is at the TAIL (`-0.710 um`); the actual wall-triangle/material-point kinematics give normal velocities of -941.45 and +520.69 mm/s and an effective restitution of 0.553 for this clean separable event.

The longest event is event 31, from 2.2981 to 2.7188 ms (420.8 us), peak 0.929 N, at the TAIL. Its contact point continues closing after the event, so it is explicitly marked non-separable and no restitution is reported. Across the catalog, contact locations are HEAD 138, TAIL 118, and BODY 7. Solver contact force remains separate from CAD geometric proximity.

## 7. Exact CAD gap reconstruction

Scientific gap authority is the frozen degree-5 CAD profile reconstructed from every RP displacement and finite rotation. The whole trajectory uses 2 us sampling; every solver-active increment uses 0.1 us local refinement. The global CAD locator has at most 2.24 um circumferential sag and every regional minimum is refined on a 2 um by 0.5 deg parameter grid, both below the accepted 5 um representation band. The trajectory minimum exact gap is `-1.403 um`.

## 8. Wall-supported fraction

Either-end support is robust across the uncertainty band: 84.045% at 15 um, 84.573% at 20 um, and 85.077% at 25 um. At the nominal 20 um threshold, HEAD support is 74.688%, TAIL support is 83.853%, and both ends are near-wall for 73.968% of the record. The four nominal states are both 73.968%, TAIL-only 9.885%, HEAD-only 0.720%, and neither 15.427%, with ten state changes on the 2 us timeline.

## 9. Wall-sector switching

The nearest supported HEAD sector spans 90.73 deg and crosses 30 deg sector bins 13 times; the TAIL sector spans 59.97 deg and crosses bins eight times. The sector timeline changes with local robot phase rather than remaining fixed at one opposing pair. HEAD/TAIL support therefore switches and moves, even though both ends are often simultaneously within the conservative 20 um proximity threshold.

## 10. Directed tilt and tumble check

Directed tilt rises from 13.34 deg to a maximum of 41.28 deg and finishes at 40.22 deg. It never crosses 90 deg and axial polarity never reverses. L2.3 therefore removes the L1.8 tumble failure over the observed quarter cycle.

## 11. Opposing bridge and jam

The longest 20 um opposing-near-wall bridge is 4.052 ms. This single Boolean condition is not a jam classification: the strong-jam gate also requires a phase plateau longer than 0.5 ms while the field continues and a locked high-tilt posture. The measured longest phase-plateau flag is only 0.0002 ms, directed posture continues to evolve, and wall sectors move. `persistent_jam` is therefore false. Relative to the long L~2.9 baseline, L2.3 removes the combined bridge-plus-phase/posture lock even though conservative geometric bridge proximity remains common.

## 12. True-axis phase

Robot phase advances `-45.879 deg` after first impact; the sign follows the chosen `(e1,e2)` handedness. Late robot and field phase rates are -15.288 and -13.864 Hz in the same convention. There is no sustained late reversal and no plateau longer than 0.5 ms. The record supports continuing phase motion, but its 8.333 ms duration cannot establish 30 Hz lock.

## 13. Torque / energy

Magnetic torque remains active: RMS `Tmag` is 0.007511 N mm before 2 ms and 0.007427 N mm afterward, so it does not collapse. Hydrodynamic power is never positive above numerical tolerance and its cumulative work is monotonic dissipative. Rigid-body kinetic energy peaks at 2.974 uJ and tracks Abaqus `ALLKE`; maximum absolute `ETOTAL` is 1.729e-10 J, only 5.81e-5 of peak kinetic energy. There is no sustained energy divergence.

## 14. Canonical translation

Canonical increasing-s is defined as FORWARD. Total displacement is `+2.5253 mm`; displacement from 2 ms to the end is `+2.3420 mm`; late median `Vt` is `+384.67 mm/s`. Both secondary propulsion gates pass.

## 15. Long vs L1800 vs L2300 GIF

The synchronized fixed-camera GIF shows the long L~2.9 case remaining opposing-wall constrained, L1.8 rotating into the tumble orientation, and L2.3 remaining wall-limited at a moderate tilt while translating and changing wall sector. The single-case GIF gives the same visual diagnosis. The observed L2.3 behavior is category C: wall-supported wobble, not free tumble and not posture/phase-locked jam.

## 16. Decision

`L2300_WALL_SUPPORTED_FORWARD_WOBBLE` is the primary label. Geometry gates pass, propulsion gates pass, and the frozen 29,141-C3D4 mesh is adequate for this rigid-body geometry screen. The result does not validate a full periodic orbit or the provisional hydrodynamic coefficients.

The requested questions resolve as follows:

1. The old gate was too strict because all four failures were threshold-edge sign/class changes smaller than the 5 um surface uncertainty.
2. The four CAD/solver pairs are listed in Section 1.
3. Yes, 5 um covers all four.
4. No non-ambiguous mismatch remains.
5. Opposing classification matches in 360/360 poses.
6. Yes, the mesh is frozen at 29,141 C3D4.
7. Datacheck passed.
8. The unique 8.333 ms run completed.
9. First contact is 0.9097 ms.
10. There are 263 solver-contact events.
11. The longest contact is 420.8 us.
12. Minimum exact CAD gap is -1.403 um.
13. Strict/nominal/loose support is 84.045%/84.573%/85.077%.
14. Both, TAIL-only, HEAD-only, and neither occupy 73.968%, 9.885%, 0.720%, and 15.427%, with ten changes.
15. Yes; supported HEAD and TAIL sectors span 90.73 and 59.97 deg with phase.
16. Maximum directed tilt is 41.28 deg.
17. It does not cross 90 deg.
18. There is no polarity reversal.
19. The longest 20 um opposing bridge is 4.052 ms.
20. There is no phase plateau longer than 0.5 ms.
21. Post-impact phase advance magnitude is 45.879 deg.
22. Late phase does not reverse relative to its post-impact/field direction.
23. `Tmag` persists.
24. Hydrodynamic work is always dissipative within tolerance.
25. Energy is bounded and `ETOTAL` is small relative to kinetic energy.
26. Total Delta s is +2.5253 mm.
27. Delta s from 2 ms to end is +2.3420 mm.
28. Late median Vt is +384.67 mm/s.
29. The GIF reads as wall-supported wobble.
30. It improves on L1.8 by eliminating the 90 deg crossing and polarity reversal.
31. It improves on L~2.9 by eliminating the combined posture/phase jam state.
32. Final classification is `L2300_WALL_SUPPORTED_FORWARD_WOBBLE`.
33. The next step is extension, not propulsion tuning or geometry modification.

## 17. Exactly one next step

Run one full 33.333 ms cycle with the identical frozen CAD, accepted mesh, contact, magnetic field, and Reduced-Hydro parameters to test whether the observed wall-supported forward wobble closes into a repeatable cycle. Do not modify geometry or tune propulsion before that extension.
