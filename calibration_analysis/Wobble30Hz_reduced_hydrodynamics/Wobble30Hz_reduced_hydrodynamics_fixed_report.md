# Wobble30Hz reduced hydrodynamics — fixed report

> **SUPERSEDED — do not use the RH-D / force-transfer diagnosis below.** The old momentum integration had an extra factor of 0.5, and the interval before the first 50-us-field-captured contact was incorrectly called contact-free. A nearest-wall projection bug also affected the archived gap values. The fresh clean-window balance closes within 0.077%. See [the hidden-impact audit](../ReducedHydro_hidden_impact_audit/ReducedHydro_force_transfer_and_hidden_impact_audit.md) for the actual mass, corrected distances and independently measured wall impulse. The original text below is retained solely as the historical audit trail, not as current findings.

## 1. Why ef4fc8b physical classification is invalid
The earlier result used a fixed Socket Z pose, global-XY robot phase, and `len(Series)` as a sample-count gate. It is retained only as a historical buggy run.

## 2. Socket Z-pose bug
The corrected bridge sends `x=-7.468174204284+U1`, `y=-3.676918015967+U2`, `z=-9.550745259298+U3`; force is N and torque N·mm. The old bridge omitted `+U3`.

## 3. Production pose regression
`reduced_hydro_fixed_pose_regression.csv` covers initial, translated, rotated and nonzero-U3 poses. Corrected B/F/T errors are exactly zero in the same production evaluator; the legacy fixed-Z force delta is nonzero for U3≠0.

## 4. Sensor identity
The deck declares named sensors RP_U1..U3, RP_UR1..UR3, RP_V1..V3, RP_VR1..VR3 and the bridge now reads them by name. Median |V−dU/dt| is 0.000874 mm/s; maximum is 292 mm/s at impact-resolution transitions. Direct VR is retained for hydro loading; relative-rotation omega is the wobble audit.

## 5. Mass/inertia audit
Robot target is 10 mg (`1e-5 kg`), density `7.80906654321e-9 tonne/mm^3`; the deck has no independent RP MASS card, so Ibody is the explicit cylindrical-envelope reconstruction in `reduced_hydro_fixed_mass_inertia_audit.csv`.

## 6. Linear momentum sanity
The pre-contact balance is not closed: relative residual = 1. Observed RP velocity reaches -490.5 to 0.2 mm/s while magnetic+hydro impulse is orders smaller. This is a hard implementation/force-transfer warning, not evidence for a physical wobble mechanism.

## 7. Canonical translation correction
COM is projected onto the global transformed authoritative centerline with continuous segment projection. `s_COM`, `Vt_direct` and `Vt_from_s` agree after the collision; the old integral-Vt curve is not authoritative.

## 8. Local true-axis phase correction
HEAD and TAIL are fixed exterior-node clusters; `a_hat=HEAD−TAIL`. Robot phase is `atan2(a_hat·e2,a_hat·e1)` in a deterministic local pipe frame, never global XY. Commanded field remains 30 Hz; local B phase is a moving-frame equivalent and is not the Socket frequency.

## 9. Classification gate correction
Late gating uses `n` (row count), not the number of summary columns. The final 2 ms requested late window is not available after the first exact contact at 6.150 ms, so no strong sustained-lock claim is made.

## 10. Contact windows and corrected 8.333 ms run
Near-wall onset (gap≤5 µm): 1.000 ms; exact gap≤0 onset: 6.150 ms; solver CPRESS>0.01 MPa onset: 6.150 ms. Minimum exact gap is -1.040 µm and aggregate wall CPRESS peak is 78.985 MPa.

## 11. Angular momentum and work
`H=I_global*omega`, `Jmag=∫Tmag dt` (N·mm·s), and `Wmag=∫Tmag·omega dt` (J) are in the CSVs. `Jnonmag_resolved=ΔH−Jmag` and `Wnonmag_resolved=ΔK−Wmag` are accounting identities, not independent wall-torque measurements. Hydro powers are non-positive.

## 12. Mechanism classification
**RH-D REDUCED_HYDRO_IMPLEMENTATION_INVALID**. The corrected pose and phase pipelines are valid, but the linear-momentum sanity check fails; therefore this run cannot be used as a physical wobble conclusion until the force/velocity transfer discrepancy is resolved. The first 6.15 ms are pre-exact-contact by geometry, yet the RP speed is already hundreds of mm/s under µN loads.

## 13. Required answers
1. Old bridge omitted `+U3`: yes. 2. Corrected pose matches production: yes. 3. B/F/T regression: 0. 4. Named sensor order: confirmed by deck and name-based reads. 5. Direct V vs dU/dt: median 0.000874, impact max 292 mm/s. 6–7. Mass/inertia source: audit CSV. 8–9. The old velocity jump is reproduced as a force-transfer inconsistency, not a canonical-arc artifact. 10. Pre-contact momentum does not close. 11–12. Global XY is invalid in a bend; local true-axis definition above. 13–14. Window rates are in `window_summary.csv`; late n is 5. 15–19. Contact timeline and jam timeline CSVs. 20. A >0.5 ms phase plateau cannot be diagnosed with a valid 2 ms late window. 21. GIF uses actual HEAD/TAIL and fixed camera. 22. CEL comparison is descriptive only. 23. The current 2.93×0.903 mm body reaches a wall event in the reduced model, but mechanism is unresolved by the momentum fail. 24. Do not change geometry yet; fix load/velocity transfer first.

## 14. Exactly one next step
**Resolve the reduced-model linear momentum discrepancy (CLOAD/socket force scaling and RP sensor/units) with a zero-cost force-transfer audit; do not tune wall or magnetic parameters and do not submit another physical run until `mΔV` closes against integrated external force.**
