# Wobble30Hz — phase-synchronised directional impulse audit

## 1. Scope and gate

This turn used the completed G6/L45 3 ms probe for zero-cost audits, then ran exactly one 3 ms Wall-ON probe after the Stage-A gates passed. No phase sweep, Gaussian/switch sweep, geometry, contact, damping, mass, frequency or cone changes were made.

## 2. Baselines

G6/L45 (`WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003`) is the reference. Its first strong CPRESS event is at 1.300 ms (node 57, 0.450 MPa); the signed exact wall gap minimum is −0.533 µm. The dynamic candidate is `Wobble_F30_PhaseDirectionalAntiImpact_WallOn_Free_003`.

## 3. Contact-gap consistency

`CPRESS` is from `CPRESS General_Contact_Domain` and is an aggregate over robot–wall and robot–CEL-fluid contact. `COPEN` is the unsupported sentinel (−3.4028e38). The exact audit uses the fixed Pipe_WALL_HELPER R3D4 triangles and the true robot exterior nodes. Therefore a positive exact wall gap and nonzero aggregate CPRESS are not contradictory: they are different contact measures/surfaces. The same caveat applies to Probe063; it must not be interpreted as wall contact solely from aggregate CPRESS.

## 4. First-impact vector geometry

The first resolved contact rows (0.8–2.0 ms) show `Vcontact·n_away` from −220 to −502 mm/s, so the robot is closing toward the wall before the 1.3 ms peak. Contact is not a single stationary node: nodes 217, 38, 34, 70, 57, 60, and 128 on wall elements 401/441/450/388/389/392/429 appear. The wall-normal/tangent projections vary (`n_away·t_centerline` −0.52 to +0.52), confirming a changing local geometry rather than one scalar normal.

## 5. Field phase map

The field remains continuous at 30 Hz. The strong CPRESS peak occurs at phase 262.05°. The selected pulse is centred at 259.05° (7° half-width), so it starts before the peak and ends at 1.67 ms. See `first_impact_phase_map.csv` and `gap_cpress_Vt_vs_field_phase.png`.

## 6. Maxwell tensor construction

The auxiliary gradient is constructed from a symmetric trace-free tensor `A` with `A*m_world=F_target` and `B_aux(COM)=0`; it is not a tangent-only projection or arbitrary body force. The selected target is `F_target = 0.400 mN*t_hat + 3.000 mN*n_away`. Frobenius norm is 3.6445 T/m, eigenvalues are −2.724, 0.324, 2.399 T/m, and reconstruction error is 1.16e−12%. At a 0.45 mm envelope, the linearized auxiliary field is about 1.64 mT. This is comparable to the previously tested local-gradient scale and remains far below a scalar 20 mT field guard.

## 7. Offline pulse estimate

The raised-cosine phase window gives positive tangential impulse 0.211 µN·s and away-from-wall impulse 1.580 µN·s. Using the measured −0.502 m/s closing speed, the predicted normal closing reduction is 31.5%, passing the Stage-A 30% gate. This is an estimate only; it does not include contact feedback.

## 8. Coordinate/sign unit test

`test_forward_force_mapping.py` passes: `[+1,0,0]·[1,0,0]=+1` and `[-1,0,0]·[1,0,0]=−1`. The authoritative initial tangent is unit length. No axis exchange was found in this synthetic test.

## 9. No-contact diagnostic

Existing Probe065 (not rerun) is diagnostic-only and gives +0.003172 mm displacement, +3.965 mm/s mean forward speed, and +19.14 mm/s peak speed with wall contact removed. Thus the magnetic bridge can accelerate in the canonical forward direction; the reverse motion in the Wall-ON probe is generated after the wall/fluid interaction is present.

## 10. Dynamic probe result

| metric | G6/L45 | directional pulse | change |
|---|---:|---:|---:|
| duration | 3.0 ms | 3.0 ms | — |
| signed arc displacement | −0.3168 mm | −0.2791 mm | 11.9% less reverse |
| final Vt | −251.1 mm/s | −225.6 mm/s | 10.1% less reverse |
| peak speed | 513.7 mm/s | 485.5 mm/s | 5.5% lower |
| aggregate CPRESS max | 0.450 MPa | 0.310 MPa | 31.1% lower |
| exact wall min gap | −0.533 µm | **+5.496 µm** | wall-gap gate improved |
| UR1 peak-to-peak | 0.847 rad | 0.749 rad | 11.6% lower |

The pulse is confirmed in the socket log and its phase window is visible in the annotated GIF. It reduced the measured reverse speed modestly and kept the exact wall gap positive, but it did not produce positive net displacement or a ≥50% Vt improvement.

## 11. Momentum accounting

For the candidate, `mΔVt = −2.256e−6 N·s`, telemetry-integrated total magnetic tangential impulse is `+0.921e−6 N·s`, of which the prescribed tensor pulse contributes `+0.259e−6 N·s`; the unresolved nonmagnetic reaction is `−3.177e−6 N·s`. The residual is identically zero because the unresolved term is defined from momentum balance. This is an accounting closure, not an independent contact-force validation. The old `Jcontact/Jmag` ratio must therefore remain retired until contact and fluid contributions are separately signed and surface-resolved.

## 12. Wobble preservation

UR1 peak-to-peak decreases only 11.6%, below the 30% preservation limit. The pulse did not cause a torque catastrophe because `B_aux(COM)=0`; the uniform rotating field remains the sole COM torque source.

## 13. Exact-wall result

The candidate minimum signed wall gap is +5.496 µm at the final frame. This passes the requested −1 µm final target and the +5 µm engineering margin, subject to the CPRESS-surface caveat in Section 3.

## 14. GIF/visual audit

`Wobble_F30_PhaseDirectionalAntiImpact_WallOn_Free_003_CEL_FSI_bend_validation.gif` uses the corrected fixed-camera Abaqus wall/centerline, robot, 3-D B vector and +s indication. Each frame is annotated with phase, pulse ON/OFF and weight, Vt proxy, exact gap, and CPRESS.

## 15. Explicit answers to the requested questions

1. **Probe063 gap/CPRESS error?** Not a physical contradiction; aggregate General Contact CPRESS and exact Pipe_WALL_HELPER gap are not the same surface quantity.
2. **Is Jcontact/Jmag trustworthy?** No; its source sign/units are inconsistent with momentum and it is retired for ranking.
3. **No-contact forward?** Yes, existing Probe065 is positive and diagnostic-only.
4. **Did the selected phase keep motion forward?** No; the single candidate remains net reverse.
5. **Was the first contact vector localized?** Yes, closing velocity is concentrated before the 1.3 ms event.
6. **Is wall normal guessed?** No; it is oriented by a positive-gap perturbation toward the authoritative centerline.
7. **Does contact normal align with centerline?** No; the dot product changes sign, demonstrating local curved geometry.
8. **What is the first-impact phase?** Approximately 262.05° at 1.300 ms.
9. **Does the drive freeze?** No; socket phase is continuous at 30 Hz.
10. **Is the tensor Maxwell-consistent?** Yes: symmetric, trace-free, zero COM auxiliary field, reconstruction error <1%.
11. **Is the pulse tangent force positive?** Yes by construction and by positive offline impulse.
12. **Is the pulse open loop?** Yes; no contact feedback or adaptive lead is used.
13. **Is the pulse comparable to tested gradients?** Yes, 3.64 T/m and ~1.64 mT over the robot envelope.
14. **Did CPRESS improve?** Aggregate maximum decreased 0.450→0.310 MPa.
15. **Did exact wall clearance improve?** Yes, −0.533→+5.496 µm.
16. **Did wobble survive?** Yes, UR1 reduction is 11.6%.
17. **Did Vt pass the 50% gate?** No; improvement is 10.1%.
18. **Did net displacement pass +0.05 mm?** No; it remains −0.279 mm.
19. **Is contact now magnetically dominated?** Not demonstrated; the unresolved nonmagnetic impulse remains larger than magnetic impulse.
20. **Was a phase sweep performed?** No; only the one gated candidate was run.
21. **Should pulse phase/amplitude be tuned after failure?** No, per the gate; do not tune this architecture.
22. **Should 8.333 ms be run now?** No; the 3 ms candidate failed the forward gate.
23. **Final decision?** **No forward phase window found** for this single directional-pulse candidate; this is not evidence that all 360° phases are impossible.

## 16. Decision and next architecture

The directional anti-impact pulse reduced early wall clearance risk and preserved wobble, but did not change the sign of net transport. No additional phase, amplitude, gradient, or Bias probes should be submitted under this architecture. The next permitted work is architecture redesign of the spatial field direction/rotation plane, preceded by a new zero-cost Maxwell/vector audit; no 8.333 ms continuation is justified from this result.

## 17. Files generated

`first_impact_vector_audit.csv`, `first_impact_phase_map.csv`, `contact_normal_tangent_geometry.csv`, `directional_tensor_design.csv`, `phase_pulse_definition.csv`, `directional_pulse_config.csv`, `phase_pulse_offline_impulse_estimate.csv`, `impulse_momentum_balance.csv`, `impulse_momentum_timeseries.csv`, dynamic motion/contact/momentum/wobble CSVs, plots, and the annotated GIF.

## 18. Reproducibility

The server change is in `J:\magpy\magpylib_socket_server.py`; the runner is `run_phase_directional_3ms.ps1`. The input deck was copied from the audited G6/L45 deck without physical parameter changes. Socket log records the pulse values, coordinate frame, and port 65498.
