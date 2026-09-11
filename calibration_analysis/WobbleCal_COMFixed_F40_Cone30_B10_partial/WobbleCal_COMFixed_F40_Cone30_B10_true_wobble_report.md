# WobbleCal_COMFixed_F40_Cone30_B10 — true long-axis / phase-lock audit

## Status

This is a partial record: ODB ends at 73.600 ms of the requested 75 ms.
The final 1.4 ms is not used to claim completion.

## Method

No explicit HEAD/TAIL sets exist in the ODB. Fixed groups were therefore defined
once from the 5% extreme projections of the reference ROBOT_SOLID mesh along its
principal long axis (19 nodes per end), then carried with the Abaqus finite
rotation-vector Rodrigues matrix. UR1–UR3 were not interpreted as Euler angles.
The field frame is reconstructed from the production 40 Hz, cone30°, Bias40°,
tangent-axis evaluator and the validated frame rotation.

## Lock results

| window | f_robot (Hz) | R_lock | phase-lag drift (rad) | mean theta_cone (deg) |
|---|---:|---:|---:|---:|
| Cycle1 | 10.427 | 0.261 | -7.662 | 42.630 |
| Cycle2 | 39.206 | 0.980 | 1.076 | 33.149 |
| PartialCycle3 | 40.265 | 1.007 | -0.452 | 31.613 |

Overall theta_cone mean is 35.871° with peak-to-peak variation 102.819°.
HEAD/TAIL flip detected: **True**. Cycle2 versus the first common partial
Cycle3 phase span gives orbit repeat error ratio **0.4188** (preferred ≤0.10,
acceptable ≤0.15).

## Socket recovery

The main→recovery gap is 74.26 microseconds; phase-gap residual is
9.669e-15 rad and B jump is 0.000e+00 T. No phase reset is evident, so the
recovery does not contaminate the Cycle2/PartialCycle3 window.

## Verdict

**40HZ_WOBBLE_NOT_LOCKED**

This verdict concerns rotational COM-fixed calibration only. It does not validate
translation, wall contact, or a complete 75 ms run. If PASS, the next step is
COM translation release with Wall-OFF + CEL-ON for one 25 ms cycle. If NOT_LOCKED,
the next frequency is a separate 30 Hz COM-fixed calibration with all other
parameters unchanged.
