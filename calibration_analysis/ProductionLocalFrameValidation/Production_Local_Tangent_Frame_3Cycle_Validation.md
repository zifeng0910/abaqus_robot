# Production local-tangent frame: three-cycle validation

## Final classifications

- Field implementation: **`PRODUCTION_LOCAL_FRAME_IMPLEMENTATION_VALID`**.
- Three-cycle dynamics: **`LOCAL_FRAME_ROTATIONAL_WOBBLE_NOT_LIMIT_CYCLE`**.

The corrected production field follows the robot center location on the tube
centerline and rotates once per 30 Hz cycle about the local tube tangent. The
100 ms Abaqus/Explicit case develops a complete, bounded, nearly phase-locked
rotational wobble in Cycles 2 and 3. It does not meet the stronger stable-wall-
wobble criterion because the Cycle 3 contact regime and axial displacement do
not repeat Cycle 2.

## 1. Robot position and field-frame definition

At every physical increment the production server performs the following map:

```text
current robot COM in Abaqus coordinates
  -> transform to the canonical tube coordinates
  -> closest continuous line-segment projection, with previous-s continuity
     and endpoint extrapolation
  -> continuous arc coordinate s and projected centerline point
  -> parallel-transported local basis (t_hat, e1, e2)
  -> transform the local frame back to Abaqus coordinates
```

Physical reversal of `s` is retained; no monotonicity constraint is imposed.
The initial tangent polarity is selected against the configured robot axis and
then propagated continuously. This avoids independent increment-by-increment
sign choices and Frenet-frame flips.

`ROBOT_LOCAL_TANGENT` always computes the COM projection, independently of
adaptive lead. The old driver position, fixed driver arc, global Z axis and
legacy bias axis do not select the frame in this mode. Historical behavior is
retained under the explicit name `LEGACY_DRIVER_FRAME`.

## 2. Independent continuous-segment regression

The numerical oracle independently reconstructs the DXF polyline, continuous
segment projection, transported frame, field and torque. It does not call the
production projector or frame builder.

| Gate | Result | Maximum discrepancy |
|---|---:|---:|
| Legacy B/F/T, sense -1 and +1 | PASS | 0 |
| Fixed-COM B components | PASS | 0 T |
| Fixed-COM B direction | PASS | 0 deg |
| Fixed-COM field magnitude | PASS | 3.47e-18 T |
| Fixed-COM cone angle | PASS | 2.13e-14 deg |
| Fixed-COM local winding | PASS | 1.000000 turn |
| Segment interior projection | PASS | 2.22e-16 mm |
| Endpoint extrapolation | PASS | -0.25 / +0.25 mm |
| Moving projected center | PASS | 1.88e-15 mm |
| Moving s/frame/B | PASS | 0 |
| Moving torque | PASS | 1.73e-18 N mm |
| G=0 force | PASS | 0 N |

The production and deployed server SHA-256 values are identical:
`776736f93268deeef7194a9ee7ffd81ad59e8b27ac51c46f344bba60ce9b0e6d`.

The historical screening wrapper is not the numerical oracle because it snaps
to the nearest resampled DXF vertex. Its previous 2.825 deg maximum field-angle
difference is the expected difference between nearest-vertex and continuous-
segment semantics, not a production defect. As descriptive evidence, the new
continuous field reduces the maximum sampled B-direction step from 6.058 deg
to 2.006 deg and the tangent step from 6.649 deg to 2.238 deg while preserving
1,153 physical reversals of `s`.

## 3. Locked dynamic case

Exactly one dynamic case was run: `PROD_LOCAL30_G0_3CYCLE`.

| Quantity | Value |
|---|---:|
| Robot | L=2.40 mm, D=0.815 mm |
| Mesh | 1,588 nodes, 7,302 C3D4 |
| Field mode | `ROBOT_LOCAL_TANGENT` |
| Field | 10 mT, 30 Hz, 30 deg cone, sense +1 |
| Bias / gradient | 0 deg / 0 |
| Ramp | 1 ms |
| Contact | mu=0.03, zeta=0.50 |
| Reduced hydrodynamics | unchanged baseline |
| Direct increment | 1e-7 s |
| Duration | 0.100 s, three complete cycles |

Abaqus completed all 1,000,000 increments successfully in 7,899 s of solver
wall time. There was no initial penetration, socket failure, time-step change
or numerical termination. The final reported total energy is approximately
4.99e-9 in the model energy units. The largest Cycle 3 `ETOTAL` excursion,
8.83e-8, is only about 0.0035% of the approximately 2.5e-3 external work and is
associated with the resolved contact transient, not energy divergence.

## 4. Cycle-resolved response

| Metric | Cycle 1 | Cycle 2 | Cycle 3 |
|---|---:|---:|---:|
| Local B winding | 0.999897 | 0.999600 | 0.999600 |
| Robot-axis winding | 1.407365 | 1.024480 | 0.986032 |
| Mean phase lag (deg) | -1.130 | -0.242 | -0.268 |
| Tilt min / mean / max (deg) | 1.41 / 26.94 / 38.52 | 27.47 / 29.77 / 32.43 | 25.91 / 29.64 / 33.75 |
| Orbit area | 0.599971 | 0.778416 | 0.771499 |
| Orbit centroid (e1, e2) | (-0.0173, 0.0143) | (-0.0002, -0.0002) | (0.0070, -0.0065) |
| Radial RMS | 0.459779 | 0.496686 | 0.494911 |
| Contact events | 25 | 17 | 59 |
| True separations | 25 | 17 | 57 |
| Both-end support fraction | 0.0389 | 0.0300 | 0.1081 |
| Longest contact (us) | 2.012 | 2.112 | 286.661 |
| Longest 20 um bridge (ms) | 0.900 | 0.700 | 0.900 |
| Delta s (mm) | +2.014010 | -0.231320 | -1.415849 |
| ALLKE range | 1.598e-3 | 9.762e-5 | 1.139e-4 |
| ETOTAL range | 4.935e-9 | 3.060e-10 | 8.831e-8 |

Cycle 1 is acquisition, as shown by its 1.407 axis turns and low initial tilt.
Cycles 2 and 3 both show one complete bounded orbit at approximately 30 deg
tilt and near-zero phase lag. The exact wall audit has a minimum sampled gap of
0.064 um. Contact and true separation are obtained independently from dense
per-increment contact-force history, so sub-0.1-ms events are not inferred from
the coarser exact-gap samples.

## 5. Cycle 2 versus Cycle 3

The geometric orbit is repeatable: axis winding differs by 0.03845 turn,
mean phase lag by 0.0257 deg, centroid by 0.00951, radial RMS by 0.00178, orbit
area by 0.00692 and the tilt envelope by 1.56 deg. The local B winding differs
by only 1.59e-7 turn.

The wall interaction is not repeatable. Contact-sequence similarity is 0.263,
event count increases by 42, true-separation count by 40, and the longest
contact increases by 284.55 us. The dominant Cycle 3 event lasts 286.661 us
from 95.4965 to 95.7814 ms near sector 344.5 deg. Both-end support fraction
increases by 0.0781. Most decisively, per-cycle axial displacement changes by
1.18453 mm, from -0.23132 mm to -1.41585 mm.

Thus the magnetic correction restores the intended rotating wobble and the
orientation orbit is already close to periodic, but the coupled wall/translation
state is still drifting. A stable limit cycle cannot be claimed from these data.

## 6. Figures and source data

- `Production_LocalTangent_Field_OneCycle.gif`: fixed-COM field verification.
- `Continuous_vs_nearest_vertex_temporal_continuity.{png,pdf,svg}`: continuous
  production frame versus descriptive nearest-vertex history.
- `Cycle1_Cycle2_Cycle3_axis_orbit.{png,pdf,svg}`: identical-axis local orbits.
- `PROD_LOCAL30_G0_3CYCLE_DualView.gif`: fixed global tube and COM-following
  local cross-section views.
- `Cycle1_vs_Cycle2_vs_Cycle3_LocalOrbit.gif`: phase-normalized comparison.
- `three_cycle_metrics.csv`, `cycle2_vs_cycle3.csv`, `contact_events.csv`,
  `exact_gap_timeseries.csv`, `three_cycle_pose.csv` and
  `three_cycle_summary.json`: traceable quantitative outputs.

Private ODB, SIM, raw telemetry, private NPZ archives, compiled user subroutines,
cache files and TIFF images are intentionally excluded from version control.

## 7. Exactly one next step

Perform a **zero-new-solve Cycle 3 long-contact transition audit** on the
existing ODB: align the 95.4965-95.7814 ms wall resultant and impulse with local
sector, exact robot/wall geometry, hydrodynamic force and `s(t)`. This isolates
why a geometrically repeatable orientation orbit changes its wall-contact and
axial-translation regime before any physical parameter is modified.
