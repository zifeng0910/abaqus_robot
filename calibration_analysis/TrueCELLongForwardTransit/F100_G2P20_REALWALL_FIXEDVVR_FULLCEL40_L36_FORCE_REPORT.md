# F100_G2P20_REALWALL_FIXEDVVR_FULLCEL40_L36

COUPLED_FLUID_WALL_LOAD

Kinematic replay: PASS

Maximum RP position error 5.90339e-07 mm; maximum orientation error 5.65169e-05 deg; p99 |delta v_s| 2.83136e-05 mm/s; p99 |delta omega_rock| 0.0030737 rad/s.

| Axial impulse (N s) | C1 | C2 | C3 | C4 |
|---|---:|---:|---:|---:|
| wall J_s no fluid | -7.8104e-08 | -4.3342e-08 | -8.2120e-08 | -1.1558e-07 |
| wall J_s with CEL | -7.3081e-06 | -4.1941e-06 | -1.0170e-05 | -1.7785e-05 |
| delta wall J_s | -7.2300e-06 | -4.1508e-06 | -1.0088e-05 | -1.7669e-05 |
| INFERRED ROBOT-CEL J_s | -2.3239e-08 | -7.0170e-07 | -1.0775e-05 | -3.1962e-06 |
| magnetic J_s | +2.7358e-07 | +2.8839e-07 | +2.9050e-07 | +2.9280e-07 |
| reaction J_s | +7.2535e-06 | +4.8524e-06 | +2.0862e-05 | +2.0865e-05 |

Force closure p99 residual / inertia p99: 0.00012.
The robot-CEL term is INFERRED ROBOT-CEL FORCE, calculated as whole-robot General Contact minus the direct robot-wall pair. It is not direct fluid traction.
HEAD/TAIL event forces are attributable only for exclusive-end contact; overlapping-end partitions remain unresolved.
Contact events and 10-us force/impulse traces are in adjacent CSV files. No restart continuation or free dynamics candidate was run.

## Contact events

- NO_FLUID_REAL_WALL HEAD: 6 sampled events; 1 exclusive-end events.
- NO_FLUID_REAL_WALL TAIL: 28 sampled events; 23 exclusive-end events.
- FULL_CEL_FIXED_REPLAY HEAD: 5 sampled events; 1 exclusive-end events.
- FULL_CEL_FIXED_REPLAY TAIL: 27 sampled events; 23 exclusive-end events.

The matched prescribed motion preserves the sampled contact-event times. The force differences describe the loads required along this trajectory; they do not establish a new free-dynamics trajectory or a change in collision timing.
