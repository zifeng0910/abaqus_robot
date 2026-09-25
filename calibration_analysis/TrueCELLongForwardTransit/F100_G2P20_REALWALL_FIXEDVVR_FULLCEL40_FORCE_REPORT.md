# F100_G2P20_REALWALL_FIXEDVVR_FULLCEL40

COUPLED_FLUID_WALL_LOAD

Kinematic replay: PASS

Maximum RP position error 5.82626e-07 mm; maximum orientation error 5.32864e-05 deg; p99 |delta v_s| 2.63688e-05 mm/s; p99 |delta omega_rock| 0.00293773 rad/s.

| Axial impulse (N s) | C1 | C2 | C3 | C4 |
|---|---:|---:|---:|---:|
| wall J_s no fluid | -7.8113e-08 | -4.3354e-08 | -8.2114e-08 | -1.1559e-07 |
| wall J_s with CEL | -5.9565e-06 | -3.4077e-06 | -8.8144e-06 | -1.2950e-05 |
| delta wall J_s | -5.8784e-06 | -3.3644e-06 | -8.7322e-06 | -1.2835e-05 |
| INFERRED ROBOT-CEL J_s | -6.0801e-08 | -6.5336e-07 | -1.6452e-05 | -1.3698e-06 |
| magnetic J_s | +2.7358e-07 | +2.8839e-07 | +2.9050e-07 | +2.9280e-07 |
| reaction J_s | +5.9389e-06 | +4.0179e-06 | +2.5184e-05 | +1.4204e-05 |

Force closure p99 residual / inertia p99: 0.000119.
The robot-CEL term is INFERRED ROBOT-CEL FORCE, calculated as whole-robot General Contact minus the direct robot-wall pair. It is not direct fluid traction.
HEAD/TAIL event forces are attributable only for exclusive-end contact; overlapping-end partitions remain unresolved.
Contact events and 10-us force/impulse traces are in adjacent CSV files. No restart continuation or free dynamics candidate was run.

## Contact events

- NO_FLUID_REAL_WALL HEAD: 5 sampled events; 1 exclusive-end events.
- NO_FLUID_REAL_WALL TAIL: 27 sampled events; 23 exclusive-end events.
- FULL_CEL_FIXED_REPLAY HEAD: 5 sampled events; 1 exclusive-end events.
- FULL_CEL_FIXED_REPLAY TAIL: 27 sampled events; 23 exclusive-end events.

The matched prescribed motion preserves the sampled contact-event times. The force differences describe the loads required along this trajectory; they do not establish a new free-dynamics trajectory or a change in collision timing.
