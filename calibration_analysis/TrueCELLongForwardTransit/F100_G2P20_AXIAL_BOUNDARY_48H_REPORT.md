# F100_G2P20 Axial Boundary Campaign

## Primary classification

`AXIAL_BOUNDARY_SENSITIVITY_UNRESOLVED`

L24 completed as a single-step, 0-40 ms, full-CEL prescribed GLOBAL V plus spatial VR replay. The only model change from L12 was the axial Eulerian extent, expanded from `[-6, 6]` mm to `[-12, 12]` mm with the physical wall helper extended to the same range.

## Verified L24 state

- Solver: completed successfully at 40.000 ms; no restart read; one dynamics run.
- Kinematic gate: PASS.
- Maximum RP position error: `6.10998e-7 mm`.
- Maximum orientation error: `5.36263e-5 deg`.
- p99 axial speed error: `2.58754e-5 mm/s`.
- p99 rocking angular-speed error: `0.00272388 rad/s`.
- Force closure p99 residual / inertia: `1.05193e-4`.
- Input SHA256: `EB2BCB1B246253E1C70540D81D86B4FE433C677522982E1935BF55D2901C23D3`.
- V/VR driver SHA256: `4EF12426EE062B576DB9AF02E12ECA750EEBC3AC8A8AE1C948C093EDAF66C8D8`.

## Verified L36 state

- Solver: completed successfully at 40.000 ms; no restart read; one dynamics run.
- Axial CEL extent: `[-18, 18]` mm, with the physical wall helper extended to the same range.
- Kinematic gate: PASS.
- Maximum RP position error: `5.90339e-7 mm`.
- Maximum orientation error: `5.65169e-5 deg`.
- p99 axial speed error: `2.83136e-5 mm/s`.
- p99 rocking angular-speed error: `0.00307370 rad/s`.
- Force closure p99 residual / inertia: `1.19734e-4`.
- Input SHA256: `E0DE815D34EAF97A08209828D70ABA232B4C813AFD1EDF2D7A18D321D8D3095B`.
- V/VR driver SHA256: `4EF12426EE062B576DB9AF02E12ECA750EEBC3AC8A8AE1C948C093EDAF66C8D8`.

## Boundary result

The inferred robot-CEL force is defined as whole-robot General Contact minus the direct robot-wall pair. Relative to L12, axial impulse changes were `D_CEL = 0.397, 0.913, 1.039, 0.023` for C1-C4 and `D_wall = 0.105, 0.130, 0.239, 0.158`. C2 and C3 therefore remain materially domain-sensitive. Near-robot and end-region fluid velocity and EVF proxies also differ between L12 and L24.

The six native field frames contain velocity, EVF, and water stress only. Pressure is reported as the explicitly labelled mean-normal-stress proxy `-trace(S_water)/3`; no native PRESS or face mass-flux history exists.

## L24 versus L36 decision

The L36 inferred robot-CEL axial impulses (C1-C4) are `-2.3239e-08`, `-7.0170e-07`, `-1.0775e-05`, and `-3.1962e-06 N s`, versus L24 values `-3.6646e-08`, `-7.5419e-06`, `+6.4695e-07`, and `-1.4021e-06 N s`. The C2 and C3 terms change by roughly an order of magnitude and C3 changes sign, so the doubled domain does not establish convergence. The L36 near-robot and end-region pressure/velocity proxies also remain strongly transient and nonuniform. No L36-to-L72 run is authorized; no wall-fluid-off Phase 2 replay is authorized because axial domain convergence was not demonstrated.

The required primary classification is therefore `AXIAL_BOUNDARY_SENSITIVITY_UNRESOLVED`.

## Last verified state

At report update, L36 was solved and extracted; the kinematic gate passed, force closure remained below `1.2e-4`, and the case identity retained `restart_read=false`. ODB/SIM/native large files remain untracked and excluded from scoped deliverables.
