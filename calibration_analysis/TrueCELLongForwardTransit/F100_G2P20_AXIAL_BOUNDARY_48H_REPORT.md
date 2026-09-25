# F100_G2P20 Axial Boundary Campaign

## Primary classification

`AXIAL_BOUNDARY_SENSITIVE_L36_REQUIRED_WALL_EXTENSION_PENDING`

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

## Boundary result

The inferred robot-CEL force is defined as whole-robot General Contact minus the direct robot-wall pair. Relative to L12, axial impulse changes were `D_CEL = 0.397, 0.913, 1.039, 0.023` for C1-C4 and `D_wall = 0.105, 0.130, 0.239, 0.158`. C2 and C3 therefore remain materially domain-sensitive. Near-robot and end-region fluid velocity and EVF proxies also differ between L12 and L24.

The six native field frames contain velocity, EVF, and water stress only. Pressure is reported as the explicitly labelled mean-normal-stress proxy `-trace(S_water)/3`; no native PRESS or face mass-flux history exists.

## Next decision

L36 is justified by the L24 sensitivity result, but it must first be prepared with a validated physical wall-helper extension to the L36 axial span. No L36 input currently exists, so no additional expensive job was launched. The current verified result is sufficient to reject an axial-domain-insensitive classification.
