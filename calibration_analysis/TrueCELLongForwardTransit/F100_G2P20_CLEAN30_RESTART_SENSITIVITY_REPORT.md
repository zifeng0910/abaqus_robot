# F100/G2.20 CLEAN30 restart sensitivity

**CLEAN_RERUN_NOT_REPRODUCIBLE**

CLEAN30 used one uninterrupted 0-30 ms Explicit step, B0=11 mT, G=2.20 mT, f=100 Hz, A_main=14.5 deg, A_cross=2.5 deg, c0=100000 mm/s, a 44x20x20 CEL mesh, scale factor 0.4, no mass scaling, one CPU and double precision. The original field interval of 8.333333 ms and per-increment RP/contact histories were retained.

The 10-ms restart audit demonstrated non-reproducible contact/CEL continuation dynamics. The old C3 delta_s=-0.023527 mm came from a 20-30 ms restart and is retained as `RESTART_CONTAMINATED_C3_RESULT`. The new single-step 0-30 ms run completed, but its 0-20 ms trajectory failed the clean baseline reproduction gate. Its 20-30 ms history is recorded but C3 metrics and physical interpretation are withheld.

The original 20 ms deck wrote restart states at 5/10/15/20 ms. Keeping `number interval=4` in the 30 ms deck moved those writes to 7.5/15/22.5/30 ms. This output-timing difference coincides with the first observed online departure near 5 ms; causality is not established by timing alone.

| Clean 0-20 ms quantity (s: mm; v_s: mm/s; angle: deg; omega: rad/s) | RMS | p99 absolute | Max absolute | First threshold crossing (ms) |
|---|---:|---:|---:|---:|
| s | 0.0274066 | 0.0964296 | 0.104882 | 5.184424 |
| v_s | 10.7948 | 36.4142 | 49.5099 | 5.113057 |
| rocking_angle | 0.690773 | 3.38529 | 3.75979 | 6.039985 |
| omega_rock | 22.3637 | 127.644 | 156.711 | 5.210514 |

| Magnetic / whole General Contact quantity (force: N; torque: N mm) | RMS | p99 absolute | Max absolute |
|---|---:|---:|---:|
| Fmag1 | 6.20827e-08 | 2.71527e-07 | 2.76015e-07 |
| Fmag2 | 3.93196e-09 | 1.7197e-08 | 1.74812e-08 |
| Fmag3 | 1.31948e-08 | 5.77094e-08 | 5.86633e-08 |
| Mmag1 | 2.51199e-05 | 8.48916e-05 | 8.57136e-05 |
| Mmag2 | 8.69215e-05 | 0.000531534 | 0.000579454 |
| Mmag3 | 0.000137663 | 0.000661518 | 0.000732214 |
| Fmag_s | 6.35911e-08 | 2.78124e-07 | 2.82721e-07 |
| CFN1 on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF | 0.000526091 | 0.00207864 | 0.00569775 |
| CFN2 on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF | 0.00176758 | 0.00856975 | 0.0439924 |
| CFN3 on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF | 0.000916042 | 0.00276101 | 0.0190604 |
| CFNM on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF | 0.0018568 | 0.00861678 | 0.0445765 |

Magnetic loads depend on the evolving pose and position; their later differences do not independently establish an excitation mismatch. The General Contact values cover the robot surface's whole contact domain, not an isolated wall-pair force. The metrics JSON records each gate. The requested GIFs show the recorded motion with a non-authoritative label; no clean C3 cycle metrics or physical interpretation were produced after this failed gate.
