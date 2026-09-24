# Magnetic-only / no-fluid F100 G2P20 control

**MAGNETIC_ONLY_STABLE_FORWARD**

Fresh t=0 to 50 ms, one uninterrupted Explicit step. The Eulerian fluid domain, fluid material, fluid/contact interactions, and all General Contact were removed. The unchanged wall geometry remains as a passive, fixed reference; contact is absent. Robot geometry, mass, initial pose, Magpylib-derived table, phase convention, B0/G/f/A, and VUAMP load calculation were preserved.

The all-rigid model required DIRECT 1e-7 s time increments because Abaqus could not derive an automatic stable increment without a deformable element. This is a numerical setting change, not a new force or damping term. The source step's bulk-viscosity keyword was left unchanged; no artificial damping was added. Fixed 5-ms restart writes were retained and never read for continuation.

| Cycle | delta_s (mm) | mean v_s (mm/s) | end v_s | minimum v_s | MAX_BACKTRACK (mm) | rocking amplitude (deg) | phase lag (deg) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | +0.127606 | +12.761 | +26.379 | +0.000 | 0.000000 | 21.177 | +7.08 |
| 2 | +0.409383 | +40.938 | +54.737 | +26.379 | 0.000000 | 25.671 | +1.66 |
| 3 | +0.694075 | +69.408 | +83.746 | +54.737 | 0.000000 | 21.931 | -8.30 |
| 4 | +0.981784 | +98.178 | +112.288 | +83.746 | 0.000000 | 21.358 | +5.86 |
| 5 | +1.271742 | +127.174 | +141.407 | +112.288 | 0.000000 | 25.204 | +2.92 |

Rocking amplitude is half the measured angle range per cycle. Phase lag is the first-harmonic angle lag relative to the unchanged 14.5-degree, 100-Hz sine command; positive means measured rocking lags the command. The first cycle includes the unchanged 1-ms magnetic ramp.

| Cycle | Mean magnetic force (N) | Peak magnetic force (N) | Mean magnetic torque (N mm) | Peak magnetic torque (N mm) |
|---:|---:|---:|---:|---:|
| 1 | 2.691e-05 | 2.921e-05 | 1.303e-03 | 2.523e-03 |
| 2 | 2.832e-05 | 2.932e-05 | 1.389e-03 | 2.564e-03 |
| 3 | 2.863e-05 | 2.966e-05 | 1.393e-03 | 2.572e-03 |
| 4 | 2.901e-05 | 2.977e-05 | 1.318e-03 | 2.505e-03 |
| 5 | 2.925e-05 | 3.044e-05 | 1.395e-03 | 2.574e-03 |

Magnetic force and torque are vector magnitudes from the unchanged VUAMP magnetic log. Their three components, RP position and velocity, rocking angle, and angular velocity are in the 10-us CSV.

| Matched comparison | FULL CEL | NO FLUID |
|---|---:|---:|
| Cycle 1 delta_s (mm) | +0.092337 | +0.127606 |
| Cycle 2 delta_s (mm) | +0.139908 | +0.409383 |
| Cycle 3 delta_s (mm) | +0.013300 | +0.694075 |
| Cycle 4 delta_s (mm) | -0.564697 | +0.981784 |
| Cycle 5 delta_s (mm) | N/A (stopped at 40.755 ms) | +1.271742 |
| MAX_BACKTRACK, 0-40 ms (mm) | 0.565250 | 0.000000 |
| Mean speed, 0-40 ms (mm/s) | -7.979 | +55.321 |

The FULL CEL Cycle 5 is incomplete; its 40-40.755 ms partial delta_s was -0.072202 mm and is excluded from the five-cycle comparison. The no-fluid 0-50 ms mean speed was +69.692 mm/s.

The magnetic-only trajectory remains directionally stable under the stated 0.010-mm within-cycle backtrack threshold. The full coupled trajectory reverses in Cycle 4. This is consistent with the removed fluid/contact physics contributing to recoil. The control also requires a fixed Explicit time increment, so this comparison alone does not uniquely attribute the difference to fluid forces, wall contact, or their interaction.

Numerical checks: Abaqus completed 50 ms; peak |ETOTAL|=1.468e-12 N mm versus peak ALLKE=5.964e-04 N mm. The four inherited distorted tetrahedra are in the fixed passive PIPE_SOLID and were also present in the FULL CEL source. RP, magnetic loads, and energy histories are finite. No restart continuation or feedback was used.
