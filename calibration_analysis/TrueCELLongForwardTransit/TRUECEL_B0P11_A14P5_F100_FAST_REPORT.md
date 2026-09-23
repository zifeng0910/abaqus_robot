# TRUE-CEL B0P11 A14P5 — single 100-Hz candidate

Classification: **F100_RECOIL_FAIL**

F120 torque-authority check: **ROCKING_CONTROL_REMAINS_COHERENT**. The measured phase-matched rocking amplitudes were 12.985° and 13.121°, fundamental shift −2.42°, RMS angle difference 2.272°. The p95 absolute instantaneous omega_rock difference was 220.43 rad/s (not small), so angular-rate variability remains visible even though the rocking amplitude and phase did not collapse.

F100: fresh start at t=0; only the forcing frequency changed from 120 to 100 Hz. F100 end time 20.000000 ms; complete cycles 2; mean speed +7.3221 mm/s; MAX_BACKTRACK 0.016973 mm; stage-1 gate FAIL/NOT PASSED.

| frequency (Hz) | Cycle 1 Δs (mm) | Cycle 2 Δs (mm) | Cycle 3 Δs (mm) | mean speed (mm/s) | MAX_BACKTRACK (mm) | Cycle 2 end v_s (mm/s) |
|---:|---:|---:|---:|---:|---:|---:|
| 120 | 0.06094774095508182 | -0.027793991839170307 | None | 1.9886538 | 0.07236198 | -28.5378 |
| 100 | 0.0850081431056382 | 0.06143387991563776 | None | 7.322101314725516 | 0.016973250847166277 | -3.463073331070662 |

Central question — does 100 Hz remove the second-cycle collapse? **It removes the cycle-integrated collapse: Cycle 2 Δs becomes positive, but sustained second-cycle recoil remains, so the two-cycle gate fails and no third cycle is run.**

Cycle metrics are in `TRUECEL_B0P11_A14P5_F100_FAST_CYCLE_SUMMARY.csv`; machine-readable results are in `TRUECEL_B0P11_A14P5_F100_FAST_METRICS.json`. Motion uses canonical +s left-to-right with tail left/head right. GIF: `TRUECEL_B0P11_A14P5_F100_FAST.gif`.
