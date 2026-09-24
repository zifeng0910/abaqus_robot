# F100 axial-bias offline screen

Classification: **OPEN_LOOP_FIRST_ORDER_G_ESTIMATE**

The actual VUAMP source computes `grad = G*grad_unit` and `F_i = ramp*sum_j(moment_j*grad_ji)`. The table tensor is `grad_unit = d(s_eff)*(c outer c)`, therefore at fixed state:

`F_mag,s = ramp*G*d(s_eff)*(moment dot c)`

This is exactly linear in G. G is the amplitude of the Gaussian gradient profile in tesla, not a spatially uniform mT/mm gradient. B0 enters the torque term and remains unchanged.

The measured baseline reaches its last maximum at 16.029561 ms. Backtracking then continues to 20 ms, where the maximum deficit is 0.016973 mm. No forward recovery occurs before the run ends. The negative-v interval is 16.029743–20.000000 ms.

| G (mT) | predicted MAX_BACKTRACK (mm) | predicted min v_s | predicted Cycle-2 end v_s | predicted Cycle-2 delta_s (mm) |
|---:|---:|---:|---:|---:|
| 2.05 | 0.012078 | -6.473 | -2.078 | +0.071733 |
| 2.10 | 0.007533 | -5.194 | -0.694 | +0.082032 |
| 2.15 | 0.003389 | -3.915 | +0.691 | +0.092331 |
| 2.20 | 0.000897 | -2.636 | +2.076 | +0.102630 |
| 2.25 | 0.000206 | -1.358 | +3.460 | +0.112929 |
| 2.30 | 0.000001 | -0.102 | +4.845 | +0.123228 |
| 2.35 | 0.000000 | +0.000 | +6.230 | +0.133527 |
| 2.40 | 0.000000 | +0.000 | +7.614 | +0.143826 |

The lowest predicted strong-pass value is G=2.15 mT. Applying the allowed one-step robustness margin selects **G=2.20 mT** as the only dynamics candidate.

Important limitation: non-magnetic loads and the baseline trajectory are frozen. The pipe-surface force is exported only as a wall-reaction proxy because the existing ODB does not contain pair-isolated robot-wall CFT; it also includes pipe-fluid contact.
