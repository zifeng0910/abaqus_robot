# F100/G2P20 rotational damping, FULL CEL

**CROT03_EARLY_RECOIL_FAIL / USER_STOPPED_AFTER_DECISIVE_RECOIL**. Branch: **ROTATIONAL_DAMPING_DOES_NOT_REMOVE_COUPLED_RECOIL**. Energy: **ENERGY_BEHAVIOR_NUMERICALLY_SUSPECT**.

Abaqus terminated on user request after completed C3 net recoil -0.417724 mm and whole-run online backtrack 0.706535 mm, above the predeclared 0.010 mm gate. This is not a completed 50-ms solve or a spontaneous solver failure. ODB RP history ends at 37.146706 ms; the online motion file ends at 37.256012 ms. C4 is partial and C5 unavailable.

0x was separately stopped by the user at 40.75 ms. The 1x and 3x cases completed 50 ms. All cases began at t=0 with one Explicit step and no restart continuation.

## 0.3x cycle evidence

| Cycle | delta s (mm) | mean v (mm/s) | end v (mm/s) | min v (mm/s) | max backtrack (mm) | fundamental amp (deg) | half range (deg) | phase lag (deg) | theta offset (deg) | residual / amp | omega peak (rad/s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C1 | +0.093665 | +9.367 | +19.641 | +0.000 | 0.000000 | 13.169 | 12.713 | +18.98 | -0.49 | 0.124 | 244.44 |
| C2 | +0.059069 | +5.907 | -5.705 | -8.011 | 0.021874 | 14.452 | 13.057 | +16.04 | +0.19 | 0.136 | 272.67 |
| C3 | -0.417724 | -41.772 | -69.632 | -69.632 | 0.417724 | 14.214 | 12.649 | +16.15 | +0.17 | 0.106 | 246.25 |
| C4 partial (37.147 ms) | -0.266681 | -37.315 | -1.416 | -138.886 | 0.266681 | NA | 7.964 | NA | NA | NA | 344.54 |
| C5 unavailable | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

## Four-level comparison

| Metric | 0x | 0.3x | 1x | 3x |
|---|---:|---:|---:|---:|
| C1 delta s (mm) | +0.092337 | +0.093665 | +0.108167 | +0.119811 |
| C2 delta s (mm) | +0.139908 | +0.059069 | +0.129046 | +0.209000 |
| C3 delta s (mm) | +0.013300 | -0.417724 | -0.049535 | -0.062934 |
| C4 delta s (mm) | -0.564697 | NA | -0.021883 | -0.161372 |
| C5 delta s (mm) | NA | NA | -0.040920 | -0.174245 |
| Max backtrack (mm) | 0.637453 | 0.706279 | 0.112509 | 0.399079 |
| Rocking amplitude spread C2-last complete (%) | 4.09 | 1.66 | 3.67 | 8.50 |
| Phase spread C2-last complete (deg) | 4.55 | 0.11 | 0.83 | 3.13 |
| Minimum wall gap (mm) | -0.007441 | -0.023590 | -0.006616 | -0.005928 |
| Negative-gap episodes | 41 | 24 | 36 | 24 |
| Negative-gap duration (ms) | 11.25 | 8.95 | 7.25 | 10.90 |
| ETOTAL at last available time (N mm) | -9.9887 | -8.2201 | -5.1574 | -0.10662 |
| ALLPW at last available time (N mm) | 10.453 | 10.273 | 5.9956 | 0.15599 |
| Last available time (ms) | 40.75 | 37.15 | 50.00 | 50.00 |

0.3x rocking stability over C2-C5 cannot be established from the stopped run. C2-C3 coherence is provisional and does not explain away axial recoil. Negative wall gap is a 0.05-ms sampled geometric overlap; whole General Contact includes robot-fluid interaction and is not a pair-isolated robot-wall force. Cycle-wise ETOTAL, ALLPW, ALLFD, ALLIE, ALLKE and minimum stable timestep are in the accompanying CSV. Partial case endpoints are not equal-duration comparisons.

Energy is marked suspect because the partial 0.3x |ETOTAL| or ALLPW is within one order of magnitude of the completed 1x case. This flags large progressive drift, without identifying its cause. The 0x energy is available only through 40.75 ms.

Rotational-damping tuning stops here. The next investigation is a separate no-fluid control with real robot-wall contact isolated; it was not started in this branch.
