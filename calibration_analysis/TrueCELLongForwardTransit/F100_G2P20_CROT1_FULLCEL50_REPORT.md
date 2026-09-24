# F100/G2P20 1x rotational damping, full CEL

**CROT1_FULLCEL50_RECOIL_CONFIRMED.** Completed one Explicit step from t=0 to 50 ms; restart snapshots were written at 5-ms intervals and never read.

| Cycle | delta s (mm) | mean v (mm/s) | end v (mm/s) | minimum v (mm/s) | max backtrack (mm) | rocking amplitude (deg) | phase lag (deg) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | +0.108167 | +10.817 | +25.040 | +0.000 | 0.000000 | 13.446 | +21.75 |
| 2 | +0.129046 | +12.905 | +0.532 | -0.424 | 0.000049 | 14.369 | +19.74 |
| 3 | -0.049535 | -4.954 | -4.745 | -13.890 | 0.049706 | 14.281 | +20.16 |
| 4 | -0.021883 | -2.188 | +2.219 | -7.249 | 0.022332 | 13.850 | +20.34 |
| 5 | -0.040920 | -4.092 | -9.424 | -11.670 | 0.042556 | 14.196 | +20.57 |

Final net displacement: +0.124874 mm; whole-run max backtrack: 0.112509 mm.
Geometric wall overlap: 36 sampled episodes, minimum signed gap -0.006616 mm. This is a sampled geometry indicator, not an isolated contact-force count.
Whole-robot General Contact CFNM peak: 0.0339733 N; this signal includes fluid contact.
End energies (N mm): ETOTAL=-5.15743, ALLPW=5.99564, ALLFD=1.59816, ALLIE=-0.982388, ALLKE=0.0285719. The large pressure-work and energy drift limit physical interpretation even though Abaqus completed.

The clean FULL CEL baseline had C3 +0.013300 mm and C4 -0.564697 mm before user termination at 40.75 ms. The magnetic-only damping screen showed more regular rocking without fluid, but passive-wall overlap was still present. This coupled case does not establish a physical wall-free rocking attractor.
