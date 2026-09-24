# F100/G2P20 3x rotational damping, full CEL

**CROT3_FULLCEL50_RECOIL_CONFIRMED.** Completed one Explicit step from t=0 to 50 ms; restart snapshots were written at 5-ms intervals and never read.

| Cycle | delta s (mm) | mean v (mm/s) | end v (mm/s) | minimum v (mm/s) | max backtrack (mm) | rocking amplitude (deg) | phase lag (deg) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | +0.119811 | +11.981 | +23.082 | +0.000 | 0.000000 | 12.718 | +32.02 |
| 2 | +0.209000 | +20.900 | -1.613 | -1.658 | 0.000529 | 12.254 | +37.49 |
| 3 | -0.062934 | -6.293 | -9.683 | -14.543 | 0.063059 | 11.917 | +38.38 |
| 4 | -0.161372 | -16.137 | -18.334 | -23.466 | 0.161372 | 11.258 | +40.62 |
| 5 | -0.174245 | -17.424 | -15.154 | -24.251 | 0.174245 | 11.426 | +39.88 |

Final net displacement: -0.069741 mm; whole-run max backtrack: 0.399079 mm.
Geometric wall overlap: 24 sampled episodes, minimum signed gap -0.005928 mm. This is a sampled geometry indicator, not an isolated contact-force count.
Whole-robot General Contact CFNM peak: 0.0210976 N; this signal includes fluid contact.

The clean FULL CEL baseline had C3 +0.013300 mm and C4 -0.564697 mm before user termination at 40.75 ms. Magnetic-only 3x had nearly constant 14.728-deg amplitude and 34.17-deg lag in C2-C5, but its passive-wall audit also found wall overlap. The coupled 3x case therefore does not establish a physical wall-free rocking attractor.

Next bounded test: reduce only rotational damping to 1x (3.5006005000754733e-6 N mm s), whose magnetic-only screen showed stable C2-C5 amplitude and phase; hold all other model terms fixed.
