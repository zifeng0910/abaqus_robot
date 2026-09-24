# F100/G2P20 CLEAN50 user-stopped audit

**CLEAN50_RECOIL_CONFIRMED**, based on the completed clean Cycle 4.

Abaqus was terminated on user request at RP t=40.754825 ms. Cycles 1-4 are complete; Cycle 5 is partial. This is not a completed 50 ms solve and is not classified as a spontaneous numerical failure.

One Explicit step from t=0; no restart read or continuation. Fixed restart writes at [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0] ms were never used for continuation. All model parameters and input differences are recorded in the setup JSON.

| Segment | Time (ms) | delta_s (mm) | mean v_s (mm/s) | end v_s (mm/s) | min v_s (mm/s) | max backtrack (mm) |
|---|---:|---:|---:|---:|---:|---:|
| Cycle 1 | 0.000-10.000 | +0.092337 | +9.234 | +15.339 | +0.000 | 0.000000 |
| Cycle 2 | 10.000-20.000 | +0.139908 | +13.991 | +4.316 | +3.147 | 0.000000 |
| Cycle 3 | 20.000-30.000 | +0.013300 | +1.330 | +2.458 | -7.333 | 0.005306 |
| Cycle 4 | 30.000-40.000 | -0.564697 | -56.470 | -93.718 | -94.313 | 0.565250 |
| Cycle 5 partial | 40.000-40.755 | -0.072202 | -95.654 | -97.937 | -99.345 | 0.072202 |

The old restart-derived C3 net recoil (-0.023527 mm) did not recur: clean C3 delta_s=+0.013300 mm. Clean C4 instead had delta_s=-0.564697 mm and max backtrack=0.565250 mm.

The 0-20 ms path diverged after the 15 ms write. Maximum baseline position difference: 0.00325071 mm; p99 velocity difference: 4.38161 mm/s. This is recorded because restart-write cadence affects the numerical path; it does not invalidate the uninterrupted CLEAN50 observation.

Minimum sampled stable increment: 6.72265e-08 s. The solver reported only the requested external termination at the end.
