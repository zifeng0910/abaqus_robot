# S4 HEAD-FORWARD Low-G Dynamic Screen

These are four new 75 ms Explicit dynamics runs. The only candidate variable is the analytic Gaussian axial gradient: `G=0`, `0.25`, `0.50`, and `1.00 mT`, all with `L=45 mm`. The robot mesh is a true rigid 180 degree flip about fixed `n_routeA`; COM, straight tube, rigid mesh topology, one-wall General Contact, `mu=0.03`, `zeta=0.50`, ReducedHydro coefficients, `B0=10 mT`, initial pose, RouteA gauge, and `dt=1e-7 s` are frozen. These are non-CEL cases.

The fixed orthographic camera is identical in every GIF: horizontal canonical `+s` left-to-right, vertical `n_routeA`, view direction `b_routeA`, and a fixed 20 mm by 5 mm display window. HEAD is the forward/right endpoint and TAIL is the rear/left endpoint. No camera following, auto-fit, camera reversal, or label exchange is used.

| case | G (mT) | delta_s (mm) | fraction v_s > 0 | longest near-zero (ms) | axial reversals | rocking min/max (deg) | max radial drift (mm) | contact count | longest contact (ms) | both-end bridge (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G0 | 0.00 | -0.033092 | 0.1634 | 61.9650 | 0 | -14.656/+14.372 | 0.243232 | 59 | 0.9796 | 0.0000 |
| G0P25 | 0.25 | +0.865311 | 1.0000 | 3.1668 | 0 | -14.914/+14.580 | 0.211527 | 52 | 1.4315 | 0.0000 |
| G0P5 | 0.50 | +1.741598 | 1.0000 | 1.8444 | 0 | -14.873/+14.404 | 0.203483 | 29 | 0.0016 | 0.0000 |
| G1P0 | 1.00 | +3.684497 | 1.0000 | 1.1658 | 0 | -15.037/+14.841 | 0.190813 | 45 | 1.4170 | 0.0000 |
| OLD_S4_REFERENCE_G6 | 6.00 | +24.652332 | 1.0000 | 0.5322 | 0 | -15.545/+14.712 | 0.184479 | 27 | 0.6952 | 1.3000 |

Per-case `v_s(t)` and derived `a_s(t)` are in each case's `*_v_s_a_s_timeseries.csv`. Contact episodes, episode region, minimum gap, peak force, and reopen flag are in each case's `*_contact_episodes.csv`. Both-end bridge is a 20 us screening-resolution exact-gap check; no publication-grade contact analysis or automatic ranking is assigned.

**USER VISUAL SELECTION REQUIRED**
