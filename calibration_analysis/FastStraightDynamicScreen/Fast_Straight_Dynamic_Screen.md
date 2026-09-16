# Fast Straight Dynamic Screen

All four cases are new short dynamics runs. They use the same non-CEL rigid model, one-wall General Contact, ReducedHydro, `G=6 mT`, `L=45 mm`, `B0=10 mT`, `dt=1e-7 s`, initial pose, and RouteA gauge.

The GIF camera is fixed and identical: horizontal canonical `+s`, vertical `n_routeA`, view direction `b_routeA`, 50 mm displayed tube length, 1 ms physical frame spacing, and 20 fps playback. Shorter 75 ms cases hold their last frame while S1 completes its final 25 ms in the four-way comparison.

| case | delta_s (mm) | fraction v_s > 0 | longest near-zero v_s (ms) | axial reversals | rocking min/max (deg) | max radial COM drift (mm) | contact count | longest contact (ms) | both-end bridge (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| S1 | +44.255487 | 1.0000 | 0.5322 | 0 | -12.928 / +12.955 | 0.013170 | 0 | 0.0000 | 0.0000 |
| S2 | +25.105890 | 1.0000 | 0.5322 | 0 | -13.298 / +13.257 | 0.005679 | 0 | 0.0000 | 0.0000 |
| S3 | +24.633356 | 1.0000 | 0.5322 | 0 | -14.987 / +14.974 | 0.574613 | 25 | 0.8637 | 6.1000 |
| S4 | +24.652332 | 1.0000 | 0.5322 | 0 | -15.545 / +14.712 | 0.184479 | 27 | 0.6952 | 1.3000 |

Near-zero means `|v_s| <= 1 mm/s`; contact uses whole-robot resultant force above `1e-8 N`. Both-end bridge is a screening-resolution exact-gap check at 0.1 ms, not a publication-grade contact audit.

No automatic ranking is assigned.

**USER VISUAL SELECTION REQUIRED**
