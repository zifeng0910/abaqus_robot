# FAST precomputed magnetic v_mean calibration

`SELECTED_G = 1.184 mT`

| G (mT) | delta_s (mm) | v_mean (mm/s) | v_final (mm/s) | rocking (deg) | TAIL contacts | HEAD contacts | peak contact (N) | rebound | runtime (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|
| 2.000 | +0.219391 | +13.1634 | +31.311 | -12.81..+13.33 | 10 | 0 | 0.1202 | True | 134.7 |
| 3.000 | +0.382423 | +22.9454 | +50.814 | -12.88..+13.33 | 14 | 0 | 0.1193 | True | 124.6 |
| 4.000 | +0.550597 | +33.0358 | +70.980 | -12.96..+13.33 | 10 | 0 | 0.1198 | True | 130.9 |
| 1.184 | +0.084468 | +5.0681 | +15.138 | -12.93..+13.33 | 12 | 0 | 0.1205 | True | 124.5 |

## Five-cycle verification

- Overall five-cycle mean: `+20.557236 mm/s`
- Cycles 2-5 mean: `+25.433288 mm/s`
- Per-cycle means: `+1.053, +9.094, +20.929, +31.316, +40.394` mm/s
- Five-cycle target gate (4.5-5.5 mm/s): **FAILED**
- Rocking: `-13.0111 to +13.5279 deg`, dominant frequency `119.9981 Hz`
- Contact safety hard gate: `True`; no tumble, no gross penetration, no BOTH-wall bridge
- Rebound diagnostic: `True`; TAIL same-end recontact remains present
- Backend: `PRECOMPUTED_TABLE`; fluid: `FAST_SURROGATE_FLUID`; socket calls: `0`; CEL: `false`

The fitted two-cycle case matched the requested mean (`+5.068059 mm/s`, `10.136x` the background flow), but the five-cycle run continued accelerating. It therefore is the bounded calibration selection, not a validated steady 5 mm/s propulsion case. No additional refinement was run because the calibration limit allowed only one fitted follow-up.
