# TRUE-CEL B0P11 A14P5 FAST full-transit screen

Classification: **FORWARD_PROGRESS_STALL**

- reached_finish: NO
- travel_distance_mm: 7.566533592
- achieved_net_displacement_mm: 0.102252568
- transit_time_ms: 13.471364
- termination_time_ms: 13.471186
- rolling_2ms_gain_at_termination_mm: 0.004839104
- mean_speed_mm_s: 7.590365
- minimum_v_s_mm_s: -10.406850
- max_backward_excursion_mm: 0.003253639
- longest_negative_velocity_duration_ms: 1.009882
- full_rocking_cycles_completed: 1
- online_event: FORWARD_PROGRESS_STALL

This is a `REDUCED_SOUND_SPEED_COARSE_CEL_FAST_SCREEN`, not final quantitative FSI validation.

## Frozen parent configuration

- parent: `TRUECEL_NO_RECOIL_B0P11_2CYCLES`
- B0: 11 mT
- G: 2 mT
- frequency: 120 Hz
- A_main: 14.8 -> 14.5 deg (only physics change)
- A_cross: 2.5 deg
- c0: 100000 mm/s
- CEL mesh: 44 x 20 x 20
- contact: parent one-wall General Contact unchanged; mu=0.03; zeta=1.0
- Explicit: automatic stable increment, scale factor 0.4, no mass scaling
- s_start / s_finish: -3.200000000 / 4.366533592 mm

## Cycle-by-cycle

| cycle | delta_s (mm) | mean v_s (mm/s) | end v_s | min v_s | max back (mm) |
|---:|---:|---:|---:|---:|---:|
| 1 | +0.060948 | +7.3137 | +13.0499 | +0.0000 | 0.000000 |
