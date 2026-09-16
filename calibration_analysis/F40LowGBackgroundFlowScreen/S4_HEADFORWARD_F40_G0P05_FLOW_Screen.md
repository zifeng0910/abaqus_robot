# F40 Low-G Background-Flow Screen

This is exactly one short straight-tube, non-CEL Explicit dynamics run. It freezes the head-forward S4 rigid robot, tube, RouteA gauge, elliptical rocking (`A_main=14.343111711438091 deg`, `A_cross=2.5 deg`), `B0=10 mT`, General Contact, `mu=0.03`, `zeta=0.50`, and ReducedHydro coefficients. The only new physics is prescribed background flow through `V_rel=V_robot-U_flow*c_hat`; this is not full FSI and does not modify the fluid field.

`U_flow=+10 mm/s` along canonical `+s` was used because no authoritative straight-tube flow value was found in the repository. It is `DIAGNOSTIC_SCREENING_FLOW_ONLY`, not a final experimental value. The magnetic gradient is `G=0.05 mT`, `L=45 mm`, and the initial `F_gradient dot c_hat` gate passed. The `dt=2e-7 s` low-precision gate and dynamic run passed; no fallback was needed.

| metric | value |
| --- | ---: |
| actual rocking frequency (Hz) | 40.039090 |
| field-to-robot phase lag (deg) | -1.024200 |
| rocking min/max (deg) | -13.8308 / +14.0451 |
| delta_s (mm) | +0.074884 |
| fraction v_s > 0 | 0.999995 |
| longest near-zero v_s (ms) | 5.675200 |
| axial reversals | 0 |
| final v_s (mm/s) | +3.197505 |
| max radial COM displacement (mm) | 0.126313 |
| HEAD contact episodes | 12 |
| TAIL contact episodes | 24 |
| BODY contact episodes | 13 |
| longest contact (ms) | 0.534200 |
| BOTH bridge duration (ms) | 0.000000 |
| magnetic forward force range (uN) | +0.000000 / +0.693631 |
| hydrodynamic axial force range (uN) | -0.051369 / +0.086243 |
| background flow (mm/s) | +10.000000 |

Flags: persistent bridge=`False`, long wall sliding=`False`, gross penetration=`False`, tumble=`False`, contact chatter=`False`.

The force/velocity plot is `S4_HEADFORWARD_F40_G0P05_FLOW_force_velocity.png`. The slow contact-inspection GIF is `S4_HEADFORWARD_F40_G0P05_FLOW_SLOW.gif`; the faster playback is `S4_HEADFORWARD_F40_G0P05_FLOW_FAST.mp4` when an ffmpeg encoder is available (`fast_mp4_status=NOT_RENDERED_FFMPEG_UNAVAILABLE`). This result is a coarse visual precursor to later true fluid simulation, not a publication-grade FSI result.

**USER VISUAL SELECTION REQUIRED**
