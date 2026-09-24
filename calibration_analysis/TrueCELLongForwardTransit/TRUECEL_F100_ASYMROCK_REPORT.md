# F100 / G2.20 asymmetric rocking screen

## Phase 0: native-history geometric audit

Classification: `WALL_EXPOSURE_GEOMETRICALLY_ASYMMETRIC`.

The existing symmetric F100 run was reconstructed from 164,986 native RP history samples spanning 0–30 ms. Each rigid robot end-region vertex was transformed with the measured RP translation and rotation vector; signed clearance was taken against the pipe polygon's inward half-planes. This is a geometric clearance, not a direct pair-isolated robot-wall force. The available whole-General-Contact resultants are not used as wall force.

| Cycle | HEAD min gap (mm) | HEAD gap < 0 (ms) | HEAD < 5 µm (ms) | HEAD < 10 µm (ms) | TAIL min gap (mm) | TAIL gap < 0 (ms) | TAIL < 5 µm (ms) | TAIL < 10 µm (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | +0.003014 | 0 | 0.1553 | 1.8298 | −0.002851 | 0.1154 | 0.2969 | 0.4641 |
| 2 | +0.004986 | 0 | 0.0045 | 0.2660 | −0.007027 | 0.3539 | 1.1704 | 1.8821 |
| 3 | +0.002897 | 0 | 0.1339 | 0.3350 | −0.004684 | 3.8927 | 6.1641 | 6.9664 |

The stationary geometry's specified touch thresholds were independently reproduced: at −14.387845° the HEAD clearance is approximately 0, and at +13.423314° the TAIL clearance is approximately 0. The static robot at −14.5° also has a negative TAIL half-plane clearance, so these thresholds are not exclusive side-contact classifiers. During free dynamics the closest TAIL approach can occur at negative actual rocking angles (for example −10.18° in C1). Thus the *threshold mapping* is confirmed, but the stronger claim “negative angle approaches only HEAD / positive angle approaches only TAIL” is false once translation and cross-rocking are included. The experiment tests waveform balancing, not a proven contact cause.

## One waveform-only candidate

`TRUECEL_B0P11_G2P20_F100_ASYMROCK_FAST` starts from the same clean t=0 deck and runs 0–30 ms. The Fortran VUAMP remains the free-dynamics magnetic torque/gradient implementation, with phase `mod(36000 t, 360)` and unchanged B0=0.011 T, G=0.0022 T, 2.5° cross component, and gradient table. No robot rotation is prescribed.

The parent table was verified against the analytic field expression to maximum component error < 1×10⁻¹². The sole intended control change is the main field-angle waveform:

`alpha_main(t) = −0.4822655° + 14.0255795° sin(2π·100 t)`.

`B_unit(t) = normalize[−c − tan(alpha_main(t)) n + tan(2.5° cos(2π·100 t)) b]`.

This gives commanded main extrema −14.507845° and +13.543314° with unchanged phase and normalized field magnitude. These are **field-angle targets**, not guaranteed actual rigid-body angle extrema. The spatial gradient rows are copied exactly from the parent. The candidate passed the exact clean-initialization audit (zero robot/fluid structural overlap and zero robot/wall penetration) and Abaqus datacheck.

## Dynamic result

Final classification: `ASYMROCK_RECOIL_FAIL`.

The single fresh-start candidate passed datacheck and Abaqus reported `THE ANALYSIS HAS COMPLETED SUCCESSFULLY`. Its online physical fail-fast ended dynamics at 18.385337 ms (`F100_RECOIL_FAIL_POSITION`): position had fallen >0.010 mm below the running maximum for 5 ms after exceeding 0.020 mm backtrack. The final native RP history sample is at 18.385518 ms. This is **not** a numerical failure and **not** a completed three-cycle run. No second candidate or five-cycle continuation was launched.

| Metric | Symmetric F100/G2.20 | Asymmetric F100/G2.20 |
|---|---:|---:|
| C1 Δs (mm) | +0.092337 | +0.096870 |
| C2 Δs (mm) | +0.143159 (complete) | −0.231405 (10–18.386 ms **partial**) |
| C3 Δs (mm) | −0.023527 | Not run |
| MAX_BACKTRACK (mm) | 0.032364 (30 ms) | 0.267091 (18.386 ms) |
| Mean speed (mm/s) | +7.0656 (30 ms) | −7.3174 (18.386 ms; unequal durations) |
| First negative-velocity cycle | C3 | C2 (11.927 ms) |
| HEAD minimum signed gap, C1 (mm) | +0.003014 | +0.002553 |
| TAIL minimum signed gap, C1 (mm) | −0.002851 | −0.002711 |
| HEAD dwell <10 µm, C1 (ms) | 1.8298 | 0.8097 |
| TAIL dwell <10 µm, C1 (ms) | 0.4641 | 1.2287 |

Only C1 and its wall-exposure rows are duration-matched in the last four comparisons. Full cycle-by-side gap/dwell data, including the partial C2 of the new run, are in the native-history CSVs. The asymmetric candidate's TAIL gap in its partial C2 reached −0.004481 mm and stayed within 10 µm of the wall for 4.6679 ms; HEAD minimum was +0.003983 mm. A brief negative speed alone was not grounds for rejection; the sustained **positional** recoil was.

The candidate peaked at +0.132557 mm at 12.141 ms, then ended at −0.134534 mm. Its C1 mean speed was +9.6870 mm/s; the truncated C2 interval averaged −27.5957 mm/s and ended at −96.4182 mm/s. The phase-synced comparison GIF stops at the new run's available end time. The single-candidate GIF depicts the same available 18.386 ms; a five-cycle GIF is inapplicable.

The field-angle bias was a controlled test of unequal geometric exposure, **not** a direct prescription of the robot angle and not proof that wall recontact causes recoil. In fact, despite a slightly better C1, TAIL near-wall dwell increased in C1 and sustained reverse motion began during C2. Thus this one simple waveform adjustment worsened recoil relative to the symmetric parent. Per the stated hard stop, the follow-on assessment is `SCALAR_AND_SIMPLE_WAVEFORM_TUNING_INSUFFICIENT`; do not run another amplitude, frequency, B0, or G scan. A cycle-to-cycle fluid-state / feedback-control investigation is the next stage, but was not automatically started.
