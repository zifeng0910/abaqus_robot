# Local Head-Tail Light-Contact Validation

## Decision

**CONTACT_DESTROYS_RECOVERED_ROCKING**

The required three-label taxonomy has no label for the observed fourth outcome: line-like rocking is preserved, but contact occurs at BODY/TAIL rather than alternating HEAD/TAIL. The failure label therefore denotes destruction of the requested HEAD-TAIL contact topology, not destruction of planar rocking.

1. The completed ROCK10 reference genuinely recovered RouteA-style planar rocking: robot range `-10.250390..10.251369 deg`, fitted amplitude `10.004025 deg`, and PCA line-likeness `0.999999999974`.
2. The light-contact robot reached `+14.577303 deg` and `-14.427106 deg`; these are robot angles, not field commands.
3. Light-contact `theta_cross` RMS/max-absolute is `0.0297923/0.0661111 deg`; PCA line-likeness is `0.999991505216`.
4. Exact positive first touch is `+13.848755976 deg`.
5. Exact negative first touch is `-13.848755976 deg`.
6. Static exact geometry says `TAIL` controls positive first touch and `TAIL` controls negative first touch. This asymmetry was retained rather than relabelled to fit the desired motion.
7. The ROCK10 no-drift dynamic gain used for design was `1.000402578470`, with phase lag `-0.001750226 deg`.
8. The target was the more conservative of a 10 um free-pose overtravel and a 0.5 deg cap on each side, then the minimum symmetric field amplitude satisfying both biased, gain-corrected peak inequalities.
9. The selected field amplitude was `14.343111711438 deg`; no other dynamic amplitude was run.
10. Force-resolved HEAD-only/TAIL-only fractions are `0.000000/0.000048`. Alternating HEAD/TAIL events: `False`.
11. There are `14` force-resolved events and `15` separate intervals; longest separated interval is `79.980395 ms`. Every event measurably reopens by more than 0.05 um before the next event/end: `True`.
12. Both-end fraction is `0.000000` and longest bridge is `0.000000 ms`.
13. Line-like rocking preserved: `True`. Minimum exact gap is `-0.003678 um`; deep penetration is `False`.
14. Final classification: **CONTACT_DESTROYS_RECOVERED_ROCKING**.
15. Exactly one next recommendation: **Before any later solve, revise the zero-solve contact target model to include the measured dynamic COM shift and the observed BODY/TAIL event sequence; do not change contact, hydro, or geometry.**

## Contact Event Audit

`CFN` and `CFS` are the Abaqus whole-robot surface history components. Event duration, impulses, and HEAD/TAIL/BODY attribution use native `1e-7 s` contact increments; attribution recomputes the exact faceted-wall gap for every active increment. The full-cycle pose/GIF grid is sampled at 0.1 ms. Relative normal and tangential velocities are rigid-surface kinematic estimates at each region's closest solver node because Abaqus did not write a direct relative-contact-velocity history.

| event | start_ms | end_ms | peak_time_ms | peak_history_index | duration_ms | state | peak_normal_force_N | peak_tangential_force_N | peak_resultant_force_N | normal_impulse_Ns | tangential_impulse_Ns | HEAD_min_gap_um | TAIL_min_gap_um | BODY_min_gap_um | max_closing_normal_velocity_mm_s | max_tangential_velocity_mm_s | separation_before_ms | separation_after_ms | maximum_reopened_gap_um | gap_reopened_after | theta_rock_deg | theta_cross_deg | alpha_B_deg | T_rock_Nmm | COM_radial_um | magnetic_torque_pushes_deeper |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 39.061 | 39.0623 | 39.0614 | 390614 | 0.00140013 | TAIL | 0.0262312 | 0.000682806 | 0.0262392 | 2.05147e-08 | 5.8352e-10 | 44.7672 | -0.00149704 | 34.8239 | 4.2426 | 2.84921 | 39.061 | 3.5094 | 12.3041 | True | 13.8708 | 6.00999e-06 | 13.5045 | -6.92752e-05 | 0.389763 | False |
| 2 | 42.5718 | 42.573 | 42.5721 | 425721 | 0.00129954 | TAIL | 0.0291066 | 0.000758175 | 0.029114 | 2.17859e-08 | 6.38541e-10 | 31.221 | -0.00167148 | 21.2731 | 4.52788 | 4.45809 | 3.5095 | 3.8714 | 17.5934 | True | 14.2889 | -0.0106764 | 13.9543 | -6.22302e-05 | 7.68916 | False |
| 3 | 46.4445 | 46.4459 | 46.4449 | 464449 | 0.00149698 | BODY | 0.043696 | 0.00116999 | 0.0437083 | 3.84955e-08 | 1.0772e-09 | 10.2202 | 12.0532 | -0.00257602 | 6.51643 | 3.05349 | 3.8715 | 2.9872 | 9.35667 | True | 14.5814 | -0.0165729 | 14.2537 | -6.05734e-05 | 24.7261 | False |
| 4 | 49.4332 | 49.4346 | 49.4336 | 494336 | 0.00149698 | BODY | 0.0534072 | 0.00142894 | 0.0534224 | 4.67109e-08 | 1.35153e-09 | 9.50509 | 13.7959 | -0.00314359 | 7.92487 | 5.92237 | 2.9873 | 3.5032 | 12.8539 | True | 14.539 | -0.0519966 | 14.3408 | -3.54022e-05 | 25.5521 | False |
| 5 | 52.9379 | 52.9392 | 52.9382 | 529382 | 0.00140013 | TAIL | 0.0468 | 0.00122201 | 0.0468116 | 3.661e-08 | 1.05306e-09 | 26.8278 | -0.00269365 | 18.5124 | 7.36663 | 4.22005 | 3.5033 | 79.9804 | 228.575 | True | 14.4243 | -0.0525382 | 14.282 | -2.48936e-05 | 9.27072 | False |
| 6 | 132.92 | 132.921 | 132.92 | 1329201 | 0.00150071 | BODY | 0.0477706 | 0.00126999 | 0.0477852 | 4.06589e-08 | 1.1831e-09 | 14.6738 | 89.1104 | -0.00279391 | 6.91714 | 4.46695 | 79.9805 | 1.02181 | 1.09468 | True | -12.1729 | -0.0593202 | -12.3273 | -3.07163e-05 | 60.5748 | True |
| 7 | 133.943 | 133.944 | 133.943 | 1339435 | 0.00150071 | BODY | 0.0258685 | 0.000687162 | 0.0258768 | 2.27571e-08 | 6.44999e-10 | 14.1987 | 82.0339 | -0.00151199 | 4.06028 | 5.70732 | 1.02191 | 0.949596 | 0.470755 | True | -12.3959 | -0.0271718 | -12.5567 | -3.16304e-05 | 57.1386 | True |
| 8 | 134.894 | 134.896 | 134.895 | 1348945 | 0.00150071 | BODY | 0.00888767 | 0.000237195 | 0.00889016 | 7.7037e-09 | 2.21546e-10 | 13.594 | 70.5032 | -0.000509749 | 1.34762 | 5.21155 | 0.949696 | 7.7466 | 18.942 | True | -12.7507 | 0.0368165 | -12.7581 | -1.57559e-06 | 51.653 | True |
| 9 | 142.642 | 142.643 | 142.643 | 1426426 | 0.0013964 | TAIL | 0.0641614 | 0.00166953 | 0.0641788 | 5.10948e-08 | 1.43458e-09 | 44.4519 | -0.00367802 | 33.9234 | 10.8946 | 5.39246 | 7.7467 | 1.323 | 2.54001 | True | -13.9116 | 0.0241028 | -13.9617 | -1.03099e-05 | 1.21537 | True |
| 10 | 143.967 | 143.968 | 143.967 | 1439670 | 0.0013964 | TAIL | 0.0471898 | 0.00122416 | 0.0472032 | 3.71792e-08 | 1.05011e-09 | 46.9644 | -0.00268317 | 36.7556 | 7.87602 | 4.82505 | 1.3231 | 0.778709 | 0.815095 | True | -13.8445 | -0.0275728 | -14.0862 | -4.77871e-05 | 1.05439 | True |
| 11 | 144.747 | 144.748 | 144.747 | 1447470 | 0.00129209 | TAIL | 0.0244778 | 0.0006369 | 0.0244842 | 1.81294e-08 | 5.34336e-10 | 42.6283 | -0.00140422 | 32.821 | 3.75443 | 5.03227 | 0.778809 | 0.793297 | 0.373602 | True | -13.9784 | -0.0575023 | -14.1482 | -3.28993e-05 | 1.71815 | True |
| 12 | 145.541 | 145.543 | 145.542 | 1455417 | 0.0014113 | TAIL | 0.00786388 | 0.000203943 | 0.0078662 | 6.21086e-09 | 1.75907e-10 | 35.2263 | -0.000446312 | 25.9532 | 1.28517 | 5.65827 | 0.793397 | 4.46649 | 13.0274 | True | -14.1988 | -0.0285953 | -14.2026 | -9.74389e-07 | 5.27982 | True |
| 13 | 150.009 | 150.011 | 150.01 | 1500096 | 0.00150071 | BODY | 0.0566933 | 0.00151362 | 0.0567106 | 4.88674e-08 | 1.43084e-09 | 7.56916 | 20.2532 | -0.00333356 | 8.08031 | 4.73885 | 4.46659 | 2.3969 | 4.49457 | True | -14.4296 | 0.00849257 | -14.3431 | 1.51478e-05 | 28.8463 | False |
| 14 | 152.408 | 152.409 | 152.408 | 1524080 | 0.00150071 | BODY | 0.0356234 | 0.000950462 | 0.0356341 | 3.06181e-08 | 8.91822e-10 | 8.10248 | 21.275 | -0.0020975 | 5.14128 | 3.83557 | 2.397 | 47.5909 | 148.681 | True | -14.372 | -0.00685963 | -14.3021 | 1.19869e-05 | 29.2632 | False |

The first force-resolved event occurs at `39.061401 ms`: `theta_rock=+13.870798 deg`, `alpha_B=+13.504503 deg`, `T_rock=-6.92752e-05 N mm`, normal/tangential impulse `2.05147e-08/5.8352e-10 N s`, closest-point closing/tangential speed `4.2426/2.84921 mm/s`, and COM radial shift `0.389763 um`. Magnetic torque pushes toward larger absolute rocking angle: `False`. The event is `TAIL`, then opens to `12.304105 um` before the next event.

## Solver And Identity Gate

- Case status: `SOLVED`; Abaqus normal completion: `True`; fatal-error scan: `False`.
- Final RP/contact/energy history times: `0.200000003`, `0.200000003`, `0.200000003 s`.
- Explicit increment: `1e-07 s`; ODB frames: `401`; socket stderr bytes: `0`.
- Live/vendored magnetic-server hashes match the frozen identity: `True`.
- Field magnitude maximum absolute error from 10 mT: `3.1189e-11 T`.
- Final robot rocking angle at 200 ms: `-0.034859 deg`; zero crossings: `1`.
- COM radial excursion: `110.546478 um`; net axial displacement: `-0.00351942193 mm`.
- Normal/tangential impulses: `4.36367e-07/1.24314e-08 N s`; peak normal/tangential/resultant force: `0.0641614/0.00166953/0.0641788 N`.
- Energy ranges (N mm): ALLKE `0..7.11687e-07`, ALLIE `0..0`, ETOTAL `-1.67951e-20..1.39934e-14`.

The complete fixed-camera GIFs, not scalar metrics alone, are the final topology gate. The legacy RouteA animation is used only for normalized rocking-phase comparison.
