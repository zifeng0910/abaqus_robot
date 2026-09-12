# Wobble30Hz Wall-ON 8.333 ms rotational-load report

## Decision
B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL

## Runtime identity
Completed at 8.333 ms, 168 ODB frames. G=6 mT=0.006 T, L=45 mm, B0=10 mT, commanded f=30 Hz, cone30, Bias40, phase248, sense+1, chi0, FORWARD, Socket update every Explicit increment. No timer/additional Job.

## First 3 ms regression
See wallon_8p333_first3ms_regression.csv. Maximum sampled dB=0.872388 mT, dF=0.879544 uN, dT=0.792835 uN mm, dRP=63.6513 um. Old/new records agree at sampled cadence; dynamics are sampled at 50 us.

## Strict orientation and omega
Formal omega is finite-rotation log of R(t+dt)R(t)^T, not gradient(UR,t). Direct VR was absent in extracted fields. Fixed directed HEAD->TAIL axis reconstruction RMS max=6.92439e-07 mm. Inertia is a 10 mg, 1.22 mm x 2.81 mm cylindrical approximation because no independent rigid-inertia history was exported.

## Phase/windows
Commanded phase remains exactly 30 Hz; local B projection phase is geometric and may differ as the Frenet basis moves. Late W4 contains 64 samples over 3.159 ms: robot phase 5.042 Hz (R=0.168), wobble RMS 128.284 rad/s, Tmag RMS ratio 1.296.

## Momentum/impulse/work
H is kg m2/s. Jmag is N mm s (multiply by 1e-3 for N m s). Jnonmag_resolved=DeltaH-Jmag is an accounting identity, not independent wall torque. Krot/Wmag are in CSVs and Wmag uses 1 N mm=1e-3 J.

## Contact/translation
Exact gap uses authoritative Pipe_WALL_HELPER R3D4 triangles; aggregate CPRESS is kept separate because it may include CEL-fluid contact. Minimum exact gap -2.318 um; aggregate CPRESS maximum 6.524 MPa. Canonical translation/Vt are recorded, not optimized.

## Mechanism and exactly one next step
Classification: B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL. Do not increase B while Tmag remains 1.296x pre-impact. If B, next step is zero-cost elliptical-polarization/contact-load-shaping audit; if C, extend same configuration to 16.667 ms; if A, extend same configuration for stability; if D, only then discuss B; if E, repair output/reconstruction.
\n## Integrated angular balance\nThe window table now contains DeltaH, integrated magnetic angular impulse Jmag, resolved Jnonmag=DeltaH-Jmag, DeltaKrot, Wmag and Wnonmag. Jnonmag is an accounting identity, not an independent contact-torque measurement.\n - W0_preimpact: DeltaH=4.969e-09 kg m2/s, |Jmag|=6.954e-06 Nmm s, |Jnonmag|=1.989e-09 kg m2/s, DeltaKrot=1.646e-06 J, Wmag=2.345e-06 J, Wnonmag=-6.992e-07 J.\n - W1_contact: DeltaH=6.535e-09 kg m2/s, |Jmag|=4.944e-06 Nmm s, |Jnonmag|=1.126e-08 kg m2/s, DeltaKrot=-1.57e-07 J, Wmag=3.088e-07 J, Wnonmag=-4.659e-07 J.\n - W2_immediate_postimpact: DeltaH=4.601e-09 kg m2/s, |Jmag|=1.553e-05 Nmm s, |Jnonmag|=1.257e-08 kg m2/s, DeltaKrot=-8.297e-08 J, Wmag=1.054e-06 J, Wnonmag=-1.137e-06 J.\n - W3_intermediate: DeltaH=2.722e-09 kg m2/s, |Jmag|=1.906e-05 Nmm s, |Jnonmag|=1.987e-08 kg m2/s, DeltaKrot=-3.77e-07 J, Wmag=9.505e-07 J, Wnonmag=-1.327e-06 J.\n - W4_late_postimpact: DeltaH=1.055e-09 kg m2/s, |Jmag|=3.086e-05 Nmm s, |Jnonmag|=3.167e-08 kg m2/s, DeltaKrot=-2.442e-07 J, Wmag=5.152e-07 J, Wnonmag=-7.595e-07 J.\n

## Cone25 decision addendum

The authorized cone25 run completed at 8.333 ms and was post-processed with the same 50 µs telemetry/ODB cadence and exact SmoothWall114 triangle audit. The comparison is in `cone25_vs_cone30_summary.csv`.

- Cone30 baseline robot arc change: **-2.0000 mm**; cone25: **-2.1000 mm** (both negative, so the motion remains reverse in the canonical arc convention).
- Mean tangential speed proxy: cone30 0.0039 mm/s; cone25 0.0037 mm/s.
- True-axis UR ranges (cone25): UR1 2.578 rad, UR2 0.726 rad, UR3 1.628 rad.
- Exact minimum gap: cone30 -2.059 µm; cone25 -2.318 µm. Both remain below zero; cone25 is not a no-penetration pass.
- Aggregate CPRESS: cone30 6.868 MPa; cone25 6.524 MPa.
- Cone25 late post-impact window: phase-equivalent rate 5.042 Hz, R=0.168, wobble RMS 128.284 rad/s; classification remains **B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL**.

Conclusion: reducing cone30 to cone25 did not remove the geometric wedge or restore forward motion. It slightly increased the reverse arc excursion and made the exact minimum gap more negative, although CPRESS decreased modestly. The zero-cost map's cylinder sanity limit (~23°) explains why 25° can still enter the constrained region once translation and transverse offset are included. No further cone-angle or contact tuning is authorized by this stage; the next action should be selected from the mechanism classification, not from a claim that cone25 solved the jam.
