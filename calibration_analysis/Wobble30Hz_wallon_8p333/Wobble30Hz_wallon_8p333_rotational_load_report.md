# Wobble30Hz Wall-ON 8.333 ms rotational-load report

## Decision
E_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS

## Runtime identity
Completed at 8.333 ms, 168 ODB frames. G=6 mT=0.006 T, L=45 mm, B0=10 mT, commanded f=30 Hz, cone30, Bias40, phase248, sense+1, chi0, FORWARD, Socket update every Explicit increment. No timer/additional Job.

## First 3 ms regression
See wallon_8p333_first3ms_regression.csv. Maximum sampled dB=0.0012921 mT, dF=0.00721493 uN, dT=0.00597947 uN mm, dRP=0.313836 um. Old/new records agree at sampled cadence; dynamics are sampled at 50 us.

## Strict orientation and omega
Formal omega is finite-rotation log of R(t+dt)R(t)^T, not gradient(UR,t). Direct VR1–VR3 is present and was extracted; its comparison with the relative-rotation reconstruction is reported below. Fixed directed HEAD->TAIL axis reconstruction RMS max=6.76374e-07 mm. Inertia is a 10 mg, 1.22 mm x 2.81 mm cylindrical approximation because no independent rigid-inertia history was exported.

## Phase/windows
Commanded phase remains exactly 30 Hz; local B projection phase is geometric and may differ as the Frenet basis moves. Late W4 contains 64 samples over 3.158 ms: robot phase 8.566 Hz (R=0.286), wobble RMS 90.6267 rad/s, Tmag RMS ratio 1.248.

## Momentum/impulse/work
H is kg m2/s. Jmag is N mm s (multiply by 1e-3 for N m s). Jnonmag_resolved=DeltaH-Jmag is an accounting identity, not independent wall torque. Krot/Wmag are in CSVs and Wmag uses 1 N mm=1e-3 J.

## Contact/translation
Exact gap uses authoritative Pipe_WALL_HELPER R3D4 triangles; aggregate CPRESS is kept separate because it may include CEL-fluid contact. Minimum exact gap -2.059 um; aggregate CPRESS maximum 6.868 MPa. Canonical translation/Vt are recorded, not optimized.

## Mechanism and exactly one next step
Classification: B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL. Do not increase B while Tmag remains 1.248x pre-impact. If B, next step is zero-cost elliptical-polarization/contact-load-shaping audit; if C, extend same configuration to 16.667 ms; if A, extend same configuration for stability; if D, only then discuss B; if E, repair output/reconstruction.
\n## Integrated angular balance\nThe window table now contains DeltaH, integrated magnetic angular impulse Jmag, resolved Jnonmag=DeltaH-Jmag, DeltaKrot, Wmag and Wnonmag. Jnonmag is an accounting identity, not an independent contact-torque measurement.\n - W0_preimpact: DeltaH=4.895e-09 kg m2/s, |Jmag|=6.871e-06 Nmm s, |Jnonmag|=1.981e-09 kg m2/s, DeltaKrot=1.597e-06 J, Wmag=2.286e-06 J, Wnonmag=-6.888e-07 J.\n - W1_contact: DeltaH=6.414e-09 kg m2/s, |Jmag|=4.8e-06 Nmm s, |Jnonmag|=1.101e-08 kg m2/s, DeltaKrot=-1.38e-07 J, Wmag=2.551e-07 J, Wnonmag=-3.93e-07 J.\n - W2_immediate_postimpact: DeltaH=3.609e-09 kg m2/s, |Jmag|=1.486e-05 Nmm s, |Jnonmag|=1.472e-08 kg m2/s, DeltaKrot=-3.12e-07 J, Wmag=1.19e-06 J, Wnonmag=-1.502e-06 J.\n - W3_intermediate: DeltaH=1.861e-09 kg m2/s, |Jmag|=1.815e-05 Nmm s, |Jnonmag|=1.787e-08 kg m2/s, DeltaKrot=-3.406e-07 J, Wmag=7.924e-07 J, Wnonmag=-1.133e-06 J.\n - W4_late_postimpact: DeltaH=1.051e-09 kg m2/s, |Jmag|=2.932e-05 Nmm s, |Jnonmag|=2.926e-08 kg m2/s, DeltaKrot=-4.709e-08 J, Wmag=5.557e-07 J, Wnonmag=-6.028e-07 J.\n\n## Integrated angular balance\nThe window table now contains DeltaH, integrated magnetic angular impulse Jmag, resolved Jnonmag=DeltaH-Jmag, DeltaKrot, Wmag and Wnonmag. Jnonmag is an accounting identity, not an independent contact-torque measurement.\n - W0_preimpact: DeltaH=4.895e-09 kg m2/s, |Jmag|=6.871e-06 Nmm s, |Jnonmag|=1.981e-09 kg m2/s, DeltaKrot=1.597e-06 J, Wmag=2.286e-06 J, Wnonmag=-6.888e-07 J.\n - W1_contact: DeltaH=6.414e-09 kg m2/s, |Jmag|=4.8e-06 Nmm s, |Jnonmag|=1.101e-08 kg m2/s, DeltaKrot=-1.38e-07 J, Wmag=2.551e-07 J, Wnonmag=-3.93e-07 J.\n - W2_immediate_postimpact: DeltaH=3.609e-09 kg m2/s, |Jmag|=1.486e-05 Nmm s, |Jnonmag|=1.472e-08 kg m2/s, DeltaKrot=-3.12e-07 J, Wmag=1.19e-06 J, Wnonmag=-1.502e-06 J.\n - W3_intermediate: DeltaH=1.861e-09 kg m2/s, |Jmag|=1.815e-05 Nmm s, |Jnonmag|=1.787e-08 kg m2/s, DeltaKrot=-3.406e-07 J, Wmag=7.924e-07 J, Wnonmag=-1.133e-06 J.\n - W4_late_postimpact: DeltaH=1.051e-09 kg m2/s, |Jmag|=2.932e-05 Nmm s, |Jnonmag|=2.926e-08 kg m2/s, DeltaKrot=-4.709e-08 J, Wmag=5.557e-07 J, Wnonmag=-6.028e-07 J.\n
## Direct VR cross-check
ODB contains VR1-VR3. Relative-rotation-log omega was cross-checked against direct VR in wallon_8p333_angular_velocity.csv: RMS absolute difference 51.392 rad/s, peak relative difference 0.992.


## Reconstruction gate correction
The ODB does contain direct VR1-VR3. Relative-rotation-log and direct VR differ by RMS 51.392 rad/s and peak relative difference 0.992, exceeding the 10% consistency gate. Therefore the formal mechanism classification is conservatively E (ambiguous), despite the provisional window ratios suggesting nonmagnetic loading. The next action is to resolve the angular-velocity convention/coordinate basis before changing any physical parameter.
