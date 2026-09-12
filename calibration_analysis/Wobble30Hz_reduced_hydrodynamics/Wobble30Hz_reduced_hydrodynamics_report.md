# Reduced-Hydro Wall-ON 8.333 ms report

## Decision
`E_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS`

## Runtime identity
One job completed at 8.333 ms with 83,330 direct increments (`dt=1e-7 s`) and 168 ODB frames (50 us). CEL fluid and CEL contact were removed; SmoothWall114 Robot-wall contact remains. External production Magpylib Socket: B0=10 mT, f=30 Hz, cone=30 deg, Bias=40 deg, phase=248 deg, sense=+1, G=6 mT=0.006 T, L=45 mm, FORWARD.

## Hydro
Body-following dissipative RP loads use provisional analytical coefficients cpar=4e-9, cperp=1.2e-8 N s/mm, kspin=1e-9, kwob=3e-9 N mm s. Force/torque powers are in `reduced_hydro_angular_dynamics.csv` and should be non-positive.

## Strict dynamics
Formal omega is finite-rotation log of consecutive orientation matrices; direct VR is independently extracted. Maximum absolute reconstruction difference is 528.927 rad/s. H is kg m2/s; Jmag is N mm s (multiply by 1e-3 for N m s); Wmag uses 1 N mm=1e-3 J. Resolved nonmagnetic terms are accounting quantities, not independent wall torque measurements.

## Windows

                 window  duration_ms   n  phase_rate_mean_Hz  R_phase_vs_30  omega_wobble_rms_rad_s  omega_wobble_peak_rad_s  Tmag_rms_Nmm  exact_gap_min_um  CPRESS_max_MPa  Wmag_end_J  Whydro_end_J  Wnonmag_end_J
           W0_preimpact       7.7001 155            4.507423       0.150247              370.032967               797.740821      0.009729         -1.151799        65.21447    0.000004 -6.694831e-09      -0.000004
             W1_contact       0.5000  11            9.215350       0.307178              104.563398               202.541998      0.010491         -1.151799        65.21447    0.000004 -6.714474e-09      -0.000004
W2_immediate_postimpact       0.1329   3           -4.181936       0.139398              213.770735               223.043667      0.010522         -0.873539         0.00000    0.000004 -6.732380e-09      -0.000004
        W3_intermediate       0.0000   1                 NaN            NaN                     NaN                      NaN           NaN               NaN             NaN         NaN           NaN            NaN
     W4_late_postimpact      -3.3671   0                 NaN            NaN                     NaN                      NaN           NaN               NaN             NaN         NaN           NaN            NaN

## Contact and translation
Exact gap is computed against authoritative Pipe_WALL_HELPER R3D4 triangles. Minimum signed gap is -1.152 um; aggregate General Contact CPRESS peak is 65.214 MPa. Because the CEL fluid was removed, this reduced run has no fluid-contact contribution; the value is from the remaining wall contact domain. Canonical Vt and displacement are in `reduced_hydro_translation.csv`.

## Interpretation
This job is a diagnostic reduced-hydrodynamics replacement for coarse CEL, not a calibrated final fluid model. Any difference from prior CEL is classified only after checking geometry/contact, field identity, direct-VR/finite-rotation consistency, and non-positive hydro power.

## Exactly one next step
Classification `E_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS`. If A, extend same configuration for stability; if B, do the approved zero-cost elliptical-polarization/contact-load-shaping audit; if C, extend same configuration to 16.667 ms; if D, revisit field amplitude only after output validation; if E, repair extraction.
