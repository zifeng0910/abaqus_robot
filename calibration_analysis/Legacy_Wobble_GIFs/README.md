# Legacy wobble animations

These GIFs preserve two earlier Abaqus visualizations from the parent working folder. They are archived as visual records, not promoted as final validation results.

## Reverse-trajectory HighEndStop partial run

File: `CEL_CurrentCenterline_Cone15_Z110_PolarityMinus_Damp055_AfterBendReverse_HighEndStop_060_FilletR010_CEL_FSI_bend_validation.gif`

- Source: `J:\abaqusfangzhen\` (original root-level GIF).
- The accompanying `partial83` report identifies frame 83 at 49.8001 ms as the last complete field output. The GIF contains 84 images indexed 0 through 83; “83-frame” refers to the last result-frame index, not the image count.
- The 60 ms job did not complete normally. The partial result also records a mismatch between the AfterBend RP location and the Fortran/Socket magnetic-field reference point, so it is for visualization and trend inspection only, not final magnetic validation.

## Prescribed head-tail rocking reference baseline

File: `Job_RouteA_CEL_SOLID_headtail_rock_probe.gif`

- Source: `J:\abaqusfangzhen\` (original root-level GIF).
- `CEL_REFERENCE_BASELINE.md` identifies this as the no-live-Magpylib head-tail alternating-rocking reference: ±10° about the transported local binormal, with no prescribed global Y/Z translation.
- The paired report and generator specify a 5 Hz physical reference, but the 0.03 s Abaqus probe compresses it to 666.667 Hz to display 20 cycles. The 121-frame GIF plays for about 7.26 s, so its animation is smoother/slower to view than the later 21-frame centerline Y-wobble probe; it is not a real-time 5 Hz playback.
- This is a kinematic/contact-shape baseline rather than magnetic-force validation or a final zero-penetration acceptance run; the solver reports CEL deep-penetration diagnostic warnings.
