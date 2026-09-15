# Legacy wobble animations

These GIFs preserve two earlier Abaqus visualizations from the parent working folder. They are archived as visual records, not promoted as final validation results.

## Reverse-trajectory HighEndStop partial run

File: `CEL_CurrentCenterline_Cone15_Z110_PolarityMinus_Damp055_AfterBendReverse_HighEndStop_060_FilletR010_CEL_FSI_bend_validation.gif`

- Source: `J:\abaqusfangzhen\` (original root-level GIF).
- The accompanying `partial83` report identifies frame 83 at 49.8001 ms as the last complete field output. The GIF contains 84 images indexed 0 through 83; “83-frame” refers to the last result-frame index, not the image count.
- The 60 ms job did not complete normally. The partial result also records a mismatch between the AfterBend RP location and the Fortran/Socket magnetic-field reference point, so it is for visualization and trend inspection only, not final magnetic validation.

## Early prescribed centerline Y-wobble probe

File: `Job_RouteA_CEL_SOLID_centerline_ywobble_probe.gif`

- Source: `J:\abaqusfangzhen\` (original root-level GIF).
- The paired report specifies a 5 Hz physical wobble (time-compressed to 666.667 Hz), 0.15 mm global Y amplitude, and a 2 ms smooth ramp. The input prescribes the motion through Abaqus amplitude boundary conditions; this is the early programmed-motion case, without Magpylib-driven actuation.
- This short probe retains centerline tangent-following rotation and is not a final zero-penetration acceptance run; the solver log reports CEL deep-penetration diagnostic warnings.
