# Centerline frame audit for COM-fixed wobble GIF

The original `curvenew_CEL_xyrot56_exact.csv` is in the Magpylib flat-DXF
frame.  The ODB pipe, robot, and fluid coordinates are in the Abaqus global
frame.  The first GIF passed the Magpylib CSV directly to the renderer, so its
black centerline was visibly offset/tilted even though the solver geometry was
unchanged.  This was a **post-processing frame-mixing error**, not evidence
that the 40 Hz Abaqus run integrated a distorted centerline.

The production server uses the shared rigid transform on every request:

```text
p_mag = R_aba_to_mag @ (p_aba - origin_aba_mm)
v_aba = R_aba_to_mag.T @ v_mag
```

The corrected plotting centerline is generated with:

```text
p_aba = R_aba_to_mag.T @ p_mag + origin_aba_mm
```

using `transform_centerline_to_abaqus.py` and
`curvenew_CEL_xyrot56_exact_abaqus.csv`.  Re-extraction and re-rendering with
that file makes the centerline follow the ODB pipe wall.  No Abaqus rerun is
required for this visualization correction.
