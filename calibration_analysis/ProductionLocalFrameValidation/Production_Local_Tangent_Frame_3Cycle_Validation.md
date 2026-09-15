# Production local-tangent frame validation

## Classification

`PRODUCTION_LOCAL_FRAME_REGRESSION_FAILED`

The candidate correctly implements the requested continuous robot-local frame,
but it does not reproduce the previously executed `FRAME_LOCAL30` screening
wrapper at the required component-wise tolerance on the moving trajectory.
Accordingly, the production file was not deployed and the 100 ms Abaqus case was
not run.

## Gate results

| Gate | Result | Maximum discrepancy |
|---|---:|---:|
| Legacy B, F, T; sense = -1 and +1 | PASS | 0 |
| Fixed-COM independent B equation | PASS | 3.47e-18 T |
| Fixed-COM field magnitude | PASS | 3.47e-18 T |
| Fixed-COM directed cone angle | PASS | 2.13e-14 deg |
| Fixed-COM winding | PASS | 1.000000 turn |
| Segment interior projection | PASS | 2.22e-16 mm |
| Endpoint extrapolation | PASS | -0.25 / +0.25 mm |
| Moving `FRAME_LOCAL30` B component | **FAIL** | 3.68e-4 T |
| Moving `FRAME_LOCAL30` B direction | **FAIL** | 2.825 deg |
| Moving `FRAME_LOCAL30` B magnitude | PASS | 5.20e-18 T |
| Moving `FRAME_LOCAL30` force, G = 0 | PASS | 0 N |
| Moving `FRAME_LOCAL30` torque component | **FAIL** | 4.19e-4 N mm |

The moving gate sampled 361 uniformly spaced states from the existing dense
`FRAME_LOCAL30` RP displacement and rotation histories. It did not run Abaqus.

## Root cause of the regression failure

The completed screening wrapper calls the historical `_project_arc_mm`, which
chooses the nearest resampled DXF vertex. The requested production correction
projects onto curve segments and evaluates a continuously transported frame at
the resulting arc coordinate. On a curved polyline these are different field
definitions: the screening frame changes at nearest-vertex boundaries, while
the candidate frame changes continuously along each segment.

As a controlled check, only the candidate's segment projector was replaced by
the historical nearest-vertex projector. With every other candidate code path
unchanged, its B, F, and torque matched the screening wrapper exactly (all
reported maximum differences were 0). This isolates projection semantics as the
cause and excludes UR mapping, phase, rotation sense, coordinate transform, and
magnetic moment direction.

## Stop decision

The specification requires both continuous segment projection and numerical
identity with the already executed nearest-vertex screening implementation.
Those requirements cannot both hold on the curved moving trajectory. The
explicit instruction for a moving B/F/T mismatch is to stop. Therefore:

- `J:\magpy\magpylib_socket_server.py` remains unchanged.
- `PROD_LOCAL30_G0_3CYCLE` was not submitted.
- No ODB, SIM, raw telemetry, private NPZ, DLL, cache, or TIFF was added.
- No dynamic classification (`LOCAL_FRAME_STABLE_WALL_WOBBLE`, transient,
  sliding, or tumble) is claimed.

## Reproduction

Run the offline gate with the configured Abaqus Python environment:

```powershell
& 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe' `
  'calibration_analysis\ProductionLocalFrameValidation\scripts\regress_production_local_frame.py'
```

The command intentionally exits nonzero while the moving screening-identity
gate fails and writes `field_regression_metrics.json` before exiting.
