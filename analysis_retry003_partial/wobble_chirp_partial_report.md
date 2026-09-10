# Retry003 partial Chirp analysis

Status: **PARTIAL_STOPPED_AT_34.2_MS**. The 40 ms run did not complete; this report uses only the preserved telemetry through **34.213 ms**.

## Data integrity

- Telemetry rows: `313,153`; median sampling interval: `0.111 µs` (p01–p99 `0.105–0.112 µs`).
- Commanded frequency: `5.000 → 34.936 Hz`; integrated phase change: `245.94°`.
- Field magnitude: `10.000000–10.000000 mT`; bad 10 mT rows: `0`.

## Motion

- Driver arc: `18.89996 → 19.10524 mm`.
- Robot arc: `13.89996 → 13.99996 mm`; Δs (driver−robot): `5.00000 → 5.10528 mm`.
- Projected robot-arc net speed: `2.923 mm/s`; 3-D RP displacement: `8.993 mm`; peak RP speed: `334.080 mm/s`.
- Tangent force mean/median: `-0.00410/-0.00423 mN`; positive fraction: `4.27%`.

## Rotation

- UR1 range: `1.0691 rad` (min `-0.9039`, max `0.1652`); `||UR||` peak `2.8955 rad`.
- Applied torque magnitude peak: `0.0115968 N·mm`.
- Welch descriptive top frequencies: `136.912;273.825;410.737;547.650;684.562;821.474;958.387;1095.299 Hz`. These are not a step-out/locking result because the record is only 34.2 ms and the drive is chirped.

## Interpretation

1. The constant-field gate passed exactly (`10 mT` throughout the preserved rows), so the previous 7.5 mT field-amplitude bug is not present in Retry003.
2. The driver advances only about `0.205 mm`, while the robot's centerline projection advances about `0.100 mm`; the projected motion is therefore not a sustained 6 mm/s forward-following trajectory in this partial window.
3. `Ft` is positive for only `4.3%` of samples and has negative mean, so this Wall-OFF diagnostic does **not** yet demonstrate a forward magnetic-drive window. The large 3-D RP motion is mostly off-centerline, not confirmed forward transport.
4. Because the run stopped without an Abaqus error message and left a lock file, treat this as an interrupted partial result. Do not infer a completed Chirp lock band or submit a Wall-ON candidate from it.

## Files

- `wobble_chirp_partial_summary.csv`
- `wobble_chirp_partial_2ms_bins.csv`
- `wobble_chirp_partial_timeseries.csv`
- `wobble_chirp_partial_overview.png`
- `wobble_chirp_partial_rotation.png`
- `wobble_chirp_partial_spectrogram.png`
