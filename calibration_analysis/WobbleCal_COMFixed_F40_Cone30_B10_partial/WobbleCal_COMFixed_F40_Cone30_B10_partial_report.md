# WobbleCal_COMFixed_F40_Cone30_B10 — partial calibration analysis

## Status

The requested duration was **75.0 ms**.  The ODB contains frames through
**73.600 ms** and the merged Socket telemetry through
**73.614 ms** (98.15% of the request).
The lock file remains present but no `explicit_dp.exe` process is visible and
no final `ANALYSIS HAS COMPLETED SUCCESSFULLY` record was written.  This is
therefore classified as **`PARTIAL_STOP_OR_STALE_LOCK`**, not as a
successful complete run.

## Drive continuity

- Linear fit to the transmitted analytic phase: **40.000 Hz**
  (configured 40 Hz).
- Median transmitted field magnitude: **0.01 T**.
- Telemetry was recovered after the original server stopped; the main and
  recovery CSV streams were merged by time and deduplicated.

## Rotation response in the available window

- UR1 range: **-2.5464 to 1.3046 rad**
- UR2 range: **-2.4737 to 4.8150 rad**
- UR3 range: **-3.2638 to 2.6614 rad**
- Peak three-axis rotation-vector norm: **4.9542 rad**
- Applied torque norm: median **0.0020504 N mm**, peak
  **0.0111753 N mm**.

The COM-fixed diagnostic intentionally suppresses translation, so this run
can only calibrate the rotational response and field continuity.  It cannot
validate forward transport, wall interaction, or a 75 ms cycle-to-cycle
steady state because the final 1.4 ms is absent and the job did not close
normally.

## Files

- `runtime_merged.csv` — merged main/recovery telemetry and ODB UR history
- `cycle_metrics.csv` — available cycle-window metrics
- `summary.json` — machine-readable summary
- `runtime_diagnostic.png` — UR, torque, and analytic phase plots

## Decision

**Do not use this partial run as a final calibration baseline.**  Preserve the
ODB and lock file for forensic inspection; after confirming the stale lock and
server shutdown cause, rerun only if a complete 75 ms record is required.
