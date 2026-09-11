# COM-fixed calibration: forensic audit of superseded runs

Date: 2026-09-11  
Production job retained and still running: `WobbleCal_COMFixed_F40_Cone30_B10`

## Scope and data-retention statement

The first two calibration attempts used the same Abaqus job name. Each later
attempt opened its telemetry CSV with write mode and Abaqus replaced the ODB
family. Consequently, the raw telemetry/ODB of the superseded attempts is no
longer present as a separate file. This report therefore records only facts
that are recoverable from the launch logs, deck metadata, and the final
snapshots observed during monitoring; it does **not** invent cycle metrics for
data that was overwritten.

## Attempt 1 (superseded, stopped early)

- Requested calibration: F40, cone half-angle 30°, B=10 mT, gradient 3 mT.
- Initial phase passed to the server: 228°.
- Cone-axis bias passed to the server: 0°.
- Bridge: original `vuforc_socket_bridge.f`, which still added the old RP
  origin `(-7.510492583, -3.6676577112, -9.52100960582) mm`.
- The deck RP had already been relocated to the volume COM, so the bridge and
  deck differed by 0.052543 mm.
- The run was intentionally interrupted at about 0.9 ms when the phase
  mismatch with the production calibration metadata was found.

This attempt is **INVALID / NOT A CALIBRATION RESULT**. It is useful only as
an implementation smoke test (socket, compilation, and early integration).

## Attempt 2 (superseded, completed but invalid)

- Initial phase corrected to the production 248°.
- Cone-axis bias was still 0° (not production Bias40°).
- The original bridge still used the old RP origin, producing the same
  0.052543 mm pose-origin mismatch.
- The run reached 75 ms, 750 output frames and about 709k increments. The
  monitor showed a stable explicit increment around 1.06e-7 s and no solver
  crash.
- The final snapshot available before replacement was approximately
  `UR=(1.762, -1.590, 3.544) rad`; the magnetic field gate remained
  `|B|=10.000 mT`.

This attempt is **INVALID / NOT PHYSICS-INTERPRETABLE**. It cannot establish
40-Hz following for the requested production field because Bias0 is a different
field orientation and the bridge did not query the COM frame consistently.

## What can and cannot be compared

The only defensible comparison is configuration-level:

| Item | Attempt 1 | Attempt 2 | Correct final run |
|---|---:|---:|---:|
| phase | 228° | 248° | 248° |
| cone-axis bias | 0° | 0° | 40° |
| bridge RP origin | old | old | COM-fixed |
| duration | ~0.9 ms, stopped | 75 ms | 75 ms target, running |
| status | invalid | invalid | valid candidate |

No frequency ratio, phase-lock score, amplitude variation, or cycle-pass label
is assigned to Attempts 1–2. The final valid run must be analyzed only after
completion using its preserved COM telemetry and ODB.

## Root causes and corrective actions

1. The calibration deck moved RP_ROBOT to the volume-weighted mesh COM, but
   the legacy bridge added a hard-coded old origin. A dedicated
   `vuforc_socket_bridge_comfixed.f` now uses the COM coordinates.
2. The server launch script initially passed phase 228° and then Bias0°.
   The final launch now passes production phase 248° and Bias40°.
3. All future calibration decks should carry a unique run identifier or write
   telemetry to a run-specific immutable path before launch, so an invalid
   rerun cannot overwrite an earlier diagnostic record.

## Decision

Attempts 1–2 are retained as invalid implementation diagnostics only. Do not
use them to infer whether the robot follows 40 Hz. The currently running
COM-fixed / Bias40 / phase248 job is the first configuration eligible for the
formal `MAGNETIC_WOBBLE_CALIBRATION_PASS` test.
