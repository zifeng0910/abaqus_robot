# F100 G2P20 CEL Boundary 12 h Campaign Report

## Primary classification

`CEL_BOUNDARY_TREATMENT_STILL_UNRESOLVED`

## Scope

This bounded campaign repaired the short 12 mm fluid-only axial reflection harness by repeating the authorized H0/H1 comparison at 24 mm. It did not run a robot model, free dynamics, restart continuation, parameter tuning, or a broad sweep.

## Runs and verification

- H0 `HARNESS_DEFAULT_FREE`: completed successfully at 0.25 ms.
- H1 `HARNESS_NONREFLECTING`: completed successfully at 0.25 ms.
- Mesh: 88 x 20 x 20 EC3D8R; axial spacing 0.272727 mm; transverse spacing 0.075 mm.
- Stable increment: approximately 4.46e-7 s for H0 and 3.73e-7 s during H1; no fatal solver error.
- Expected 24 mm center-to-end travel time: 120 us; expected round trip: 240 us.

## 24 mm result

| Case | Probe | R_pressure | R_velocity |
|---|---|---:|---:|
| H0 free | low quarter | 0.344 | 0.463 |
| H0 free | high quarter | 0.352 | 0.540 |
| H1 nonreflecting | low quarter | 0.489 | 0.578 |
| H1 nonreflecting | high quarter | 0.506 | 0.590 |

H1 returns remain substantial and are not reduced by approximately 50%; the selected H1 ratios are higher than H0. Arrival times are in the expected 120-240 us scale, so the result is not explained by the 12 mm short-domain timing ambiguity. The pressure quantity is the explicitly labelled mean-normal-stress proxy `-trace(S_water)/3`; native `PRESS` and face mass-flux histories were unavailable.

## Decision

Do not run H2 `ZERO PRESSURE`, because the branch condition requiring clear H1 suppression was not met. Do not run another harness length, a full robot replay, L36/L48/L72/L100, wall-fluid-off replay, restart continuation, or tuning in this campaign. The boundary treatment remains unresolved and any inferred robot-CEL residual remains distinct from native fluid traction.

## Artifacts

Lightweight 24 mm CSV, JSON, and PNG outputs are under `boundary_harness_results_24mm/`. Abaqus ODB, SIM, restart, and native large files are intentionally excluded from scoped deliverables.
